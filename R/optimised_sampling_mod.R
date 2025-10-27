#' @title Boosted Sampling Optimization (Modernized)
#'
#' @description
#' Performs multi-objective sampling optimization by generating multiple 
#' permutations of random sampling points and evaluating each configuration
#' based on weighted scores for:
#' \itemize{
#'  \item NDVI variance (environmental heterogeneity)
#'  \item Mean ignorance value (prioritizing understudied areas)
#'  \item Spatial regularity/dispersion (using Ripley's L-function penalty)
#' }
#' Returns the best non-overlapping configuration.
#'
#' @param ndvi A `terra SpatRaster` object representing NDVI or other 
#'  environmental index.
#' @param ignorance A `terra SpatRaster` object representing the Map of 
#'  Relative Floristic Ignorance (MRFI).
#' @param boundary An `sf` polygon object defining the study area boundary.
#' @param samp_strategy Character. Sampling strategy for `sf::st_sample()`.
#'  Options: "random", "regular", "hexagonal". Default = "random".
#' @param nplot Integer. Number of sampling plots per permutation.
#' @param areaplot Numeric. Area of each circular plot in square meters,
#'  used to calculate buffer radius for overlap detection.
#' @param perm Integer. Number of permutations (sampling configurations) to test.
#' @param ndvi.weight Numeric. Weight for NDVI variance score (default = 1).
#' @param igno.weight Numeric. Weight for mean ignorance score (default = 1).
#' @param dist.weight Numeric. Weight for spatial dispersion score (default = 1).
#' @param verbose Logical. Print progress messages (default = TRUE).
#' @param output_dir Character. Directory for saving output files.
#' @param output_prefix Character. Prefix for output filenames (default = "SampleBoost").
#'
#' @return A list containing the run results, including the best solution points and plots.
#' 
#' @details
#' The spatial score uses Ripley's L-function to penalize clustering. The score
#' is the **maximum observed deviation above Complete Spatial Randomness (CSR)**,
#' meaning a lower score indicates better dispersion. This score is then 
#' inverted (multiplied by -1) before normalization so that lower clustering 
#' results in a higher final score.
#' 
#' @importFrom sf st_transform st_crs st_sample st_as_sf st_buffer st_intersects st_coordinates st_geometry st_bbox
#' @importFrom terra extract plot crs project resample
#' @importFrom tidyterra geom_spatraster
#' @importFrom dplyr bind_rows sample_n
#' @importFrom ggplot2 ggplot geom_sf scale_fill_distiller ggtitle theme_minimal aes geom_density theme_classic labs coord_sf geom_point
#' @importFrom plot3D scatter3D
#' @importFrom spatstat.geom as.owin ppp 
#' @importFrom spatstat.explore Lest
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom stats aggregate na.omit var
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.draw
#' @importFrom gridExtra grid.arrange tableGrob
#' @export
#'
sampleboost_mod <- function(ndvi, 
                            ignorance, 
                            boundary, 
                            samp_strategy = "random", 
                            nplot, 
                            areaplot, 
                            perm, 
                            ndvi.weight = 1, 
                            igno.weight = 1, 
                            dist.weight = 1,
                            verbose = TRUE,
                            output_dir = getwd(),
                            output_prefix = "SampleBoost") {
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  ## ---- 1. Input Validation ----
  msg("Validating inputs...")
  
  if (missing(ndvi) || missing(ignorance) || missing(boundary) || missing(nplot) || missing(areaplot) || missing(perm)) {
    stop("Missing required arguments (ndvi, ignorance, boundary, nplot, areaplot, perm).")
  }
  if (nplot < 2) stop("'nplot' must be at least 2")
  if (areaplot <= 0) stop("'areaplot' must be positive")
  if (perm < 1) stop("'perm' must be at least 1")
  if (ndvi.weight < 0 || igno.weight < 0 || dist.weight < 0) stop("All weights must be non-negative")
  if (ndvi.weight == 0 && igno.weight == 0 && dist.weight == 0) stop("At least one weight must be greater than zero")
  valid_strategies <- c("random", "regular", "hexagonal", "Fibonacci")
  if (!samp_strategy %in% valid_strategies) stop("'samp_strategy' must be one of: ", paste(valid_strategies, collapse = ", "))
  
  ## ---- 2. CRS & Extent Alignment ----
  msg("Checking and aligning coordinate reference systems...")
  
  # 2.1 Raster Alignment
  if (!identical(terra::crs(ndvi), terra::crs(ignorance))) {
    msg("Rasters have different CRS. Reprojecting 'ignorance' to match NDVI CRS...")
    ignorance <- terra::project(ignorance, ndvi)
  }
  if (!all(terra::res(ndvi) == terra::res(ignorance))) {
    msg("Rasters have different resolutions. Resampling 'ignorance' to NDVI resolution...")
    ignorance <- terra::resample(ignorance, ndvi, method = "bilinear")
  }
  
  # 2.2 Boundary Alignment
  if (inherits(boundary, "Spatial")) boundary <- sf::st_as_sf(boundary)
  if (is.na(sf::st_crs(boundary))) {
    msg("Boundary has no CRS; assuming EPSG:4326")
    sf::st_crs(boundary) <- 4326
  }
  boundary_proj <- sf::st_transform(boundary, terra::crs(ndvi))
  
  # 2.3 Spatial Metrics Setup
  boundary_owin <- spatstat.geom::as.owin(boundary_proj)
  plot_radius <- sqrt(areaplot / pi)
  msg(paste0("Plot radius: ", round(plot_radius, 2), " meters"))
  
  ## ---- 3. Normalization Function ----
  normalize <- function(x) {
    x_clean <- x[!is.na(x)]
    if(length(x_clean) == 0) return(x)
    rng <- max(x_clean) - min(x_clean)
    if(rng == 0) return(rep(0.5, length(x)))
    (x - min(x_clean)) / rng
  }
  
  ## ---- 4. Generate and Evaluate Permutations (L-FUNCTION) ----
  msg(paste0("Generating and evaluating ", perm, " permutations..."))
  
  permutation_results <- vector("list", perm)
  clustering_penalty <- numeric(perm) # Stores max(L(r) - r)
  has_intersection <- logical(perm)
  
  pb <- utils::txtProgressBar(min = 0, max = perm, style = 3)
  
  for (i in 1:perm) {
    sample_points <- sf::st_sample(boundary_proj, size = nplot, type = samp_strategy, iter = 10)
    points_sf <- sf::st_as_sf(sample_points)
    
    point_buffers <- sf::st_buffer(points_sf, dist = plot_radius)
    intersection_matrix <- sf::st_intersects(point_buffers, sparse = FALSE)
    diag(intersection_matrix) <- FALSE
    has_intersection[i] <- any(intersection_matrix)
    
    ndvi_values <- terra::extract(ndvi, points_sf, ID = FALSE)[,1]
    ignorance_values <- terra::extract(ignorance, points_sf, ID = FALSE)[,1]
    
    coords <- sf::st_coordinates(points_sf)
    
    # Calculate L-function for spatial dispersion metric
    pp <- spatstat.geom::ppp(coords[,1], coords[,2], window = boundary_owin)
    L_function <- spatstat.explore::Lest(pp, correction = "isotropic")
    
    # Spatial Score: Maximum deviation above CSR (L(r) - r > 0)
    # Higher values mean worse clustering/less spreading (PENALTY).
    clustering_penalty[i] <- max(L_function$iso - L_function$r) 
    
    permutation_results[[i]] <- data.frame(
      x = coords[,1],
      y = coords[,2],
      ndvi = ndvi_values,
      ignorance = ignorance_values,
      has_intersection = has_intersection[i],
      permutation_id = i,
      stringsAsFactors = FALSE
    )
    
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ## ---- 5. Aggregate Results ----
  msg("Aggregating and scoring permutations...")
  full_matrix <- dplyr::bind_rows(permutation_results)
  full_matrix$permutation_id <- as.factor(full_matrix$permutation_id)
  
  ndvi_variance <- stats::aggregate(full_matrix$ndvi, by = list(permutation_id = full_matrix$permutation_id), FUN = stats::var)
  mean_ignorance <- stats::aggregate(full_matrix$ignorance, by = list(permutation_id = full_matrix$permutation_id), FUN = mean)
  
  aggregated_scores <- data.frame(
    permutation_id = 1:perm,
    ndvi_variance = ndvi_variance$x,
    # Clustering penalty is stored here
    clustering_penalty = clustering_penalty, 
    mean_ignorance = mean_ignorance$x,
    has_intersection = has_intersection,
    stringsAsFactors = FALSE
  )
  
  aggregated_scores <- stats::na.omit(aggregated_scores)
  
  # 5.1 Apply Weights and Normalize
  
  # Maximize objectives: NDVI Variance and Mean Ignorance
  aggregated_scores$ndvi_weighted <- aggregated_scores$ndvi_variance * ndvi.weight
  aggregated_scores$igno_weighted <- aggregated_scores$mean_ignorance * igno.weight
  
  # Minimize objective: Clustering Penalty (L-function)
  # INVERSION: Multiply by -1 so that low clustering (low penalty) becomes a high value
  aggregated_scores$spatial_weighted <- aggregated_scores$clustering_penalty * -1 * dist.weight
  
  # Normalize (all objectives are now "higher is better")
  aggregated_scores$ndvi_norm <- normalize(aggregated_scores$ndvi_weighted)
  aggregated_scores$igno_norm <- normalize(aggregated_scores$igno_weighted)
  aggregated_scores$spatial_norm <- normalize(aggregated_scores$spatial_weighted) 
  
  # Calculate final combined score
  aggregated_scores$final_score <- aggregated_scores$ndvi_norm + aggregated_scores$igno_norm + aggregated_scores$spatial_norm
  
  ## ---- 6. Select Best Non-Intersecting Solution ----
  msg("Selecting best non-overlapping configuration...")
  valid_solutions <- aggregated_scores[!aggregated_scores$has_intersection, ]
  
  if (nrow(valid_solutions) == 0) {
    warning("No valid solutions found. Try increasing 'perm', decreasing 'areaplot', or decreasing 'nplot'.")
    return(list(full_matrix = full_matrix, aggregated_scores = aggregated_scores, best_solution = NULL, best_solution_sf = NULL))
  }
  
  best <- valid_solutions[order(valid_solutions$final_score, decreasing = TRUE),][1,]
  best_solution <- full_matrix[full_matrix$permutation_id == best$permutation_id, ]
  best_solution_sf <- sf::st_as_sf(best_solution, coords = c("x","y"), crs = sf::st_crs(boundary_proj))
  best_buffers_sf <- sf::st_buffer(best_solution_sf, dist = plot_radius)
  
  ## ---- 7. Statistics ----
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  statistics <- data.frame(
    Statistic = c("Started", "Finished", "Elapsed_time_secs", "CRS_EPSG", "Sampling_strategy",
                  "Number_of_plots", "Plot_area_m2", "Plot_radius_m", "Permutations_tested",
                  "Valid_solutions", "NDVI_weight", "Ignorance_weight", "Distance_weight",
                  "Best_permutation_id", "Best_final_score", "Best_NDVI_variance", "Best_mean_ignorance",
                  "Best_L_Clustering_Penalty"),
    Value = c(as.character(start_time), as.character(end_time), round(elapsed_time, 2),
              as.character(sf::st_crs(boundary_proj)$epsg), samp_strategy, nplot, areaplot,
              round(plot_radius,2), perm, nrow(valid_solutions), ndvi.weight, igno.weight, dist.weight,
              best$permutation_id, round(best$final_score, 4), round(best$ndvi_variance,4),
              round(best$mean_ignorance,2), 
              round(best$clustering_penalty,4)),
    stringsAsFactors = FALSE
  )
  
  ## ---- 8. Plotting (Modernized & Corrected) ----
  msg("Generating plots and output files...")
  
  # Plot 1: Best solution on NDVI map
  plot_ndvi <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = ndvi) +
    ggplot2::geom_sf(data = boundary_proj, fill = NA, color="black", linewidth=1) +
    ggplot2::geom_sf(data = best_buffers_sf, fill = NA, color="red", linetype="dashed") +
    ggplot2::geom_sf(data = best_solution_sf, color = "darkred", size = 3, shape = 21, fill = "white", stroke = 1.5) +
    ggplot2::scale_fill_distiller(palette = "YlGn", name = "NDVI", na.value = "transparent", direction = 1) +
    ggplot2::ggtitle("Best Sampling Configuration on NDVI") +
    ggplot2::coord_sf() + ggplot2::theme_minimal()
  
  # Plot 2: Best solution on Ignorance map
  plot_ignorance <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = ignorance) +
    ggplot2::geom_sf(data = boundary_proj, fill = NA, color="black", linewidth=1) +
    ggplot2::geom_sf(data = best_buffers_sf, fill = NA, color="blue", linetype="dashed") +
    ggplot2::geom_sf(data = best_solution_sf, color = "darkblue", size = 3, shape = 21, fill = "white", stroke = 1.5) +
    ggplot2::scale_fill_distiller(palette = "RdBu", name = "Ignorance", na.value = "transparent", direction = -1) +
    ggplot2::ggtitle("Best Sampling Configuration on Ignorance Map") +
    ggplot2::coord_sf() + ggplot2::theme_minimal()
  
  # Plot 3: Scatter plot of permutation scores (normalized)
  plot_scores <- ggplot2::ggplot(valid_solutions, ggplot2::aes(x=ndvi_norm, y=igno_norm, size=spatial_norm, color=final_score)) +
    ggplot2::geom_point(alpha=0.7) +
    ggplot2::geom_point(data = best, color = "red", size = 6, shape = 18) +
    ggplot2::scale_color_distiller(palette = "Spectral", name = "Final Score") +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("Multi-Objective Optimization Space") +
    ggplot2::xlab("NDVI Variance (Normalized)") +
    ggplot2::ylab("Ignorance (Normalized)") +
    ggplot2::labs(size="Dispersion (Normalized)")
  
  # Plot 4: Histogram of L-function penalties
  plot_spatial <- ggplot2::ggplot(aggregated_scores, ggplot2::aes(x=clustering_penalty)) +
    ggplot2::geom_histogram(alpha=0.6, fill="#FF6666", bins=30) +
    ggplot2::geom_vline(xintercept = best$clustering_penalty, color="red", linetype="dashed", linewidth=1) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("Clustering Penalty Distribution (L-function)") +
    ggplot2::xlab("Max L(r) - r (Clustering Penalty, lower is better)") +
    ggplot2::ylab("Frequency") +
    ggplot2::labs(subtitle="Red line indicates the best solution's penalty score.")
  
  # Save Outputs
  csv_best <- file.path(output_dir, paste0(output_prefix, "_best-solution.csv"))
  csv_all <- file.path(output_dir, paste0(output_prefix, "_all-scores.csv"))
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_plots.pdf"))
  
  utils::write.csv(best_solution, csv_best, row.names = FALSE)
  utils::write.csv(aggregated_scores, csv_all, row.names = FALSE)
  
  # Save PDF
  grDevices::pdf(pdf_path, width = 10, height = 8, onefile = TRUE)
  print(plot_ndvi)
  print(plot_ignorance)
  print(plot_scores)
  print(plot_spatial)
  
  # Note: 3D plot is not automatically saved in the PDF, it requires manual handling
  # The statistics table is saved:
  grid::grid.draw(gridExtra::grid.arrange(top="Summary Statistics", gridExtra::tableGrob(statistics)))
  
  grDevices::dev.off()
  
  msg(paste0("Done! Files saved to: ", output_dir))
  
  ## ---- 9. Return ----
  return(list(
    full_matrix = full_matrix,
    aggregated_scores = aggregated_scores,
    best_solution = best_solution,
    best_solution_sf = best_solution_sf,
    best_variance = best$ndvi_variance,
    best_ignorance = best$mean_ignorance,
    best_clustering_penalty = best$clustering_penalty,
    best_final_score = best$final_score,
    plots = list(ndvi = plot_ndvi, ignorance = plot_ignorance, scores = plot_scores, spatial = plot_spatial),
    statistics = statistics
  ))
}
