#' @title Virtual Floristic List (Modernized & Fixed)
#'
#' @description
#' Computes a Virtual Floristic List (VFL): taxa potentially occurring within
#' a study site, with probabilities based on spatial uncertainty and temporal
#' decay using inclusion-exclusion principle.
#'
#' @param data_flor Data frame with columns: 'Taxon', 'Long', 'Lat',
#'   'uncertainty' (radius in meters), 'year'.
#' @param site sf polygon or SpatialPolygonsDataFrame of study area.
#' @param year_study Numeric year of analysis (default = current year).
#' @param excl_areas Optional sf polygon(s) of unsuitable areas to exclude.
#' @param CRS.new Numeric EPSG code for projected CRS (default = 3035).
#' @param tau Percent taxa loss per 100 years (0 ≤ tau < 100).
#' @param upperlimit Max records per taxon in inclusion-exclusion (default=20).
#' @param verbose Logical; print progress messages (default = TRUE).
#' @param output_dir Directory for output files (default = working directory).
#' @param output_prefix Filename prefix (default = "VFL").
#' @param use_parallel Use parallel processing for large datasets (default=TRUE).
#' @param n_cores Number of cores (default = half available cores).
#'
#' @return List with:
#' \itemize{
#'   \item{\code{VFL}}{Data frame: Taxon, probability, records, max, min.}
#'   \item{\code{Statistics}}{Metadata table.}
#'   \item{\code{Plots}}{Named list of ggplot objects.}
#'   \item{\code{spatial_data}}{sf objects for further analysis (optional).}
#' }
#'
#' @importFrom sf st_as_sf st_make_valid st_crs st_transform st_buffer st_intersection st_union st_area st_geometry st_difference st_coordinates
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point geom_sf theme_classic xlab ylab ggtitle coord_sf
#' @importFrom utils write.csv txtProgressBar setTxtProgressBar
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.draw
#' @importFrom gridExtra grid.arrange tableGrob
#' @export
#'
#' @examples
#' \dontrun{
#' # Using sf objects
#' data(floratus)
#' data(park)
#' park_sf <- sf::st_as_sf(park)
#'
#' # Short example
#' set.seed(123)
#' vfl <- virtual_list_mod(
#'   data_flor = floratus[sample(nrow(floratus), 2000), ],
#'   site = park_sf,
#'   tau = 30,
#'   upperlimit = 25
#' )
#'
#' # View results
#' head(vfl$VFL)
#' vfl$Statistics
#' print(vfl$Plots$prob_histogram)
#' }

virtual_list_mod <- function(data_flor, site, year_study = NULL,
                             excl_areas = NULL, CRS.new = 3035,
                             tau, upperlimit = 20, verbose = TRUE,
                             output_dir = getwd(),
                             output_prefix = "VFL",
                             use_parallel = TRUE,
                             n_cores = NULL) {
  
  # Helper for conditional messages
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  ## ---- 1. Input Validation ----
  msg("Validating inputs...")
  
  # Set year_study to current year if not provided
  if (is.null(year_study)) {
    year_study <- as.numeric(format(Sys.Date(), "%Y"))
  }
  
  # Validate tau parameter
  if (tau < 0 || tau >= 100) {
    stop("tau must satisfy: 0 <= tau < 100")
  }
  
  # Check required columns
  req_cols <- c("Taxon", "Long", "Lat", "uncertainty", "year")
  missing_cols <- setdiff(req_cols, names(data_flor))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for future dates
  if (max(data_flor$year, na.rm = TRUE) > year_study) {
    warning("Some occurrence dates are more recent than year_study")
  }
  
  # Validate upperlimit
  if (upperlimit < 1) {
    stop("upperlimit must be at least 1")
  }
  if (upperlimit > 50) {
    warning("Large upperlimit (>50) may cause memory issues and slow computation")
  }
  
  # Store initial counts for statistics
  total_initial_records <- nrow(data_flor)
  total_initial_taxa <- length(unique(data_flor$Taxon))
  
  ## ---- 2. CRS Setup and Reprojection ----
  msg(paste("Reprojecting to EPSG:", CRS.new))
  
  crs_sf <- sf::st_crs(CRS.new)
  crs_terra <- paste0("EPSG:", CRS.new)
  
  # Handle site polygon
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  if (is.na(sf::st_crs(site))) {
    msg("Site has no CRS; assuming EPSG:4326")
    sf::st_crs(site) <- 4326
  }
  site_proj <- sf::st_transform(sf::st_make_valid(site), crs_sf)
  
  # Handle exclusion areas if provided
  excl_proj <- NULL
  if (!is.null(excl_areas)) {
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    if (is.na(sf::st_crs(excl_areas))) {
      msg("Exclusion areas have no CRS; assuming EPSG:4326")
      sf::st_crs(excl_areas) <- 4326
    }
    excl_proj <- sf::st_union(
      sf::st_transform(sf::st_make_valid(excl_areas), crs_sf)
    )
    msg("Exclusion areas provided and will be removed from buffers")
  }
  
  # Handle floristic data
  if (inherits(data_flor, "sf")) {
    pts_sf <- data_flor
  } else if (inherits(data_flor, "Spatial")) {
    pts_sf <- sf::st_as_sf(data_flor)
  } else {
    pts_sf <- sf::st_as_sf(
      data_flor, 
      coords = c("Long", "Lat"), 
      crs = 4326, 
      remove = FALSE
    )
  }
  
  pts_proj <- sf::st_transform(pts_sf, crs_sf)
  
  ## ---- 3. Create Buffers ----
  msg("Creating uncertainty buffers...")
  
  # Add unique record ID before any operations
  pts_proj$record_id <- seq_len(nrow(pts_proj))
  
  # Create buffers with uncertainty radius
  buffers_sf <- sf::st_buffer(pts_proj, dist = pts_proj$uncertainty)
  
  # Calculate buffer areas BEFORE intersection
  buffers_sf$area_buffer <- as.numeric(sf::st_area(buffers_sf))
  
  # Remove exclusion areas from buffers if provided
  if (!is.null(excl_proj)) {
    msg("Removing exclusion areas from buffers...")
    buffers_sf <- sf::st_make_valid(
      sf::st_difference(buffers_sf, excl_proj)
    )
    # Recalculate areas after exclusion
    buffers_sf$area_buffer <- as.numeric(sf::st_area(buffers_sf))
  }
  
  ## ---- 4. Intersection with Study Area ----
  msg("Computing buffer intersections with study area...")
  
  # Find which buffers intersect the site
  intersects_idx <- sf::st_intersects(buffers_sf, site_proj, sparse = FALSE)[, 1]
  
  if (sum(intersects_idx) == 0) {
    stop("No occurrence buffers intersect the study area. Check CRS and coordinates.")
  }
  
  # Filter to only intersecting buffers
  buffers_intersecting <- buffers_sf[intersects_idx, ]
  
  msg(paste("Retained", nrow(buffers_intersecting), "buffers that intersect site"))
  
  # Perform intersection (this preserves area_buffer from buffers_sf)
  inter_sf <- sf::st_intersection(buffers_intersecting, site_proj)
  
  # Calculate intersection areas
  inter_sf$area_intersection <- as.numeric(sf::st_area(inter_sf))
  
  ## ---- 5. Calculate Probabilities ----
  msg("Calculating spatial and temporal probabilities...")
  
  # Spatial probability: intersection area / buffer area
  inter_sf$p_spatial <- inter_sf$area_intersection / inter_sf$area_buffer
  
  # Temporal probability: exponential decay based on tau
  inter_sf$p_temporal <- (1 - (tau / 100))^((year_study - inter_sf$year) / 100)
  
  # Combined spatio-temporal probability
  inter_sf$p_st <- inter_sf$p_spatial * inter_sf$p_temporal
  
  # Validate probability ranges
  if (any(inter_sf$p_st > 1.0, na.rm = TRUE)) {
    warning("Some probabilities exceed 1.0; capping at 1.0")
    inter_sf$p_st <- pmin(inter_sf$p_st, 1.0)
  }
  if (any(inter_sf$p_st < 0, na.rm = TRUE)) {
    warning("Some probabilities are negative; setting to 0")
    inter_sf$p_st <- pmax(inter_sf$p_st, 0)
  }
  
  ## ---- 6. Aggregate by Taxon (Inclusion-Exclusion Principle) ----
  msg("Aggregating probabilities by taxon using inclusion-exclusion...")
  
  # Check if parallel processing is requested and is worth the overhead
  is_parallel_run <- FALSE
  if (isTRUE(use_parallel) && n_taxa > 10) {
    if (is.null(n_cores)) {
      # Calculate cores only if needed
      n_cores <- max(1, floor(parallel::detectCores() / 2))
    }
    
    # Attempt to set up parallel processing with a fallback
    tryCatch({
      future::plan(future::multisession, workers = n_cores)
      msg(paste("Using parallel processing with", n_cores, "cores"))
      is_parallel_run <- TRUE
    }, error = function(e) {
      msg(paste("Parallel setup failed (Error:", e$message, "). Falling back to sequential mode."))
      future::plan(future::sequential)
      # Set n_cores to 1 for statistics consistency if it was NULL
      if (is.null(n_cores)) n_cores <- 1
    })
  }
  
  # If parallel plan was not set or failed, force sequential
  if (!is_parallel_run) {
    future::plan(future::sequential)
    msg("Running in sequential mode.")
    # Set n_cores to 1 for statistics consistency if it was NULL
    if (is.null(n_cores)) n_cores <- 1
  }
  
  # Progress bar setup
  pb <- utils::txtProgressBar(min = 0, max = n_taxa, style = 3)
  
  # Process each taxon
  results <- future.apply::future_lapply(
    seq_along(taxa_list),
    function(i) {
      taxon <- taxa_list[i]
      # Extract probabilities for this taxon
      taxon_probs <- inter_sf$p_st[inter_sf$Taxon == taxon]
      # Compute aggregated probability
      result <- compute_taxon_prob(taxon_probs, upperlimit)
      # Update progress bar (only prints cleanly in sequential mode)
      if (!is_parallel_run) utils::setTxtProgressBar(pb, i)
      result
    },
    future.seed = TRUE
  )
  
  close(pb)
  future::plan(future::sequential) # ALWAYS reset the plan back to sequential after the loop
  
  ## ---- 7. Build VFL Data Frame ----
  msg("Building Virtual Floristic List...")
  
  VFL <- data.frame(
    Taxon = taxa_list,
    Estimated_Spatiotemporal_probability = sapply(results, `[`, 1),
    Number_of_records = sapply(results, `[`, 2),
    Max_probability = sapply(results, `[`, 3),
    Min_probability = sapply(results, `[`, 4),
    stringsAsFactors = FALSE
  )
  
  # Filter out taxa with zero probability
  VFL <- VFL[VFL$Estimated_Spatiotemporal_probability > 0, ]
  
  # Convert to percentages and round
  VFL$Estimated_Spatiotemporal_probability <- 
    round(VFL$Estimated_Spatiotemporal_probability * 100, 1)
  VFL$Max_probability <- round(VFL$Max_probability * 100, 1)
  VFL$Min_probability <- round(VFL$Min_probability * 100, 1)
  
  # Sort by probability (descending), then alphabetically
  VFL <- VFL[order(-VFL$Estimated_Spatiotemporal_probability, VFL$Taxon), ]
  rownames(VFL) <- NULL
  
  ## ---- 8. Compile Statistics ----
  msg("Compiling statistics...")
  end_time <- Sys.time()
  
  stats <- data.frame(
    Statistic = c(
      "Started",
      "Finished",
      "Elapsed_time_secs",
      "CRS_EPSG",
      "Year_study",
      "Tau_percent",
      "Upperlimit",
      "N_cores_used",
      "Total_initial_records",
      "Total_initial_taxa",
      "Records_intersecting_site",
      "Taxa_in_VFL",
      "Median_uncertainty_m",
      "Median_year"
    ),
    Value = c(
      as.character(start_time),
      as.character(end_time),
      round(as.numeric(difftime(end_time, start_time, units = "secs")), 2),
      as.character(CRS.new),
      year_study,
      tau,
      upperlimit,
      n_cores,
      total_initial_records,
      total_initial_taxa,
      nrow(inter_sf),
      nrow(VFL),
      round(median(inter_sf$uncertainty, na.rm = TRUE)),
      round(median(inter_sf$year, na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )
  
  ## ---- 9. Generate Diagnostic Plots ----
  msg("Generating diagnostic plots...")
  
  # Plot 1: Probability distribution histogram
  p_prob_hist <- ggplot2::ggplot(
    VFL, 
    ggplot2::aes(x = Estimated_Spatiotemporal_probability)
  ) +
    ggplot2::geom_histogram(
      fill = "#FF6666", 
      alpha = 0.7, 
      bins = 30,
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Spatio-temporal probability (%)") +
    ggplot2::ylab("Number of taxa") +
    ggplot2::ggtitle("Distribution of Occurrence Probabilities")
  
  # Plot 2: Probability vs number of records
  p_prob_vs_records <- ggplot2::ggplot(
    VFL,
    ggplot2::aes(
      x = Number_of_records,
      y = Estimated_Spatiotemporal_probability
    )
  ) +
    ggplot2::geom_point(alpha = 0.6, color = "darkblue", size = 2) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Number of records per taxon") +
    ggplot2::ylab("Probability (%)") +
    ggplot2::ggtitle("Probability vs Record Count")
  
  # Plot 3: Max vs min probability scatter
  p_max_vs_min <- ggplot2::ggplot(
    VFL,
    ggplot2::aes(x = Max_probability, y = Min_probability)
  ) +
    ggplot2::geom_point(alpha = 0.6, color = "darkred", size = 2) +
    ggplot2::geom_abline(
      intercept = 0, 
      slope = 1, 
      linetype = "dashed", 
      color = "gray50"
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Max probability (%)") +
    ggplot2::ylab("Min probability (%)") +
    ggplot2::ggtitle("Range of Individual Record Probabilities")
  
  # Plot 4: Spatial visualization (buffers and site)
  p_spatial <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = site_proj,
      fill = NA,
      color = "black",
      linewidth = 1
    ) +
    ggplot2::geom_sf(
      data = buffers_intersecting,
      fill = rgb(0, 0, 1, 0.1),
      color = rgb(0, 0, 1, 0.3),
      linewidth = 0.2
    ) +
    ggplot2::coord_sf() +
    ggplot2::theme_classic() +
    ggplot2::ggtitle("Uncertainty Buffers Intersecting Study Area") +
    ggplot2::xlab("Longitude") +
    ggplot2::ylab("Latitude")
  
  # Plot 5: Temporal distribution
  p_temporal <- ggplot2::ggplot(
    as.data.frame(inter_sf),
    ggplot2::aes(x = year)
  ) +
    ggplot2::geom_histogram(
      fill = "#66B2FF",
      alpha = 0.7,
      bins = 30,
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Year of occurrence") +
    ggplot2::ylab("Number of records") +
    ggplot2::ggtitle("Temporal Distribution of Records")
  
  # Combine plots in a named list
  plots_list <- list(
    prob_histogram = p_prob_hist,
    prob_vs_records = p_prob_vs_records,
    max_vs_min = p_max_vs_min,
    spatial_buffers = p_spatial,
    temporal_distribution = p_temporal
  )
  
  ## ---- 10. Save Outputs ----
  msg("Saving output files...")
  
  csv_path <- file.path(output_dir, paste0(output_prefix, "_VFL.csv"))
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_VFL.pdf"))
  
  # Save CSV
  utils::write.csv(VFL, csv_path, row.names = FALSE)
  
  # Save multi-page PDF
  grDevices::pdf(pdf_path, width = 10, height = 8, onefile = TRUE)
  print(p_prob_hist)
  print(p_prob_vs_records)
  print(p_max_vs_min)
  print(p_spatial)
  print(p_temporal)
  grid::grid.draw(
    gridExtra::grid.arrange(
      top = "Summary Statistics",
      gridExtra::tableGrob(stats)
    )
  )
  grDevices::dev.off()
  
  msg(paste0("Done! Files saved to: ", output_dir))
  msg(paste0("  - CSV: ", basename(csv_path)))
  msg(paste0("  - PDF: ", basename(pdf_path)))
  
  ## ---- 11. Display Plots ----
  if (verbose) {
    print(p_prob_hist)
    print(p_prob_vs_records)
  }
  
  ## ---- 12. Return Results ----
  return(list(
    VFL = VFL,
    Statistics = stats,
    Plots = plots_list,
    spatial_data = list(
      site = site_proj,
      buffers = buffers_intersecting,
      intersections = inter_sf
    )
  ))
}