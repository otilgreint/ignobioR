#' @title Boosted Sampling Optimization
#' @description
#' Performs a multi-objective sampling optimization by generating 'perm' 
#' permutations of 'nplot' random points. It evaluates each set based on
#' weighted scores for NDVI variance, map ignorance, and spatial distance.
#'
#' @param ndvi A `terra SpatRaster` object for NDVI (or other spectral index).
#' @param ignorance A `terra SpatRaster` object for the ignorance map.
#' @param boundary An `sf` object (polygon) defining the study area.
#' @param samp_strategy Character. The sampling strategy for `sf::st_sample` 
#' (e.g., 'random', 'regular').
#' @param nplot Integer. The number of plots (points) to generate per permutation.
#' @param areaplot Numeric. The area of each plot (in m^2) used to calculate
#' a buffer radius for intersection checks.
#' @param perm Integer. The number of permutations (sampling sets) to try.
#' @param ndvi.weight Numeric. The weight for the NDVI variance score.
#' @param igno.weight Numeric. The weight for the mean ignorance score.
#' @param dist.weight Numeric. The weight for the minimum distance score.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item{`Full matrix`}{A data.frame of all points from all permutations.}
#'   \item{`Aggregated matrix`}{A data.frame summarizing scores for each permutation.}
#'   \item{`Best`}{A data.frame of the points from the best-scoring permutation.}
#'   \item{`Variance of sampling points`}{NDVI variance for the best solution.}
#'   \item{`Mean Ignorance`}{Weighted ignorance score for the best solution.}
#'   \item{`Spatial Median of Distance`}{Minimum distance for the best solution.}
#'   \item{`Final score`}{The final combined score for the best solution.}
#'   \item{`Plot_NDVI`}{A ggplot object showing the best solution on the NDVI map.}
#'   \item{`Plot_Ignorance`}{A ggplot object showing the best solution on the ignorance map.}
#'   \item{`Plot_Density`}{A ggplot object showing the density of NDVI values.}
#' }
#'
#' @importFrom sf st_transform st_crs st_sample st_as_sf st_buffer st_intersects st_coordinates st_geometry
#' @importFrom terra plot extract
#' @importFrom tidyterra geom_spatraster
#' @importFrom dplyr bind_rows sample_n
#' @importFrom ggplot2 ggplot geom_sf scale_fill_distiller ggtitle theme_minimal aes geom_density theme
#' @importFrom plot3D scatter3D
#' @importFrom spatstat.geom as.owin ppp pairdist.ppp
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats aggregate na.omit var
#' @export

sampleboost <- function(ndvi, ignorance, boundary, samp_strategy, nplot, areaplot, perm, ndvi.weight, igno.weight, dist.weight){
  
  # Normalization function (made safer for divide-by-zero)
  normalize <- function(x) {
    x_na <- x[!is.na(x)]
    if (length(x_na) == 0) return(x)
    val_range <- max(x_na) - min(x_na)
    if (val_range == 0) return(rep(0.5, length(x))) # All values are the same
    return ((x - min(x_na)) / val_range)
  }
  
  start_time <- Sys.time() # Measure time: start
  result <- list()
  distanze <- matrix(ncol = 1, nrow = perm)
  check <- c()
  
  # Project boundary to match ignorance raster CRS using sf::st_transform()
  boundary <- sf::st_transform(boundary, sf::st_crs(ignorance))
  
  # Convert sf boundary to spatstat 'owin' for distance calculation window
  boundary_owin <- spatstat.geom::as.owin(boundary)
  
  pb <- utils::txtProgressBar(min = 0, max = perm, style = 3)
  for (i in 1:perm){
    
    # Use sf::st_sample() instead of spsample()
    punti_random <- sf::st_sample(boundary, size = nplot, type = samp_strategy, iter = 5)
    
    # Convert points to sf data.frame
    spdf <- sf::st_as_sf(punti_random)
    
    # Use sf::st_buffer() instead of gBuffer()
    spdf_buffer <- sf::st_buffer(spdf, dist = sqrt(areaplot / pi))
    
    # Test self intersection (much faster using sf matrix)
    int_matrix <- sf::st_intersects(spdf_buffer, sparse = FALSE)
    diag(int_matrix) <- FALSE # Ignore self-intersections
    check[[i]] <- any(int_matrix)
    
    # Use terra::extract() instead of raster::extract()
    # It works directly with sf objects and returns a data.frame
    spectral_values <- terra::extract(ndvi, spdf)[, 2] # Col 2 is the data
    igno_values <- terra::extract(ignorance, spdf)[, 2]
    
    ## Calculate distances with metric CRS
    xy <- sf::st_coordinates(spdf) # Get coordinates from sf object
    
    # 1. Get a ppp object from the data (using pre-calculated owin)
    m <- spatstat.geom::ppp(xy[, 1], xy[, 2], window = boundary_owin)
    # 2. Calculate a Euclidean distance matrix
    pairwise_distances <- spatstat.geom::pairdist.ppp(m)
    
    pairwise_distances[pairwise_distances == 0] <- NA # replace 0s on the diagonal with NA
    distanze[[i]] <- min(pairwise_distances, na.rm = TRUE) 
    distance_values <- rep(distanze[[i]], nplot)
    
    
    estratti <- data.frame(xy, spectral_values, igno_values, distance_values)
    names(estratti) <- c("x", "y", "ndvi", "ignorance", "distances")
    
    estratti$INTERSECTION <- check[[i]]
    result[[i]] <- data.frame(estratti)
    
    utils::setTxtProgressBar(pb, i)
    
  }
  
  # Use dplyr::bind_rows() instead of plyr::ldply()
  new_mat <- dplyr::bind_rows(result)
  new_mat$try <- as.factor(rep(1:perm, each = nplot))
  
  agg1 <- stats::aggregate(new_mat$ndvi, by = list(new_mat$try), FUN = stats::var)
  agg_igno <- stats::aggregate(new_mat$ignorance, by = list(new_mat$try), FUN = mean)
  
  
  agg2 <- data.frame(agg1, distanze, agg_igno[[2]], unlist(check))
  colnames(agg2) <- c('Try', 'Variance', 'Mean Dist', 'Mean Ignorance', "INTERSECTION")
  agg2 <- stats::na.omit(agg2)
  
  agg2$ndvi_score <- agg2$Variance * ndvi.weight # apply the weighted value
  agg2$igno_score <- agg2$`Mean Ignorance` * igno.weight
  agg2$spatial_score <- agg2$`Mean Dist` * dist.weight
  
  agg2$ndvi_norm <- normalize(agg2$ndvi_score) # normalize
  agg2$igno_norm <- normalize(agg2$igno_score)
  agg2$spatial_norm <- normalize(agg2$spatial_score)
  
  agg2$FINAL_SCORE <- agg2$ndvi_norm + agg2$igno_norm + agg2$spatial_norm
  
  # remove spatial configurations with intersection
  agg2 <- agg2[agg2$INTERSECTION == "FALSE",] 
  
  # --- Robustness Check ---
  # Check if any valid solutions remain
  if (nrow(agg2) == 0) {
    close(pb)
    warning("No valid (non-intersecting) solutions found. Try increasing 'perm' or decreasing 'areaplot'.")
    return(list("Full matrix" = new_mat, "Aggregated matrix" = agg2, "Best" = NA))
  }
  # ---
  
  ordered_solutions <- agg2[order(agg2[, 'FINAL_SCORE'], decreasing = TRUE),]
  Index <- as.numeric(ordered_solutions[1, 1])
  sol <- subset(new_mat[new_mat$try %in% Index,])
  sol2 <- subset(agg2[agg2$Try %in% Index,])
  
  ## Plot best solution
  
  # Create sf objects from the best solution data
  out1_points <- sf::st_as_sf(sol, coords = c("x", "y"), crs = sf::st_crs(boundary))
  out1_buffers <- sf::st_buffer(out1_points, dist = sqrt(areaplot / pi))
  
  # Reproject site to match raster CRS for plotting
  site2 <- sf::st_transform(site, newproj)
  
  # Plot using sf/terra
  plot(terra::rast(ndvi_map, 1), main = "Best Solution") # Plot first layer of ndvi as basemap
  plot(sf::st_geometry(site2), add = TRUE, border = "blue", lwd = 2)
  plot(sf::st_geometry(out1_buffers), add = TRUE, border = "red", col = NA)
  
  
  # --- Modern Plotting with ggplot2 + tidyterra ---
  
  # Plot 1: NDVI
  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = ndvi) +
    ggplot2::geom_sf(data = out1_points, color = "black", size = 2) +
    ggplot2::scale_fill_distiller(palette = "YlGn", name = "NDVI", na.value = "transparent") +
    ggplot2::ggtitle("Best Solution on NDVI") +
    ggplot2::theme_minimal()
  
  # Plot 2: Ignorance
  p1 <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = ignorance) +
    ggplot2::geom_sf(data = out1_points, color = "darkgray", size = 1.5) +
    ggplot2::scale_fill_distiller(palette = "RdBu", name = "Ignorance", na.value = "transparent") +
    ggplot2::ggtitle("Best Solution on Ignorance Map") +
    ggplot2::theme_minimal()
  
  
  if (nrow(new_mat) > 1000) { new_mat <- dplyr::sample_n(new_mat, 1000) }
  
  # Plot 3: Density
  p2 <-  ggplot2::ggplot(new_mat, ggplot2::aes(x = ndvi, group = try)) +
    ggplot2::geom_density(colour = "lightgrey") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_density(data = sol, ggplot2::aes(x = ndvi), colour = "red", size = 1) +
    ggplot2::ggtitle("Distribution of NDVI values")
  
  
  # Plot 4: 3D Scatter (unchanged, but index is now safer)
  index_graph <- match(max(agg2$FINAL_SCORE), agg2$FINAL_SCORE)
  
  plot3D::scatter3D(agg2$ndvi_norm, agg2$igno_norm, agg2$spatial_norm, 
                    bty = "b2", colvar = agg2$FINAL_SCORE, 
                    xlab = "NDVI", ylab = "IGNORANCE", zlab = "SPACE", 
                    clab = c("Multiobjective", "Sampling Optimisation"),
                    main = "Optimization Solution Space")
  
  plot3D::scatter3D(x = agg2$ndvi_norm[index_graph], y = agg2$igno_norm[index_graph], z = agg2$spatial_norm[index_graph], 
                    add = TRUE, colkey = FALSE, 
                    pch = 18, cex = 3, col = "black")
  
  end_time <- Sys.time()
  message()
  message(paste0("Total runtime: ", round(end_time - start_time, 2), " seconds"))
  
  
  return(list("Full matrix" = new_mat, "Aggregated matrix" = agg2, "Best" = sol, "Variance of sampling points" = sol2[, 'Variance'],
              "Mean Ignorance" = sol2[, 'igno_score'],
              "Spatial Median of Distance" = sol2[, 'Mean Dist'], "Final score" = sol2[, 'FINAL_SCORE'], 
              "Plot_NDVI" = p, "Plot_Ignorance" = p1, "Plot_Density" = p2))
  
  
}