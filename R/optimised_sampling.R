#' @title Boosted Sampling Optimization
#'
#' @description
#' Performs multi-objective sampling optimization by generating multiple 
#' permutations of random sampling points and evaluating each configuration
#' based on weighted scores for:
#' \itemize{
#'  \item NDVI variance (environmental heterogeneity, between and within plots)
#'  \item Mean ignorance value (prioritizing understudied areas)
#'  \item Spatial dispersion (using mean nearest neighbor distance)
#'  \item DEM variance (optional, elevation heterogeneity)
#' }
#' Returns the best non-overlapping configuration with field-ready outputs.
#'
#' @param ndvi A `terra SpatRaster` object representing NDVI or other environmental index.
#' @param ignorance A `terra SpatRaster` object representing the Map of Relative Floristic Ignorance (MRFI).
#' @param site An `sf` polygon object defining the study area boundary.
#' @param excl_areas Optional `sf` object (or 'SpatialPolygonsDataFrame') to 
#'   delimit unsuitable areas (e.g., water bodies, inaccessible areas) to be 
#'   excluded from sampling. If no CRS is set, assumes EPSG:4326.
#' @param CRS.new Numeric EPSG code for projected CRS used in calculations (must 
#'   be in meters). Default = 3035 (ETRS89-LAEA Europe).
#' @param samp_strategy Character. Sampling strategy. Options: "random", "regular", "hexagonal". Default = "random".
#' @param nplot Integer. Number of sampling plots per permutation.
#' @param plot_radius Numeric. Radius of circular plots in meters (e.g., 20 for 20m radius plots, giving ~1,257m² area).
#' @param perm Integer. Number of permutations (sampling configurations) to test.
#' @param ndvi.weight Numeric. Weight for NDVI variance score (default = 1).
#' @param igno.weight Numeric. Weight for mean ignorance score (default = 1).
#' @param dist.weight Numeric. Weight for spatial dispersion score (default = 1).
#' @param dem A `terra SpatRaster` DEM (optional). If provided, elevation variance included in optimization.
#' @param dem.weight Numeric. Weight for DEM variance (default = 1, ignored if dem = NULL).
#' @param slope_max Numeric. Maximum slope in degrees for accessibility (e.g., 30). Requires dem.
#' @param extract_method Character. "auto" (default), "point", or "area" for raster value extraction.
#' @param within_var_weight Numeric. Weight for within-plot variance (0-1, default = 0.5).
#' @param verbose Logical. Print progress messages (default = TRUE).
#' @param output_dir Character. Directory for output files (default = working directory).
#' @param output_prefix Character. Prefix for output filenames (default = "SampleBoost").
#'
#' @return A list with full_matrix, aggregated_scores, best_solution, best_solution_sf,
#' best_variance, best_ignorance, best_dispersion, best_final_score, plots, and statistics.
#' 
#' @details
#' Uses mean nearest neighbor distance for spatial dispersion (100x faster than L-function).
#' 
#' **Extraction Method**: Automatically determined based on plot size relative to raster resolution:
#' \itemize{
#'   \item{Area-based extraction}{ Used when plot area > 4× raster cell area. Extracts all pixels 
#'   whose centroids fall within the circular plot buffer (centroid method). This avoids including 
#'   pixels that barely touch the plot edge, providing accurate representation of plot contents. 
#'   Calculates mean and variance from these pixels for within-plot heterogeneity assessment.}
#'   \item{Point extraction}{ Used for plots smaller than 4× raster cell area. Extracts single value 
#'   at plot center, with within-plot variance set to 0.}
#' }
#' 
#' For 90m resolution NDVI (8,100m² cells), area extraction is used for plots >32,400m² (radius >102m).
#' For 10m resolution data (100m² cells), area extraction is used for plots >400m² (radius >11m).
#' 
#' **Resolution Handling**: If DEM and NDVI have different resolutions, DEM is automatically 
#' resampled to match NDVI resolution using bilinear interpolation. This ensures fair comparison 
#' of variances and consistent pixel alignment across all rasters.
#' 
#' All rasters are cropped and masked to the site boundary.
#' 
#' **Boundary Constraint**: Plot centers are constrained to be at least plot_radius distance from 
#' the study area boundary, ensuring all plot buffers remain entirely within the study area. This 
#' is achieved through negative buffering of the sampling area.
#' 
#' **Adaptive Map Layout**: The function automatically detects study area shape and optimizes layout:
#' \itemize{
#'   \item{E-W oriented areas (aspect ratio > 1.3)}{ Maps stacked vertically (one above the other) 
#'   to use full width of landscape page}
#'   \item{N-S or square areas (aspect ratio ≤ 1.3)}{ Maps arranged side-by-side for efficient space use}
#' }
#' 
#' PDF Output Structure:
#' 
#' **Without DEM (3 pages, A4 landscape):**
#' \itemize{
#'   \item{Page 1}{ Best Sampling Configuration - NDVI and Ignorance maps (adaptive layout)}
#'   \item{Page 2}{ Optimization Diagnostics - Objective contribution bar chart showing why the 
#'   best solution won, plus multi-objective optimization space (NDVI vs Ignorance)}
#'   \item{Page 3}{ Summary Statistics - Complete parameter table}
#' }
#' 
#' **With DEM (4 pages, A4 landscape):**
#' \itemize{
#'   \item{Page 1}{ Best Sampling Configuration - NDVI and Ignorance maps (adaptive layout)}
#'   \item{Page 2}{ Best Sampling Configuration - DEM and Slope maps (adaptive layout)}
#'   \item{Page 3}{ Optimization Diagnostics - Objective contribution bar chart showing why the 
#'   best solution won, plus two multi-objective spaces (Environmental: NDVI vs Ignorance; 
#'   Topographic: DEM vs Ignorance)}
#'   \item{Page 4}{ Summary Statistics - Complete parameter table}
#' }
#' 
#' @importFrom sf st_transform st_crs st_sample st_as_sf st_buffer st_intersects st_coordinates st_geometry st_difference st_union st_make_valid st_area
#' @importFrom terra extract crs project resample res terrain
#' @importFrom tidyterra geom_spatraster
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot geom_sf scale_fill_distiller ggtitle theme_minimal aes theme_classic labs coord_sf geom_point geom_col geom_text scale_y_continuous xlab ylab theme element_text element_blank
#' @importFrom scales percent
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom stats aggregate na.omit var dist
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.draw textGrob gpar
#' @importFrom gridExtra grid.arrange tableGrob
#' @export
#' 
#' @examples
#' \dontrun{
#' # Load example data
#' ndvi <- load_ndvi_example()
#' mrfi <- load_mrfi_example()
#' data(park)
#' data(unsuitablezone)
#' 
#' # Basic usage
#' result <- sampleboost(
#'   ndvi = ndvi,
#'   ignorance = mrfi,
#'   site = park,
#'   nplot = 50,
#'   plot_radius = 5.64,  # ~100m² plots
#'   perm = 100
#' )
#' 
#' # With exclusion areas (water bodies, inaccessible areas)
#' result <- sampleboost(
#'   ndvi = ndvi,
#'   ignorance = mrfi,
#'   site = park,
#'   excl_areas = unsuitablezone,
#'   nplot = 50,
#'   plot_radius = 5.64,
#'   perm = 100
#' )
#' 
#' # With DEM and slope constraints
#' dem <- load_dem_example()
#' result <- sampleboost(
#'   ndvi = ndvi,
#'   ignorance = mrfi,
#'   site = park,
#'   excl_areas = unsuitablezone,
#'   nplot = 50,
#'   plot_radius = 20,  # ~1,257m² plots
#'   perm = 100,
#'   dem = dem,
#'   slope_max = 30
#' )
#' 
#' # With custom CRS (e.g., UTM zone 32N for Italy)
#' result <- sampleboost(
#'   ndvi = ndvi,
#'   ignorance = mrfi,
#'   site = park,
#'   CRS.new = 32632,  # EPSG:32632 = WGS 84 / UTM zone 32N
#'   nplot = 50,
#'   plot_radius = 10,  # ~314m² plots
#'   perm = 100
#' )
#' 
#' # Access results
#' plot(result$plots$ndvi)
#' write.csv(result$best_solution, "best_sampling_points.csv")
#' }

sampleboost <- function(ndvi, ignorance, site, excl_areas = NULL,
                        CRS.new = 3035,
                        samp_strategy = "random", 
                        nplot, plot_radius, perm, 
                        ndvi.weight = 1, igno.weight = 1, dist.weight = 1,
                        dem = NULL, dem.weight = 1, slope_max = NULL,
                        extract_method = "auto", within_var_weight = 0.5,
                        verbose = TRUE, output_dir = getwd(), 
                        output_prefix = "SampleBoost") {
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  # Calculate plot area from radius
  areaplot <- pi * plot_radius^2
  
  normalize <- function(x) {
    x_clean <- x[!is.na(x)]
    if(length(x_clean) == 0) return(x)
    rng <- max(x_clean) - min(x_clean)
    if(rng == 0) return(rep(0.5, length(x)))
    (x - min(x_clean)) / rng
  }
  
  # ============================================================================
  # SECTION 1: INPUT VALIDATION
  # ============================================================================
  
  msg("Validating inputs...")
  
  if (missing(ndvi) || missing(ignorance) || missing(site) || 
      missing(nplot) || missing(plot_radius) || missing(perm)) {
    stop("Missing required arguments.")
  }
  
  if (nplot < 2) stop("'nplot' must be at least 2")
  if (plot_radius <= 0) stop("'plot_radius' must be positive")
  if (perm < 1) stop("'perm' must be at least 1")
  if (within_var_weight < 0 || within_var_weight > 1) stop("'within_var_weight' must be 0-1")
  if (!samp_strategy %in% c("random", "regular", "hexagonal")) stop("Invalid samp_strategy")
  if (!extract_method %in% c("auto", "point", "area")) stop("Invalid extract_method")
  if (!is.null(slope_max) && is.null(dem)) stop("'slope_max' requires 'dem'")
  
  msg("Inputs validated.")
  
  # ============================================================================
  # SECTION 2: CRS AND EXTENT ALIGNMENT
  # ============================================================================
  
  msg(paste("Reprojecting inputs to EPSG:", CRS.new, "..."))
  
  # Define target CRS
  crs_sf <- sf::st_crs(CRS.new)
  crs_terra <- paste0("EPSG:", CRS.new)
  
  # Reproject all rasters to target CRS
  if (!identical(terra::crs(ndvi), crs_terra)) {
    msg("  Reprojecting NDVI...")
    ndvi <- terra::project(ndvi, crs_terra)
  }
  
  if (!identical(terra::crs(ignorance), crs_terra)) {
    msg("  Reprojecting ignorance...")
    ignorance <- terra::project(ignorance, crs_terra)
  }
  
  # Resample ignorance to match NDVI resolution if needed
  if (!all(terra::res(ndvi) == terra::res(ignorance))) {
    msg("  Resampling ignorance to match NDVI resolution...")
    ignorance <- terra::resample(ignorance, ndvi, method = "bilinear")
  }
  
  use_dem <- !is.null(dem)
  if (use_dem) {
    msg("  Processing DEM...")
    if (!identical(terra::crs(dem), crs_terra)) {
      msg("    Reprojecting DEM...")
      dem <- terra::project(dem, crs_terra)
    }
    if (!all(terra::res(dem) == terra::res(ndvi))) {
      dem_res <- round(sqrt(prod(terra::res(dem))), 1)
      ndvi_res <- round(sqrt(prod(terra::res(ndvi))), 1)
      msg(paste0("    DEM resolution (", dem_res, "m) differs from NDVI (", ndvi_res, 
                 "m) - resampling DEM to match NDVI..."))
      dem <- terra::resample(dem, ndvi, method = "bilinear")
    }
    if (!is.null(slope_max)) slope <- terra::terrain(dem, "slope", unit = "degrees")
  }
  
  # Handle site (boundary) projection
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  if (is.na(sf::st_crs(site))) {
    msg("  Site has no CRS; assuming EPSG:4326")
    sf::st_crs(site) <- 4326
  }
  site_proj <- sf::st_transform(sf::st_make_valid(site), crs_sf)
  
  # ============================================================================
  # SECTION 2B: CREATE SAMPLING AREA WITH NEGATIVE BUFFER
  # ============================================================================
  
  msg(paste0("Creating sampling area with ", round(plot_radius, 1), "m buffer from boundary..."))
  
  # Create negative buffer to ensure plot buffers don't extend outside site
  site_buffered <- sf::st_buffer(site_proj, dist = -plot_radius)
  
  # Check if negative buffer resulted in valid geometry
  if (inherits(site_buffered, "sfc_GEOMETRY") && length(site_buffered) == 0) {
    stop(paste0("Negative buffer of ", round(plot_radius, 1), 
                "m resulted in empty geometry. Plot radius is too large for this study area. ",
                "Try reducing plot_radius or increasing study area size."))
  }
  
  # Check for area loss
  original_area <- as.numeric(sf::st_area(site_proj))
  buffered_area <- as.numeric(sf::st_area(site_buffered))
  area_loss_pct <- (1 - buffered_area / original_area) * 100
  
  if (area_loss_pct > 15) {
    warning(paste0("Negative buffer removes ", round(area_loss_pct, 1), 
                   "% of study area to prevent plots extending beyond boundary. ",
                   "Consider reducing plot_radius if this is too much area loss."))
  }
  
  msg(paste0("  Sampling area: ", round(buffered_area / 1e6, 2), " km² (", 
             round(area_loss_pct, 1), "% reduction from negative buffer)"))
  
  # Handle exclusion areas
  excl_proj <- NULL
  has_exclusions <- !is.null(excl_areas)
  
  if (has_exclusions) {
    msg("  Processing exclusion areas...")
    
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    
    if (is.na(sf::st_crs(excl_areas))) {
      msg("  Exclusion areas have no CRS; assuming EPSG:4326")
      sf::st_crs(excl_areas) <- 4326
    }
    
    excl_proj <- sf::st_union(sf::st_transform(sf::st_make_valid(excl_areas), crs_sf))
    
    # Create sampling region: buffered site minus exclusion areas
    msg("  Creating sampling region (buffered site minus exclusions)...")
    site_sampling <- sf::st_difference(site_buffered, excl_proj)
    site_sampling <- sf::st_make_valid(site_sampling)
  } else {
    site_sampling <- site_buffered
  }
  
  # Calculate plot area and check raster resolution
  raster_cell_area <- prod(terra::res(ndvi))
  plot_area <- areaplot  # Already calculated from plot_radius at function start
  
  if (extract_method == "auto") {
    # Use area extraction if plot covers more than 4 raster cells
    use_area_extraction <- plot_area > raster_cell_area * 4
  } else {
    use_area_extraction <- (extract_method == "area")
  }
  
  msg(paste0("  Plot radius: ", round(plot_radius, 2), "m (area: ", round(areaplot, 1), "m²)"))
  msg(paste0("  Raster resolution: ", round(sqrt(raster_cell_area), 1), "m (cell area: ", round(raster_cell_area, 1), "m²)"))
  if (use_area_extraction) {
    msg("  Extraction method: Area-based (centroid method - only pixels whose centers fall within plots)")
  } else {
    msg("  Extraction method: Point-based (single value at plot center)")
  }
  if (has_exclusions) {
    msg("  Exclusion areas will be avoided during sampling")
  }
  
  # ============================================================================
  # SECTION 2C: MASK RASTERS TO SITE BOUNDARY
  # ============================================================================
  
  msg("Masking rasters to site boundary...")
  
  # Create site vector for masking (use original site_proj, not buffered)
  site_vect <- terra::vect(site_proj)
  
  # Mask and crop NDVI - this removes data outside site AND crops extent
  ndvi <- terra::crop(ndvi, site_vect)
  ndvi <- terra::mask(ndvi, site_vect)
  
  # Mask and crop ignorance
  ignorance <- terra::crop(ignorance, site_vect)
  ignorance <- terra::mask(ignorance, site_vect)
  
  # Mask and crop DEM if provided
  if (use_dem) {
    dem <- terra::crop(dem, site_vect)
    dem <- terra::mask(dem, site_vect)
    
    # Calculate slope for visualization (always, regardless of slope_max)
    slope_vis <- terra::terrain(dem, "slope", unit = "degrees")
  }
  
  # ============================================================================
  # SECTION 3: GENERATE AND EVALUATE PERMUTATIONS
  # ============================================================================
  
  msg(paste0("Generating ", perm, " permutations..."))
  
  permutation_results <- vector("list", perm)
  mean_nn_dist <- numeric(perm)
  has_intersection <- logical(perm)
  points_in_exclusion <- logical(perm)
  
  pb <- utils::txtProgressBar(min = 0, max = perm, style = 3)
  
  for (i in 1:perm) {
    # Sample from site_sampling (which excludes excl_areas if present)
    sample_points <- sf::st_sample(site_sampling, size = nplot, type = samp_strategy, iter = 10)
    points_sf <- sf::st_as_sf(sample_points)
    
    # Additional check: filter out any points that somehow ended up in exclusion areas
    if (has_exclusions) {
      in_excl <- sf::st_intersects(points_sf, excl_proj, sparse = FALSE)[, 1]
      if (any(in_excl)) {
        points_in_exclusion[i] <- TRUE
        points_sf <- points_sf[!in_excl, ]
        if (nrow(points_sf) < nplot) {
          # Not enough valid points, mark as invalid
          has_intersection[i] <- TRUE
        }
      }
    }
    
    point_buffers <- sf::st_buffer(points_sf, dist = plot_radius)
    
    # Check for overlapping buffers
    intersection_matrix <- sf::st_intersects(point_buffers, sparse = FALSE)
    diag(intersection_matrix) <- FALSE
    has_intersection[i] <- has_intersection[i] || any(intersection_matrix)
    
    if (use_area_extraction) {
      # Extract using centroid method: only include pixels whose centroids fall within plot
      # This avoids including pixels that barely touch the circular plot edge
      ndvi_extract <- terra::extract(ndvi, point_buffers, fun = NULL, ID = TRUE, touches = FALSE)
      ndvi_stats <- stats::aggregate(ndvi_extract[[2]], by = list(ID = ndvi_extract$ID), 
                                     FUN = function(x) {
                                       x <- x[!is.na(x)]
                                       if(length(x) == 0) return(c(mean = NA, var = 0))
                                       c(mean = mean(x), var = stats::var(x))
                                     })
      ndvi_means <- ndvi_stats$x[,1]
      ndvi_within_vars <- ndvi_stats$x[,2]
      
      igno_extract <- terra::extract(ignorance, point_buffers, fun = NULL, ID = TRUE, touches = FALSE)
      igno_stats <- stats::aggregate(igno_extract[[2]], by = list(ID = igno_extract$ID), 
                                     FUN = function(x) {
                                       x <- x[!is.na(x)]
                                       if(length(x) == 0) return(c(mean = NA, var = 0))
                                       c(mean = mean(x), var = stats::var(x))
                                     })
      igno_means <- igno_stats$x[,1]
      
      if (use_dem) {
        dem_extract <- terra::extract(dem, point_buffers, fun = NULL, ID = TRUE, touches = FALSE)
        dem_stats <- stats::aggregate(dem_extract[[2]], by = list(ID = dem_extract$ID), 
                                      FUN = function(x) {
                                        x <- x[!is.na(x)]
                                        if(length(x) == 0) return(c(mean = NA, var = 0))
                                        c(mean = mean(x), var = stats::var(x))
                                      })
        dem_means <- dem_stats$x[,1]
        dem_within_vars <- dem_stats$x[,2]
        
        if (!is.null(slope_max)) {
          slope_values <- terra::extract(slope, point_buffers, fun = max, na.rm = TRUE, ID = FALSE, touches = FALSE)
          if (any(slope_values > slope_max, na.rm = TRUE)) has_intersection[i] <- TRUE
        }
      }
    } else {
      ndvi_means <- terra::extract(ndvi, points_sf, ID = FALSE)[,1]
      ndvi_within_vars <- rep(0, length(ndvi_means))
      igno_means <- terra::extract(ignorance, points_sf, ID = FALSE)[,1]
      
      if (use_dem) {
        dem_means <- terra::extract(dem, points_sf, ID = FALSE)[,1]
        dem_within_vars <- rep(0, length(dem_means))
        if (!is.null(slope_max)) {
          slope_values <- terra::extract(slope, points_sf, ID = FALSE)[,1]
          if (any(slope_values > slope_max, na.rm = TRUE)) has_intersection[i] <- TRUE
        }
      }
    }
    
    coords <- sf::st_coordinates(points_sf)
    dist_matrix <- as.matrix(stats::dist(coords))
    diag(dist_matrix) <- Inf
    mean_nn_dist[i] <- mean(apply(dist_matrix, 1, min))
    
    permutation_results[[i]] <- data.frame(
      x = coords[,1], y = coords[,2],
      ndvi_mean = ndvi_means, ndvi_within_var = ndvi_within_vars,
      ignorance_mean = igno_means,
      elevation = if(use_dem) dem_means else rep(NA, nrow(points_sf)),
      elevation_within_var = if(use_dem) dem_within_vars else rep(0, nrow(points_sf)),
      has_intersection = has_intersection[i], permutation_id = i,
      stringsAsFactors = FALSE
    )
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  if (has_exclusions && any(points_in_exclusion)) {
    msg(paste0("  Note: ", sum(points_in_exclusion), " permutations had points in exclusion areas (filtered)"))
  }
  
  # ============================================================================
  # SECTION 4: AGGREGATE AND SCORE
  # ============================================================================
  
  msg("Scoring permutations...")
  
  full_matrix <- dplyr::bind_rows(permutation_results)
  full_matrix$permutation_id <- as.factor(full_matrix$permutation_id)
  
  ndvi_between <- stats::aggregate(full_matrix$ndvi_mean, by = list(permutation_id = full_matrix$permutation_id), FUN = stats::var)
  ndvi_within <- stats::aggregate(full_matrix$ndvi_within_var, by = list(permutation_id = full_matrix$permutation_id), FUN = mean)
  mean_ignorance <- stats::aggregate(full_matrix$ignorance_mean, by = list(permutation_id = full_matrix$permutation_id), FUN = mean)
  
  if (use_dem) {
    dem_between <- stats::aggregate(full_matrix$elevation, by = list(permutation_id = full_matrix$permutation_id), FUN = stats::var, na.rm = TRUE)
    dem_within <- stats::aggregate(full_matrix$elevation_within_var, by = list(permutation_id = full_matrix$permutation_id), FUN = mean)
  } else {
    dem_between <- data.frame(permutation_id = 1:perm, x = 0)
    dem_within <- data.frame(permutation_id = 1:perm, x = 0)
  }
  
  aggregated_scores <- data.frame(
    permutation_id = 1:perm, ndvi_between_var = ndvi_between$x, ndvi_within_var = ndvi_within$x,
    dem_between_var = dem_between$x, dem_within_var = dem_within$x,
    mean_nn_dist = mean_nn_dist, mean_ignorance = mean_ignorance$x,
    has_intersection = has_intersection, stringsAsFactors = FALSE
  )
  
  aggregated_scores <- stats::na.omit(aggregated_scores)
  
  aggregated_scores$ndvi_combined <- aggregated_scores$ndvi_between_var + within_var_weight * aggregated_scores$ndvi_within_var
  aggregated_scores$ndvi_weighted <- aggregated_scores$ndvi_combined * ndvi.weight
  aggregated_scores$igno_weighted <- aggregated_scores$mean_ignorance * igno.weight
  aggregated_scores$spatial_weighted <- aggregated_scores$mean_nn_dist * dist.weight
  
  if (use_dem) {
    aggregated_scores$dem_combined <- aggregated_scores$dem_between_var + within_var_weight * aggregated_scores$dem_within_var
    aggregated_scores$dem_weighted <- aggregated_scores$dem_combined * dem.weight
  } else {
    aggregated_scores$dem_weighted <- 0
  }
  
  aggregated_scores$ndvi_norm <- normalize(aggregated_scores$ndvi_weighted)
  aggregated_scores$igno_norm <- normalize(aggregated_scores$igno_weighted)
  aggregated_scores$spatial_norm <- normalize(aggregated_scores$spatial_weighted)
  aggregated_scores$dem_norm <- if(use_dem) normalize(aggregated_scores$dem_weighted) else 0
  
  aggregated_scores$final_score <- aggregated_scores$ndvi_norm + aggregated_scores$igno_norm + 
    aggregated_scores$spatial_norm + aggregated_scores$dem_norm
  
  # ============================================================================
  # SECTION 5: SELECT BEST SOLUTION
  # ============================================================================
  
  msg("Selecting best configuration...")
  
  valid_solutions <- aggregated_scores[!aggregated_scores$has_intersection, ]
  if (nrow(valid_solutions) == 0) {
    warning("No valid solutions. Increase perm or decrease plot_radius/nplot.")
    return(list(full_matrix = full_matrix, aggregated_scores = aggregated_scores, 
                best_solution = NULL, best_solution_sf = NULL))
  }
  
  best <- valid_solutions[order(valid_solutions$final_score, decreasing = TRUE),][1,]
  best_solution <- full_matrix[full_matrix$permutation_id == best$permutation_id, ]
  best_solution_sf <- sf::st_as_sf(best_solution, coords = c("x","y"), crs = sf::st_crs(site_proj))
  best_buffers_sf <- sf::st_buffer(best_solution_sf, dist = plot_radius)
  
  # ============================================================================
  # SECTION 6: STATISTICS
  # ============================================================================
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  statistics <- data.frame(
    Statistic = c("Started", "Finished", "Elapsed time (secs)", "CRS (EPSG)", "Sampling strategy",
                  "Exclusion areas applied",
                  "Number of plots", "Plot area (m²)", "Plot radius (m)", "Permutations tested", "Valid solutions",
                  "Extraction method", "Raster resolution (m)", "Within-plot var weight",
                  "DEM used", "Slope threshold (°)", "NDVI weight", "Ignorance weight", "Distance weight", "DEM weight",
                  "Best permutation ID", "Best final score", "Best NDVI variance", "Best mean ignorance",
                  "Best mean NN distance (m)", "Best elevation variance"),
    Value = c(as.character(start_time), as.character(end_time), round(elapsed, 2),
              as.character(sf::st_crs(site_proj)$epsg), samp_strategy,
              if (has_exclusions) "Yes" else "No",
              nplot, round(areaplot, 2), round(plot_radius, 2),
              perm, nrow(valid_solutions), ifelse(use_area_extraction, "Area", "Point"),
              round(sqrt(raster_cell_area), 1), within_var_weight, ifelse(use_dem, "Yes", "No"),
              ifelse(is.null(slope_max), "None", slope_max), ndvi.weight, igno.weight, dist.weight,
              ifelse(use_dem, dem.weight, "N/A"), best$permutation_id, round(best$final_score, 4),
              round(best$ndvi_combined, 4), round(best$mean_ignorance, 2), round(best$mean_nn_dist, 2),
              ifelse(use_dem, round(best$dem_combined, 2), "N/A")),
    stringsAsFactors = FALSE
  )
  
  # ============================================================================
  # SECTION 7: PLOTS AND OUTPUT
  # ============================================================================
  
  msg("Creating outputs...")
  
  # Helper function to add boundaries and exclusion areas to plots
  add_boundaries <- function(p) {
    p <- p + ggplot2::geom_sf(data = site_proj, fill = NA, color = "black", linewidth = 1, inherit.aes = FALSE)
    if (has_exclusions) {
      p <- p + ggplot2::geom_sf(data = excl_proj, fill = "lightgray", color = "gray50", 
                                linewidth = 0.5, linetype = "dotted", inherit.aes = FALSE, alpha = 0.3)
    }
    return(p)
  }
  
  plot_ndvi <- ggplot2::ggplot() + tidyterra::geom_spatraster(data = ndvi) +
    ggplot2::geom_sf(data = best_buffers_sf, fill = NA, color = "white", linewidth = 1.5) +
    ggplot2::geom_sf(data = best_solution_sf, color = "white", size = 2.5, shape = 21, fill = "black", stroke = 1.2) +
    ggplot2::scale_fill_distiller(palette = "YlGn", name = "NDVI", na.value = "transparent", direction = 1) +
    ggplot2::ggtitle("Best Sampling Configuration on NDVI") + ggplot2::coord_sf() + ggplot2::theme_minimal()
  plot_ndvi <- add_boundaries(plot_ndvi)
  
  plot_ignorance <- ggplot2::ggplot() + tidyterra::geom_spatraster(data = ignorance) +
    ggplot2::geom_sf(data = best_buffers_sf, fill = NA, color = "white", linewidth = 1.5) +
    ggplot2::geom_sf(data = best_solution_sf, color = "white", size = 2.5, shape = 21, fill = "black", stroke = 1.2) +
    ggplot2::scale_fill_distiller(palette = "Spectral", name = "Ignorance", na.value = "transparent", direction = -1) +
    ggplot2::ggtitle("Best Sampling Configuration on Ignorance Map") + ggplot2::coord_sf() + ggplot2::theme_minimal()
  plot_ignorance <- add_boundaries(plot_ignorance)
  
  # Create DEM and Slope plots if DEM is used
  if (use_dem) {
    plot_dem <- ggplot2::ggplot() + tidyterra::geom_spatraster(data = dem) +
      ggplot2::geom_sf(data = best_buffers_sf, fill = NA, color = "white", linewidth = 1.5) +
      ggplot2::geom_sf(data = best_solution_sf, color = "white", size = 2.5, shape = 21, fill = "black", stroke = 1.2) +
      ggplot2::scale_fill_distiller(palette = "YlOrBr", name = "Elevation (m)", na.value = "transparent", direction = 1) +
      ggplot2::ggtitle("Best Sampling Configuration on DEM") + ggplot2::coord_sf() + ggplot2::theme_minimal()
    plot_dem <- add_boundaries(plot_dem)
    
    plot_slope <- ggplot2::ggplot() + tidyterra::geom_spatraster(data = slope_vis) +
      ggplot2::geom_sf(data = best_buffers_sf, fill = NA, color = "white", linewidth = 1.5) +
      ggplot2::geom_sf(data = best_solution_sf, color = "white", size = 2.5, shape = 21, fill = "black", stroke = 1.2) +
      ggplot2::scale_fill_distiller(palette = "YlOrBr", name = "Slope (°)", na.value = "transparent", direction = 1) +
      ggplot2::ggtitle("Best Sampling Configuration on Slope") + ggplot2::coord_sf() + ggplot2::theme_minimal()
    plot_slope <- add_boundaries(plot_slope)
  }
  
  plot_scores <- ggplot2::ggplot(valid_solutions, ggplot2::aes(x=ndvi_norm, y=igno_norm, size=spatial_norm, color=final_score)) +
    ggplot2::geom_point(alpha=0.7) + ggplot2::geom_point(data = best, color = "black", size = 6, shape = 18) +
    ggplot2::scale_color_distiller(palette = "Spectral", name = "Final Score") + ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold")) +
    ggplot2::ggtitle("Multi-Objective Optimization Space") + ggplot2::xlab("NDVI Variance (Normalized)") +
    ggplot2::ylab("Ignorance (Normalized)") + ggplot2::labs(size="Dispersion (Normalized)")
  
  # Add DEM-focused optimization plot if DEM is used
  if (use_dem) {
    plot_scores_dem <- ggplot2::ggplot(valid_solutions, ggplot2::aes(x=dem_norm, y=igno_norm, size=ndvi_norm, color=final_score)) +
      ggplot2::geom_point(alpha=0.7) + ggplot2::geom_point(data = best, color = "black", size = 6, shape = 18) +
      ggplot2::scale_color_distiller(palette = "Spectral", name = "Final Score") + ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold")) +
      ggplot2::ggtitle("Topographic Optimization Space") + ggplot2::xlab("DEM Variance (Normalized)") +
      ggplot2::ylab("Ignorance (Normalized)") + ggplot2::labs(size="NDVI Variance (Normalized)")
  }
  
  # Create objective contribution bar chart for best solution
  contributions <- data.frame(
    Objective = c("NDVI Variance", "Ignorance", "Spatial Dispersion", if(use_dem) "DEM Variance" else NULL),
    Value = c(best$ndvi_norm, best$igno_norm, best$spatial_norm, if(use_dem) best$dem_norm else NULL),
    stringsAsFactors = FALSE
  )
  contributions$Objective <- factor(contributions$Objective, levels = rev(contributions$Objective))
  
  plot_contributions <- ggplot2::ggplot(contributions, ggplot2::aes(x = Objective, y = Value)) +
    ggplot2::geom_col(fill = "#4CAF50", alpha = 0.8, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f (%.0f%%)", Value, Value*100)), 
                       hjust = -0.1, size = 4, fontface = "bold") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25), labels = scales::percent) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 11, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = sprintf("Best Solution #%d - Objective Contributions (Final Score: %.2f)", 
                      best$permutation_id, best$final_score),
      x = NULL,
      y = "Normalized Value (% of Maximum)"
    )
  
  # Note: Spatial dispersion is shown via point size in scatter plots
  # and included in the objective contribution bars
  
  csv_best <- file.path(output_dir, paste0(output_prefix, "_best-solution.csv"))
  csv_all <- file.path(output_dir, paste0(output_prefix, "_all-scores.csv"))
  csv_field <- file.path(output_dir, paste0(output_prefix, "_field-sheet.csv"))
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_plots.pdf"))
  
  utils::write.csv(best_solution, csv_best, row.names = FALSE)
  utils::write.csv(aggregated_scores, csv_all, row.names = FALSE)
  
  field_sheet <- data.frame(
    Plot_ID = sprintf("PLOT_%03d", 1:nrow(best_solution)),
    Longitude_WGS84 = round(sf::st_coordinates(sf::st_transform(best_solution_sf, 4326))[,1], 6),
    Latitude_WGS84 = round(sf::st_coordinates(sf::st_transform(best_solution_sf, 4326))[,2], 6),
    NDVI = round(best_solution$ndvi_mean, 3), Ignorance = round(best_solution$ignorance_mean, 2),
    stringsAsFactors = FALSE
  )
  if (use_dem) field_sheet$Elevation_m <- round(best_solution$elevation, 0)
  field_sheet$Notes <- ""
  utils::write.csv(field_sheet, csv_field, row.names = FALSE)
  
  # Detect study area aspect ratio for optimal layout
  site_bbox <- sf::st_bbox(site_proj)
  bbox_width <- site_bbox$xmax - site_bbox$xmin
  bbox_height <- site_bbox$ymax - site_bbox$ymin
  aspect_ratio <- bbox_width / bbox_height
  
  # Determine layout: E-W areas (wide) use vertical stacking, N-S/square use horizontal
  if (aspect_ratio > 1.3) {
    map_layout <- 1  # ncol=1, vertical stacking for wide areas
    layout_type <- "vertical (one above the other)"
    msg(paste0("  Study area is E-W oriented (aspect ratio: ", round(aspect_ratio, 2), 
               ") - stacking maps vertically for better visibility"))
  } else {
    map_layout <- 2  # ncol=2, side-by-side for tall/square areas
    layout_type <- "horizontal (side-by-side)"
    msg(paste0("  Study area is N-S/square (aspect ratio: ", round(aspect_ratio, 2), 
               ") - arranging maps side-by-side"))
  }
  
  # Generate PDF with adaptive layout
  grDevices::pdf(pdf_path, width = 11.69, height = 8.27, onefile = TRUE)
  
  # Page 1: NDVI + Ignorance (adaptive layout)
  gridExtra::grid.arrange(plot_ndvi, plot_ignorance, ncol = map_layout,
                          top = grid::textGrob("Page 1: Best Sampling Configuration (Environmental Layers)",
                                               gp = grid::gpar(fontsize = 14, fontface = "bold")))
  
  if (use_dem) {
    # Page 2: DEM + Slope (adaptive layout)
    gridExtra::grid.arrange(plot_dem, plot_slope, ncol = map_layout,
                            top = grid::textGrob("Page 2: Best Sampling Configuration (Topographic Layers)",
                                                 gp = grid::gpar(fontsize = 14, fontface = "bold")))
    
    # Page 3: Optimization Diagnostics (with DEM: contribution bars + 2 optimization plots)
    gridExtra::grid.arrange(
      plot_contributions,
      gridExtra::arrangeGrob(plot_scores, plot_scores_dem, ncol = 2),
      ncol = 1,
      heights = c(0.5, 1),
      top = grid::textGrob("Page 3: Optimization Diagnostics",
                           gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
    
    # Page 4: Summary Statistics
    grid::grid.draw(gridExtra::grid.arrange(top = grid::textGrob("Page 4: Summary Statistics",
                                                                 gp = grid::gpar(fontsize = 14, fontface = "bold")),
                                            gridExtra::tableGrob(statistics)))
  } else {
    # Page 2: Optimization Diagnostics (no DEM)
    gridExtra::grid.arrange(
      plot_contributions,
      plot_scores,
      ncol = 1,
      heights = c(0.5, 1),
      top = grid::textGrob("Page 2: Optimization Diagnostics",
                           gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
    
    # Page 3: Summary Statistics (no DEM)
    grid::grid.draw(gridExtra::grid.arrange(top = grid::textGrob("Page 3: Summary Statistics",
                                                                 gp = grid::gpar(fontsize = 14, fontface = "bold")),
                                            gridExtra::tableGrob(statistics)))
  }
  
  grDevices::dev.off()
  
  msg(paste0("Done! Best permutation: #", best$permutation_id, " (score: ", round(best$final_score, 2), ")"))
  msg(paste0("Files saved to: ", output_dir))
  if (use_dem) {
    msg(paste0("PDF structure: 4 pages (", layout_type, " layout for maps)"))
  } else {
    msg(paste0("PDF structure: 3 pages (", layout_type, " layout for maps)"))
  }
  
  # ============================================================================
  # SECTION 8: RETURN
  # ============================================================================
  
  return(list(
    full_matrix = full_matrix, aggregated_scores = aggregated_scores,
    best_solution = best_solution, best_solution_sf = best_solution_sf,
    best_variance = best$ndvi_combined, best_ignorance = best$mean_ignorance,
    best_dispersion = best$mean_nn_dist, best_final_score = best$final_score,
    plots = if(use_dem) {
      list(ndvi = plot_ndvi, ignorance = plot_ignorance, dem = plot_dem, slope = plot_slope, 
           scores = plot_scores, scores_dem = plot_scores_dem, contributions = plot_contributions)
    } else {
      list(ndvi = plot_ndvi, ignorance = plot_ignorance, scores = plot_scores, contributions = plot_contributions)
    },
    statistics = statistics
  ))
}