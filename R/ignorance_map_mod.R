#' @title Map of Relative Floristic Ignorance (Modernized)
#'
#' @description Computes a Map of Relative Floristic Ignorance (MRFI) using
#' modern `sf` and `terra` packages for high performance. Aligned with MRFI 
#' philosophy: shows all areas where data exists, letting IRFI values indicate 
#' sampling quality rather than excluding poorly-sampled areas.
#'
#' @param data_flor A data.frame with 5 columns: 'Taxon' (species identity), 
#' 'Long' (longitude), 'Lat' (latitude), 'uncertainty' (spatial uncertainty 
#' radius in meters), and 'year' (year of occurrence record).
#' @param site An `sf` object (or 'SpatialPolygonsDataFrame') representing the 
#' study area. If no CRS is set, assumes EPSG:4326 (WGS84).
#' @param year_study Numeric year of analysis (e.g., 2025). Defaults to current 
#' system year if not specified.
#' @param excl_areas Optional `sf` object (or 'SpatialPolygonsDataFrame') to 
#' delimit unsuitable areas (e.g., marine surfaces for terrestrial flora) to be 
#' excluded from calculations. If no CRS is set, assumes EPSG:4326.
#' @param CRS.new Numeric EPSG code for projected CRS used in calculations (must 
#' be in meters). Default = 3035 (ETRS89-LAEA Europe).
#' @param tau Numeric. Percentual value of taxa loss over 100 years (0 <= tau < 100). 
#' Represents the temporal decay of floristic records.
#' @param cellsize Numeric. Resolution of output raster in meters.
#' @param verbose Logical. If TRUE (default), prints progress messages during execution.
#' @param check_overlap Logical. If TRUE (default), checks and plots spatial overlap 
#' between occurrence points and study area.
#' @param output_dir Character. Directory path for output files. Defaults to working directory.
#' @param output_prefix Character. Prefix for output filenames. Default = "MRFI".
#' @param site_buffer Logical. If TRUE, expands the study area for both processing 
#' and final masking. If FALSE (default), uses original boundary. Buffer width 
#' controlled by buffer_width parameter.
#' @param buffer_width Numeric. Width of buffer in meters. If NULL (default), uses 
#' cellsize when site_buffer=TRUE. If specified, overrides default. Must be > 0. 
#' Buffer ensures complete raster coverage and avoids edge artifacts. Buffer will 
#' NOT extend into exclusion areas.
#' @param mask_method Character. Method for final raster masking. One of:
#'   \itemize{
#'     \item{"touches"}{ Include any cell intersecting site boundary (default, MRFI-aligned)}
#'     \item{"none"}{ Include entire raster extent, no masking (shows all data)}
#'   }
#' @param use_coverage_weighting Logical. If TRUE (default), weights ignorance 
#' contribution by fraction of cell covered by each buffer (matches original algorithm, 
#' slower but more accurate). If FALSE, uses faster binary approach where any buffer 
#' touch contributes full ignorance value (faster but may underestimate ignorance).
#' @param color_scale Character. Color scale method for MRFI plot. One of:
#'   \itemize{
#'     \item{"continuous"}{ Continuous color gradient (default)}
#'     \item{"quantile"}{ Colors based on quantile breaks (highlights sampling differences)}
#'   }
#' @param n_quantiles Numeric. Number of quantile breaks when color_scale="quantile". 
#' Default = 8 (octiles). Ignored if color_scale="continuous".
#'
#' @return A list with 4 objects:
#' \itemize{
#'   \item{MRFI}{ A `terra SpatRaster` of the Map of Relative Floristic Ignorance}
#'   \item{RICH}{ A `terra SpatRaster` of species richness computed without uncertainties}
#'   \item{Uncertainties}{ A data.frame of uncertainty and year values for all records 
#'   used in computation}
#'   \item{Statistics}{ A data.frame summarizing settings and results}
#' }
#'
#' @importFrom sf st_as_sf st_make_valid st_crs st_transform st_buffer st_intersects st_intersection st_union st_bbox st_coordinates st_geometry st_difference st_sf
#' @importFrom terra rast values ncell cellFromXY mask crop crs vect rasterize extract ext global writeRaster app
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom ggplot2 ggplot aes geom_tile geom_sf coord_equal theme_classic labs theme xlab ylab scale_fill_distiller scale_fill_gradientn guide_legend guide_colorbar ggtitle guides geom_histogram coord_cartesian scale_y_continuous after_stat element_text
#' @importFrom grid unit grid.draw
#' @importFrom gridExtra grid.arrange tableGrob
#' @importFrom grDevices pdf dev.off
#' @importFrom stats quantile
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales rescale
#' 
#' @export
#' 
#' @examples \dontrun{
#' # Load example data (requires 'ignobioR' package data)
#' data(floratus)
#' data(park)
#' data(unsuitablezone)
#'
#' # Example 1: Default settings (no buffer, continuous scale, accurate)
#' mrfi <- ignorance_map_mod(
#'   data_flor = floratus,
#'   site = park,
#'   excl_areas = unsuitablezone,
#'   tau = 20,
#'   cellsize = 2000
#' )
#'
#' # Example 2: With buffer (defaults to cellsize width)
#' mrfi_buffered <- ignorance_map_mod(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 20,
#'   cellsize = 2000,
#'   site_buffer = TRUE  # Uses cellsize (2000m) as buffer
#' )
#'
#' # Example 3: Custom buffer distance (5km)
#' mrfi_5km <- ignorance_map_mod(
#'   data_flor = floratus,
#'   site = park,
#'   excl_areas = unsuitablezone,
#'   tau = 20,
#'   cellsize = 2000,
#'   site_buffer = TRUE,
#'   buffer_width = 5000
#' )
#'
#' # Example 4: Quantile color scale to highlight sampling differences
#' mrfi_quantile <- ignorance_map_mod(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 20,
#'   cellsize = 2000,
#'   color_scale = "quantile"  # Uses 8 quantiles by default (octiles)
#' )
#'
#' # Example 5: Fast mode without coverage weighting
#' mrfi_fast <- ignorance_map_mod(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 20,
#'   cellsize = 2000,
#'   use_coverage_weighting = FALSE
#' )
#'
#' # Example 6: No buffer, no masking (maximum extent)
#' mrfi_full <- ignorance_map_mod(
#'   data_flor = floratus,
#'   site = park,
#'   tau = 20,
#'   cellsize = 2000,
#'   site_buffer = FALSE,
#'   mask_method = "none"
#' )
#'
#' # View results
#' terra::plot(mrfi$MRFI)
#' terra::plot(mrfi$RICH)
#' print(mrfi$Statistics)
#' }

ignorance_map_mod <- function(data_flor, site, year_study = NULL, excl_areas = NULL,
                              CRS.new = 3035, tau, cellsize, verbose = TRUE,
                              check_overlap = TRUE, output_dir = getwd(),
                              output_prefix = "MRFI", site_buffer = FALSE,
                              buffer_width = NULL,
                              mask_method = "touches",
                              use_coverage_weighting = TRUE,
                              color_scale = "continuous",
                              n_quantiles = 8) {
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  # Conditional message printing based on verbose setting
  msg <- function(...) if (verbose) message(...)
  
  # Record start time for performance tracking
  start_time <- Sys.time()
  
  # ============================================================================
  # SECTION 1: SETTINGS AND INPUT VALIDATION
  # ============================================================================
  
  msg("Checking settings and inputs...")
  
  # Set year_study to current year if not specified
  if (is.null(year_study)) year_study <- as.numeric(format(Sys.Date(), "%Y"))
  
  # Validate year_study range
  if (year_study < 1800 || year_study > 2100) {
    warning("year_study seems unusual. Double-check the value.")
  }
  
  # Check for future dates in occurrence data
  if (max(data_flor$year) > year_study) {
    warning("Some occurrence dates are more recent than year_study.")
  }
  
  # Validate tau parameter (percentage loss rate)
  if (tau < 0 || tau >= 100) stop("0 <= tau < 100 is required.")
  
  # Validate cellsize parameter
  if (cellsize <= 0) stop("cellsize must be positive.")
  
  # Validate CRS parameter
  if (CRS.new <= 0 || !is.numeric(CRS.new)) {
    stop("CRS.new must be a valid positive EPSG code.")
  }
  
  # Validate and process site_buffer and buffer_width parameters
  if (!is.logical(site_buffer)) {
    stop("site_buffer must be TRUE or FALSE.")
  }
  
  if (site_buffer) {
    # Buffer enabled
    if (is.null(buffer_width)) {
      # Default: use cellsize as buffer width
      buffer_distance <- cellsize
      msg(paste0("Buffer enabled: using default width = cellsize (", cellsize, "m)..."))
    } else {
      # Custom buffer width specified
      if (!is.numeric(buffer_width) || buffer_width <= 0) {
        stop("buffer_width must be a positive numeric value (in meters).")
      }
      buffer_distance <- buffer_width
      msg(paste0("Buffer enabled: using custom width = ", buffer_distance, "m..."))
    }
    use_buffer <- TRUE
  } else {
    # No buffer requested
    if (!is.null(buffer_width)) {
      warning("buffer_width specified but site_buffer=FALSE. Buffer width will be ignored.")
    }
    buffer_distance <- 0
    use_buffer <- FALSE
  }
  
  # Validate mask_method parameter (only touches and none)
  valid_methods <- c("touches", "none")
  if (!mask_method %in% valid_methods) {
    stop("mask_method must be one of: ", paste(valid_methods, collapse = ", "), ".")
  }
  
  # Validate use_coverage_weighting parameter
  if (!is.logical(use_coverage_weighting)) {
    stop("use_coverage_weighting must be TRUE or FALSE.")
  }
  
  if (use_coverage_weighting) {
    msg("Coverage weighting: ENABLED (accurate, matches original algorithm, slower)...")
  } else {
    msg("Coverage weighting: DISABLED (faster, binary touch method)...")
  }
  
  # Validate color_scale parameter
  valid_color_scales <- c("continuous", "quantile")
  if (!color_scale %in% valid_color_scales) {
    stop("color_scale must be one of: ", paste(valid_color_scales, collapse = ", "), ".")
  }
  
  # Validate n_quantiles parameter
  if (color_scale == "quantile") {
    if (!is.numeric(n_quantiles) || n_quantiles < 2 || n_quantiles > 20) {
      stop("n_quantiles must be a numeric value between 2 and 20.")
    }
  }
  
  # Check required column names in data_flor
  req_cols <- c("Taxon", "Long", "Lat", "uncertainty", "year")
  if (!all(req_cols %in% names(data_flor))) {
    stop("data_flor must contain columns: ", paste(req_cols, collapse = ", "), ".")
  }
  
  # Store initial record count for statistics
  total_initial_records <- nrow(data_flor)
  
  # Warn if uncertainties are too small relative to cell size
  if (any(2 * data_flor$uncertainty < (cellsize / 20))) {
    warning("Some records have uncertainty very small relative to cellsize. ",
            "They may not be well-captured by the raster grid.")
  }
  
  msg("Inputs validated.")
  
  # ============================================================================
  # SECTION 2: COORDINATE SYSTEM SETUP AND REPROJECTION
  # ============================================================================
  
  msg(paste("Reprojecting inputs to EPSG:", CRS.new, "..."))
  
  # Define CRS in both sf and terra formats
  crs_sf <- sf::st_crs(CRS.new)
  crs_terra <- paste0("EPSG:", CRS.new)
  
  # Handle site polygon: convert from sp if needed
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  
  # Assume WGS84 if no CRS is set
  if (is.na(sf::st_crs(site))) {
    msg("Input 'site' has no CRS. Assuming EPSG:4326.")
    sf::st_crs(site) <- 4326
  }
  
  # Reproject site to target CRS and fix any geometry issues
  site_proj_original <- sf::st_transform(sf::st_make_valid(site), crs_sf)
  
  # Handle exclusion areas FIRST (before buffering site)
  excl_proj <- NULL
  has_exclusions <- !is.null(excl_areas)
  
  if (has_exclusions) {
    msg("Processing exclusion areas...")
    
    # Convert from sp if needed
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    
    # Assume WGS84 if no CRS is set
    if (is.na(sf::st_crs(excl_areas))) {
      msg("Input 'excl_areas' has no CRS. Assuming EPSG:4326.")
      sf::st_crs(excl_areas) <- 4326
    }
    
    # Reproject and union into single geometry
    excl_proj <- sf::st_union(sf::st_transform(sf::st_make_valid(excl_areas), crs_sf))
  }
  
  # ============================================================================
  # SECTION 3: SITE BUFFER LOGIC
  # ============================================================================
  
  # Create expanded boundary, but NOT into exclusion zones
  if (use_buffer) {
    msg(paste0("Creating buffered site boundary (buffer = ", buffer_distance, "m)..."))
    
    # Buffer the site
    site_proj_buffered <- sf::st_buffer(site_proj_original, dist = buffer_distance)
    
    # Remove exclusion areas from buffered site if they exist
    if (has_exclusions) {
      msg("  Removing exclusion areas from buffered site...")
      site_proj_buffered <- sf::st_difference(site_proj_buffered, excl_proj)
    }
    
    site_proj_buffered <- sf::st_make_valid(site_proj_buffered)
    
    # For processing: use buffered site
    # For masking: also use buffered site (consistency!)
    site_proj_processing <- site_proj_buffered
    site_proj_mask <- site_proj_buffered
    
  } else {
    msg("Using original site boundary (no buffer)...")
    site_proj_processing <- site_proj_original
    site_proj_mask <- site_proj_original
    site_proj_buffered <- NULL  # Set to NULL when no buffer is used
  }
  
  # Convert to terra vectors for rasterization
  site_vect_processing <- terra::vect(site_proj_processing)
  site_vect_mask <- terra::vect(site_proj_mask)
  
  # ============================================================================
  # SECTION 4: FLORISTIC DATA PREPARATION
  # ============================================================================
  
  msg("Standardizing input: data_flor...")
  
  # Handle different input types for floristic data
  if (inherits(data_flor, "sf")) {
    msg("data_flor is already an sf object.")
    pts_sf <- data_flor
  } else if (inherits(data_flor, "Spatial")) {
    msg("Converting sp object to sf...")
    pts_sf <- sf::st_as_sf(data_flor)
  } else if (is.data.frame(data_flor)) {
    msg("Converting data.frame to sf (using Long/Lat)...")
    required_cols <- c("Long", "Lat")
    if (!all(required_cols %in% names(data_flor))) {
      stop("If 'data_flor' is a data.frame, it must include 'Long' and 'Lat' columns.")
    }
    pts_sf <- sf::st_as_sf(data_flor, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)
  } else {
    stop("Unsupported data_flor type: must be data.frame, sf, or Spatial object.")
  }
  
  # Reproject points to target CRS
  pts_proj <- sf::st_transform(pts_sf, crs_sf)
  
  # Calculate points within ORIGINAL site for statistics
  pts_in_site <- sf::st_intersection(pts_proj, site_proj_original)
  
  # ============================================================================
  # SECTION 5: OPTIONAL OVERLAP CHECK AND VISUALIZATION
  # ============================================================================
  
  if (check_overlap) {
    # Test spatial intersection between points and site
    overlap_check <- sf::st_intersects(pts_proj, site_proj_processing, sparse = FALSE)
    n_overlap <- sum(overlap_check)
    msg(paste("Number of occurrence points overlapping the processing area:", n_overlap))
    
    if (n_overlap == 0) {
      warning("No points overlap the study area. Check projections or coordinates.")
    } else {
      # Create diagnostic plot showing boundaries
      plot(sf::st_geometry(site_proj_original), col = NA, border = "black", 
           main = "Points vs Site Boundaries", lwd = 2)
      
      if (use_buffer) {
        plot(sf::st_geometry(site_proj_buffered), col = NA, border = "black", 
             add = TRUE, lty = 2, lwd = 1)
      }
      
      if (has_exclusions) {
        plot(sf::st_geometry(excl_proj), col = "lightgray", border = "black", 
             add = TRUE, lty = 3, lwd = 0.5)
      }
      
      points(sf::st_coordinates(pts_proj), col = "red", pch = 20, cex = 0.5)
      
      # Build legend dynamically
      legend_items <- c("Original site")
      legend_cols <- c("black")
      legend_lty <- c(1)
      legend_lwd <- c(2)
      
      if (use_buffer) {
        legend_items <- c(legend_items, "Buffered site (for analysis)")
        legend_cols <- c(legend_cols, "black")
        legend_lty <- c(legend_lty, 2)
        legend_lwd <- c(legend_lwd, 1)
      }
      
      if (has_exclusions) {
        legend_items <- c(legend_items, "Excluded areas")
        legend_cols <- c(legend_cols, "black")
        legend_lty <- c(legend_lty, 3)
        legend_lwd <- c(legend_lwd, 0.5)
      }
      
      legend("topright", legend = legend_items, col = legend_cols, 
             lty = legend_lty, lwd = legend_lwd, cex = 0.8)
    }
  }
  
  # ============================================================================
  # SECTION 6: OPTIMIZED POINT FILTERING AND BUFFER CREATION
  # ============================================================================
  
  msg("Filtering records and creating buffers (optimized)...")
  
  # Get maximum uncertainty to create search buffer around site
  max_uncertainty <- max(pts_proj$uncertainty)
  
  # Create expanded search area
  site_search_area <- sf::st_buffer(site_proj_processing, dist = max_uncertainty)
  
  # Pre-filter: keep only points whose locations could potentially intersect site
  potential_pts_idx <- sf::st_intersects(pts_proj, site_search_area, sparse = FALSE)[, 1]
  pts_filtered <- pts_proj[potential_pts_idx, ]
  
  msg(paste("Pre-filtered to", nrow(pts_filtered), "potentially relevant points..."))
  
  # Create buffers only for filtered points (memory optimization)
  buffers_filtered_sf <- sf::st_buffer(pts_filtered, dist = pts_filtered$uncertainty)
  
  # Check which buffers actually intersect the processing site
  intersects_idx <- sf::st_intersects(buffers_filtered_sf, site_proj_processing, sparse = FALSE)[, 1]
  
  # Keep only points and buffers that intersect
  pts_computed <- pts_filtered[intersects_idx, ]
  buffers_computed_sf <- buffers_filtered_sf[intersects_idx, ]
  
  # Stop if no records remain
  if (nrow(pts_computed) == 0) stop("No occurrence buffers intersect the study area.")
  
  # Remove unsuitable areas from buffers if provided
  if (has_exclusions) {
    msg("Removing exclusion areas from occurrence buffers...")
    buffers_computed_sf <- sf::st_make_valid(
      sf::st_difference(buffers_computed_sf, excl_proj)
    )
  }
  
  msg(paste("Retained", nrow(pts_computed), "records for computation."))
  
  # ============================================================================
  # SECTION 7: RASTER TEMPLATE CREATION WITH GUARANTEED COVERAGE
  # ============================================================================
  
  msg("Creating template raster with guaranteed boundary coverage...")
  
  # Get bounding box from processing site
  site_bbox <- sf::st_bbox(site_proj_processing)
  
  # Expand bbox by cellsize to ensure complete coverage of all corners/edges
  # This guarantees at least one cell center will be inside even small protrusions
  bbox_expanded <- site_bbox
  bbox_expanded["xmin"] <- site_bbox["xmin"] - cellsize
  bbox_expanded["xmax"] <- site_bbox["xmax"] + cellsize
  bbox_expanded["ymin"] <- site_bbox["ymin"] - cellsize
  bbox_expanded["ymax"] <- site_bbox["ymax"] + cellsize
  
  # Create empty raster template with specified resolution
  r_template <- terra::rast(
    extent = terra::ext(bbox_expanded),
    resolution = cellsize,
    crs = crs_terra
  )
  terra::values(r_template) <- NA
  
  msg(paste0("  Template extent expanded by ", cellsize, 
             "m on all sides to ensure complete coverage."))
  
  # ============================================================================
  # SECTION 8: SPECIES RICHNESS MAP CALCULATION
  # ============================================================================
  
  msg("Calculating species richness map (RICH)...")
  
  # Convert points to terra vector format
  pts_vect <- terra::vect(pts_computed)
  
  # Rasterize directly counting unique taxa per cell
  r_rich <- terra::rasterize(
    pts_vect,
    r_template,
    field = "Taxon",
    fun = function(x) length(unique(x))
  )
  
  # Set NA cells to 0 (cells with no observations)
  r_rich[is.na(r_rich)] <- 0
  
  # ============================================================================
  # SECTION 9: SPATIO-TEMPORAL SCORE CALCULATION
  # ============================================================================
  
  msg("Calculating spatio-temporal scores...")
  
  # Convert buffers to terra vector format
  v_buff <- terra::vect(buffers_computed_sf)
  
  # Extract cells covered by each buffer
  cell_data <- terra::extract(r_template, v_buff, cells = TRUE, ID = TRUE)
  
  # Count how many cells each buffer covers (for spatial uncertainty)
  counts_per_id <- table(cell_data$ID)
  
  # Initialize spatial count vector
  spatial_count <- rep(1, nrow(pts_computed))
  spatial_count[as.integer(names(counts_per_id))] <- as.numeric(counts_per_id)
  
  # Calculate spatial score (inverse of cells covered)
  pts_computed$spatial_score <- 1 / spatial_count
  
  # Calculate temporal score using exponential decay
  pts_computed$time_score <- (1 - (tau / 100))^((year_study - pts_computed$year) / 100)
  
  # Combine spatial and temporal scores
  pts_computed$st_ignorance <- pts_computed$spatial_score * pts_computed$time_score
  
  # ============================================================================
  # SECTION 10: PER-TAXON IGNORANCE RASTERIZATION
  # ============================================================================
  
  # Get unique taxa list
  taxa_list <- unique(pts_computed$Taxon)
  
  if (use_coverage_weighting) {
    # METHOD 1: Coverage-weighted (ACCURATE, matches original algorithm)
    msg("Drafting Map of Relative Floristic Ignorance (processing by taxon)...")
    msg("  Using coverage-weighted rasterization (accurate, slower)...")
    
    # Initialize progress bar
    pb <- utils::txtProgressBar(min = 0, max = length(taxa_list), style = 3)
    
    raster_list <- suppressWarnings(
      lapply(seq_along(taxa_list), function(i) {
        tname <- taxa_list[i]
        
        # Filter for this taxon
        taxon_idx <- pts_computed$Taxon == tname
        pts_taxon <- pts_computed[taxon_idx, ]
        bufs_taxon_sf <- buffers_computed_sf[taxon_idx, ]
        
        # Initialize empty raster for this taxon
        tax_r <- r_template
        terra::values(tax_r) <- 0
        
        # Process each buffer for this taxon individually
        for (j in 1:nrow(pts_taxon)) {
          # Create single buffer
          single_buf <- sf::st_sf(
            geometry = sf::st_geometry(bufs_taxon_sf[j, ])
          )
          
          # Rasterize to get coverage fraction (0-1)
          coverage_r <- terra::rasterize(
            terra::vect(single_buf),
            r_template,
            cover = TRUE
          )
          
          # Replace NA with 0
          coverage_r[is.na(coverage_r)] <- 0
          
          # Weight the ignorance value by coverage fraction
          weighted_r <- coverage_r * pts_taxon$st_ignorance[j]
          
          # Take maximum value per cell (if multiple buffers of same taxon overlap)
          tax_r <- max(tax_r, weighted_r, na.rm = TRUE)
        }
        
        utils::setTxtProgressBar(pb, i)
        return(tax_r)
      })
    )
    
    close(pb)
    
  } else {
    # METHOD 2: Binary touch method (FAST, but less accurate)
    msg("Drafting Map of Relative Floristic Ignorance (processing by taxon)...")
    msg("  Using binary touch method (fast, no coverage weighting)...")
    
    # Initialize progress bar
    pb <- utils::txtProgressBar(min = 0, max = length(taxa_list), style = 3)
    
    raster_list <- suppressWarnings(
      lapply(seq_along(taxa_list), function(i) {
        tname <- taxa_list[i]
        
        # Filter for this taxon
        taxon_idx <- pts_computed$Taxon == tname
        pts_taxon <- pts_computed[taxon_idx, ]
        bufs_taxon_sf <- buffers_computed_sf[taxon_idx, ]
        
        # Create clean sf object
        bufs_taxon_clean <- sf::st_sf(
          st_ignorance = pts_taxon$st_ignorance,
          geometry = sf::st_geometry(bufs_taxon_sf)
        )
        
        # Rasterize without coverage weighting (faster)
        tax_r <- terra::rasterize(
          terra::vect(bufs_taxon_clean),
          r_template,
          field = "st_ignorance",
          fun = "max",
          touches = TRUE
        )
        
        tax_r[is.na(tax_r)] <- 0
        
        utils::setTxtProgressBar(pb, i)
        return(tax_r)
      })
    )
    
    close(pb)
  }
  
  # ============================================================================
  # SECTION 11: MRFI RASTER FINALIZATION
  # ============================================================================
  
  msg("Finalizing rasters...")
  
  # Keep only valid SpatRaster objects
  raster_list_valid <- raster_list[sapply(raster_list, function(x) inherits(x, "SpatRaster"))]
  
  # Sum ignorance scores across all taxa
  if (length(raster_list_valid) > 0) {
    raster_stack <- terra::rast(raster_list_valid)
    raster_sum <- terra::app(raster_stack, fun = sum, na.rm = TRUE)
  } else {
    raster_sum <- r_template
    terra::values(raster_sum) <- 0
  }
  
  # Rescale to MRFI: maximum possible ignorance minus actual ignorance
  rmax <- max(terra::values(raster_sum), na.rm = TRUE)
  if (!is.finite(rmax)) rmax <- 0
  mrfi_r <- rmax - raster_sum
  
  # ============================================================================
  # SECTION 12: MASKING APPLICATION
  # ============================================================================
  
  msg(paste0("Applying mask with method: '", mask_method, "'..."))
  
  if (mask_method == "touches") {
    # Include any cell touched by the polygon (MRFI-aligned, inclusive)
    site_mask_r <- terra::rasterize(
      site_vect_mask,
      r_template,
      field = 1,
      touches = TRUE
    )
    site_mask_r[site_mask_r == 0] <- NA
    
    # Apply mask
    mrfi_final <- terra::mask(terra::crop(mrfi_r, site_mask_r), site_mask_r)
    rich_final <- terra::mask(terra::crop(r_rich, site_mask_r), site_mask_r)
    
  } else if (mask_method == "none") {
    # No masking - show entire extent (most MRFI-aligned)
    msg("  No masking applied - showing entire raster extent.")
    mrfi_final <- mrfi_r
    rich_final <- r_rich
  }
  
  # ============================================================================
  # SECTION 13: EXCLUSION AREA REMOVAL
  # ============================================================================
  
  if (has_exclusions) {
    msg("Removing excluded areas from final output...")
    
    # Rasterize exclusion areas (binary: 1 = excluded, NA = not excluded)
    excl_mask_r <- terra::rasterize(
      terra::vect(excl_proj),
      r_template,
      field = 1
    )
    
    # Set any cell overlapping exclusion areas to NA (binary exclusion)
    mrfi_final[!is.na(excl_mask_r)] <- NA
    rich_final[!is.na(excl_mask_r)] <- NA
  }
  
  # ============================================================================
  # SECTION 14: STATISTICS COMPILATION
  # ============================================================================
  
  msg("Compiling statistics...")
  end_time <- Sys.time()
  
  # Calculate excluded records count
  total_used_records <- nrow(pts_computed)
  excluded_records_count <- total_initial_records - total_used_records
  
  # Create statistics summary
  statistics_df <- data.frame(
    Statistic = c(
      "Started",
      "Finished",
      "Elapsed time (secs)",
      "Site buffer applied (m)",
      "Processing & masking boundary",
      "Final masking method",
      "Coverage weighting",
      "Color scale",
      if (color_scale == "quantile") "Number of quantiles" else NULL,
      "Exclusion areas applied",
      "CRS (EPSG code)",
      "Cell size (m)",
      "100 years % loss ratio (tau)",
      "Total initial records",
      "Total occurrences within original site",
      "Total occurrences computed (buffers)",
      "Records excluded from analysis",
      "Occ. uncertainty (median, m)",
      "Occ. dates (median, year)"
    ),
    Value = c(
      as.character(start_time),
      as.character(end_time),
      round(as.numeric(end_time - start_time, units = "secs")),
      if (use_buffer) round(buffer_distance, 1) else "None",
      if (use_buffer) paste0("Original site + ", round(buffer_distance, 1), "m buffer") else "Original site only",
      mask_method,
      if (use_coverage_weighting) "Enabled (accurate)" else "Disabled (fast)",
      color_scale,
      if (color_scale == "quantile") n_quantiles else NULL,
      if (has_exclusions) "Yes" else "No",
      as.character(CRS.new),
      cellsize,
      tau,
      total_initial_records,
      nrow(pts_in_site),
      nrow(pts_computed),
      excluded_records_count,
      round(median(pts_computed$uncertainty, na.rm = TRUE)),
      round(median(pts_computed$year, na.rm = TRUE))
    )
  )
  
  # ============================================================================
  # SECTION 15: PLOT GENERATION
  # ============================================================================
  
  msg("Generating plots and output files...")
  
  # Clip exclusion areas to site vicinity for plotting
  excl_plot <- NULL
  if (has_exclusions) {
    excl_plot <- sf::st_intersection(excl_proj, site_proj_processing)
    # Handle case where intersection might be empty
    if (length(excl_plot) == 0 || all(sf::st_is_empty(excl_plot))) {
      excl_plot <- NULL
    }
  }
  
  # Convert rasters to data frames for ggplot
  mrfi_df <- as.data.frame(mrfi_final, xy = TRUE)
  colnames(mrfi_df) <- c("x", "y", "value")
  
  rich_df <- as.data.frame(rich_final, xy = TRUE)
  colnames(rich_df) <- c("x", "y", "value")
  
  # Get max values
  mrfi_max_val <- terra::global(mrfi_final, "max", na.rm = TRUE)$max
  rich_max_val <- terra::global(rich_final, "max", na.rm = TRUE)$max
  
  # Create breaks based on color_scale setting
  if (color_scale == "quantile") {
    # Quantile breaks: highlights sampling differences
    msg("  Creating quantile breaks for MRFI plot...")
    
    # Get non-NA values from raster
    mrfi_values <- terra::values(mrfi_final, na.rm = TRUE)
    
    # Check if we have enough unique values for quantiles
    n_unique <- length(unique(mrfi_values))
    
    if (n_unique < n_quantiles) {
      warning("Not enough unique IRFI values (", n_unique, ") for ", n_quantiles, 
              " quantiles. Using ", n_unique, " breaks instead.")
      n_quantiles_actual <- n_unique
    } else {
      n_quantiles_actual <- n_quantiles
    }
    
    # Calculate quantile breaks (boundary values for legend)
    if (n_unique > 1) {
      mrfi_breaks <- stats::quantile(mrfi_values, 
                                     probs = seq(0, 1, length.out = n_quantiles_actual + 1), 
                                     na.rm = TRUE)
      mrfi_breaks <- unique(mrfi_breaks)  # Remove duplicates
    } else {
      # Only one unique value
      mrfi_breaks <- c(min(mrfi_values), max(mrfi_values))
    }
    
  } else {
    # Continuous breaks: evenly spaced (9 values for consistency)
    mrfi_breaks <- seq(0, mrfi_max_val, length.out = 9)
  }
  
  # Richness breaks: 9 values (consistent with MRFI continuous)
  rich_breaks <- seq(0, rich_max_val, length.out = 9)
  
  # --------------------------------------------------------------------------
  # Plot 1: Map of Relative Floristic Ignorance
  # --------------------------------------------------------------------------
  
  p1 <- ggplot2::ggplot(mrfi_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "IRFI") +
    ggplot2::theme(
      legend.position = "right",
      legend.direction = 'vertical',
      legend.key.width = grid::unit(0.6, "cm")
    ) +
    ggplot2::xlab("Longitude") +
    ggplot2::ylab("Latitude")
  
  # Add appropriate color scale based on method
  if (color_scale == "quantile") {
    # Quantile-based color scale with 9 colors
    p1 <- p1 + ggplot2::scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
      values = scales::rescale(mrfi_breaks),
      breaks = mrfi_breaks,
      labels = round(mrfi_breaks, 0),
      limits = c(min(mrfi_breaks), max(mrfi_breaks)),
      na.value = "transparent",
      guide = ggplot2::guide_legend(
        title = "IRFI",
        keyheight = grid::unit(1.2, "lines"),
        keywidth = grid::unit(1.2, "lines"),
        label.theme = ggplot2::element_text(size = 9),
        title.theme = ggplot2::element_text(size = 10, face = "bold"),
        override.aes = list(alpha = 1)
      )
    )
  } else {
    # Continuous color scale with 9 colors
    p1 <- p1 + ggplot2::scale_fill_distiller(
      palette = "Spectral",
      direction = -1,
      limits = c(min(mrfi_breaks), max(mrfi_breaks)),
      breaks = mrfi_breaks,
      labels = round(mrfi_breaks, 0)
    )
  }
  
  # Add raster tiles
  p1 <- p1 +
    ggplot2::geom_tile(
      mapping = ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
      alpha = 0.8,
      colour = "black",
      linewidth = 0.1
    )
  
  # Add original site boundary (solid black, thick)
  p1 <- p1 + ggplot2::geom_sf(
    data = site_proj_original,
    fill = NA,
    color = "black",
    linewidth = 1,
    inherit.aes = FALSE
  )
  
  # Add buffered boundary if buffer is used (dashed black, same thickness)
  if (use_buffer && !is.null(site_proj_buffered)) {
    p1 <- p1 + ggplot2::geom_sf(
      data = site_proj_buffered,
      fill = NA,
      color = "black",
      linewidth = 1,
      linetype = "dashed",
      inherit.aes = FALSE
    )
  }
  
  # Add exclusion areas (dotted, thin)
  if (has_exclusions && !is.null(excl_plot)) {
    p1 <- p1 + ggplot2::geom_sf(
      data = excl_plot,
      fill = NA,
      color = "black",
      linewidth = 0.5,
      linetype = "dotted",
      inherit.aes = FALSE
    )
  }
  
  # Add title
  p1 <- p1 + ggplot2::ggtitle(paste0("Map of Relative Floristic Ignorance (MRFI) - ", 
                                     ifelse(color_scale == "quantile", "Quantile scale", "Continuous scale")))
  
  # --------------------------------------------------------------------------
  # Plot 2: Species Richness Map
  # --------------------------------------------------------------------------
  
  p2 <- ggplot2::ggplot(rich_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      legend.direction = 'vertical',
      legend.key.width = grid::unit(0.6, "cm")
    ) +
    ggplot2::geom_tile(
      mapping = ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
      alpha = 0.8,
      colour = "black",
      linewidth = 0.1
    ) +
    ggplot2::scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(9, "Spectral"),
      limits = c(0, rich_max_val),
      breaks = rich_breaks,
      labels = round(rich_breaks, 0),
      guide = ggplot2::guide_legend(title = "Value")
    ) +
    ggplot2::ggtitle("Species richness map (without uncertainties)") +
    ggplot2::xlab("Longitude") +
    ggplot2::ylab("Latitude")
  
  # Add original site boundary (solid black, thick)
  p2 <- p2 + ggplot2::geom_sf(
    data = site_proj_original,
    fill = NA,
    color = "black",
    linewidth = 1,
    inherit.aes = FALSE
  )
  
  # Add buffered boundary if buffer is used (dashed black, same thickness)
  if (use_buffer && !is.null(site_proj_buffered)) {
    p2 <- p2 + ggplot2::geom_sf(
      data = site_proj_buffered,
      fill = NA,
      color = "black",
      linewidth = 1,
      linetype = "dashed",
      inherit.aes = FALSE
    )
  }
  
  # Add exclusion areas (dotted, thin)
  if (has_exclusions && !is.null(excl_plot)) {
    p2 <- p2 + ggplot2::geom_sf(
      data = excl_plot,
      fill = NA,
      color = "black",
      linewidth = 0.5,
      linetype = "dotted",
      inherit.aes = FALSE
    )
  }
  
  # --------------------------------------------------------------------------
  # Plot 3: Temporal Uncertainty (Occurrence Year Distribution)
  # --------------------------------------------------------------------------
  
  p3 <- ggplot2::ggplot(pts_computed) +
    ggplot2::aes(x = .data$year, y = ggplot2::after_stat(density)) +
    ggplot2::geom_histogram(
      alpha = 0.6,
      fill = "#FF6666",
      binwidth = diff(range(pts_computed$year)) / 30
    ) +
    ggplot2::coord_cartesian(xlim = c(min(pts_computed$year), year_study)) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    ggplot2::ggtitle("Occurrence date") +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Frequency") +
    ggplot2::theme_classic()
  
  # --------------------------------------------------------------------------
  # Plot 4: Spatial Uncertainty Distribution
  # --------------------------------------------------------------------------
  
  p4 <- ggplot2::ggplot(pts_computed) +
    ggplot2::aes(x = .data$uncertainty, y = ggplot2::after_stat(density)) +
    ggplot2::geom_histogram(
      alpha = 0.6,
      fill = "#FF6666",
      binwidth = diff(range(pts_computed$uncertainty)) / 30
    ) +
    ggplot2::coord_cartesian(
      xlim = c(min(pts_computed$uncertainty), max(pts_computed$uncertainty))
    ) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    ggplot2::ggtitle("Occurrence spatial uncertainty") +
    ggplot2::xlab("Uncertainty (m)") +
    ggplot2::ylab("Frequency") +
    ggplot2::theme_classic()
  
  # ============================================================================
  # SECTION 16: FILE OUTPUT
  # ============================================================================
  
  # Create file paths with custom prefix
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_output.pdf"))
  tif_path <- file.path(output_dir, paste0(output_prefix, "_map.tif"))
  csv_path <- file.path(output_dir, paste0(output_prefix, "_taxa.csv"))
  
  # Save multi-page PDF with all plots and statistics
  grDevices::pdf(pdf_path, onefile = TRUE)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  grid::grid.draw(
    gridExtra::grid.arrange(
      top = "Summary statistics",
      gridExtra::tableGrob(statistics_df)
    )
  )
  grDevices::dev.off()
  
  # Save MRFI raster as GeoTIFF
  terra::writeRaster(mrfi_final, filename = tif_path, overwrite = TRUE)
  
  # Save list of taxa considered
  utils::write.csv(taxa_list, row.names = FALSE, csv_path)
  
  msg(paste0("Done! Files saved to: ", output_dir))
  
  # ============================================================================
  # SECTION 17: CONSOLE DISPLAY
  # ============================================================================
  
  # Display plots in console
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  grid::grid.draw(
    gridExtra::grid.arrange(
      top = "Summary statistics",
      gridExtra::tableGrob(statistics_df)
    )
  )
  
  # ============================================================================
  # SECTION 18: RETURN RESULTS
  # ============================================================================
  
  return(list(
    MRFI = mrfi_final,
    RICH = rich_final,
    Uncertainties = data.frame(
      uncertainty = pts_computed$uncertainty,
      year = pts_computed$year,
      Taxon = pts_computed$Taxon
    ),
    Statistics = statistics_df
  ))
}