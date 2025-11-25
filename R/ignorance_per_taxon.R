#' @title Calculate Per-Taxon Ignorance Rasters (Internal)
#'
#' @description
#' Internal function that calculates ignorance contribution for each taxon
#' separately. Used by compare_mrfi() for incremental MRFI updates.
#' Results are saved to temporary files to minimize memory usage.
#'
#' @param data_flor Data frame with floristic data (Taxon, Long, Lat, uncertainty, year)
#' @param site An `sf` object representing the study area
#' @param year_study Numeric year for temporal score calculation
#' @param excl_areas Optional `sf` object for exclusion areas
#' @param CRS.new Numeric EPSG code for projected CRS (default = 3035)
#' @param tau Numeric. Temporal decay parameter (% taxa loss per 100 years)
#' @param cellsize Numeric. Raster cell size in meters
#' @param template_raster Optional SpatRaster to use as template (for consistency)
#' @param use_coverage_weighting Logical. Use coverage-weighted rasterization (default = TRUE)
#' @param temp_dir Character. Directory for temporary files (default = tempdir())
#' @param keep_in_memory Logical. If FALSE (default), save to temp files; if TRUE, keep in memory
#' @param verbose Logical. Print progress messages (default = TRUE)
#' @param site_buffer Logical. Expand site boundary (default = FALSE)
#' @param buffer_width Numeric. Buffer width in meters if site_buffer = TRUE
#' @param mask_method Character. Masking method: "touches" or "none" (default = "touches")
#'
#' @return A list with:
#' \itemize{
#'   \item{taxon_files}{ Character vector of paths to temporary raster files}
#'   \item{taxon_names}{ Character vector of taxon names (matches file order)}
#'   \item{taxon_observations}{ sf object with observation data and calculated scores}
#'   \item{baseline_year}{ Year used for calculation}
#'   \item{template_raster}{ SpatRaster template used}
#'   \item{parameters}{ List of calculation parameters}
#'   \item{temp_dir}{ Directory containing temporary files}
#' }
#'
#' @importFrom sf st_as_sf st_make_valid st_crs st_transform st_buffer st_intersects st_intersection st_union st_bbox st_coordinates st_geometry st_difference
#' @importFrom terra rast values mask crop vect rasterize extract ext writeRaster
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @keywords internal
#' @noRd

ignorance_per_taxon <- function(
    data_flor,
    site,
    year_study,
    excl_areas = NULL,
    CRS.new = 3035,
    tau,
    cellsize,
    template_raster = NULL,
    use_coverage_weighting = TRUE,
    temp_dir = tempdir(),
    keep_in_memory = FALSE,
    verbose = TRUE,
    site_buffer = FALSE,
    buffer_width = NULL,
    mask_method = "touches"
) {
  
  msg <- function(...) if (verbose) message(...)
  
  msg("  [PER-TAXON IGNORANCE CALCULATION]")
  
  # ==========================================================================
  # SECTION 1: SETUP AND VALIDATION (simplified from ignorance_map)
  # ==========================================================================
  
  if (is.null(year_study)) year_study <- as.numeric(format(Sys.Date(), "%Y"))
  
  # Validate inputs
  req_cols <- c("Taxon", "Long", "Lat", "uncertainty", "year")
  if (!all(req_cols %in% names(data_flor))) {
    stop("data_flor must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  if (any(is.na(data_flor[, c("uncertainty", "Lat", "Long", "year")]))) {
    stop("Missing values (NA) found in critical columns")
  }
  
  # ==========================================================================
  # SECTION 2: COORDINATE SYSTEM SETUP
  # ==========================================================================
  
  crs_sf <- sf::st_crs(CRS.new)
  crs_terra <- paste0("EPSG:", CRS.new)
  
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  if (is.na(sf::st_crs(site))) sf::st_crs(site) <- 4326
  
  site_proj_original <- sf::st_transform(sf::st_make_valid(site), crs_sf)
  
  # Handle exclusions
  excl_proj <- NULL
  has_exclusions <- !is.null(excl_areas)
  
  if (has_exclusions) {
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    if (is.na(sf::st_crs(excl_areas))) sf::st_crs(excl_areas) <- 4326
    excl_proj <- sf::st_union(sf::st_transform(sf::st_make_valid(excl_areas), crs_sf))
  }
  
  # ==========================================================================
  # SECTION 3: SITE BUFFER LOGIC
  # ==========================================================================
  
  if (site_buffer) {
    if (is.null(buffer_width)) buffer_width <- cellsize
    site_proj_buffered <- sf::st_buffer(site_proj_original, dist = buffer_width)
    
    if (has_exclusions) {
      site_proj_buffered <- sf::st_difference(site_proj_buffered, excl_proj)
    }
    
    site_proj_processing <- sf::st_make_valid(site_proj_buffered)
  } else {
    site_proj_processing <- site_proj_original
  }
  
  site_vect_processing <- terra::vect(site_proj_processing)
  
  # ==========================================================================
  # SECTION 4: FLORISTIC DATA PREPARATION
  # ==========================================================================
  
  if (inherits(data_flor, "sf")) {
    pts_sf <- data_flor
  } else if (inherits(data_flor, "Spatial")) {
    pts_sf <- sf::st_as_sf(data_flor)
  } else if (is.data.frame(data_flor)) {
    pts_sf <- sf::st_as_sf(data_flor, coords = c("Long", "Lat"), 
                           crs = 4326, remove = FALSE)
  } else {
    stop("Unsupported data_flor type")
  }
  
  pts_proj <- sf::st_transform(pts_sf, crs_sf)
  
  # Remove records with NA uncertainty
  if (any(is.na(pts_proj$uncertainty))) {
    n_removed <- sum(is.na(pts_proj$uncertainty))
    msg(paste("  Removing", n_removed, "records with NA uncertainty..."))
    pts_proj <- pts_proj[!is.na(pts_proj$uncertainty), ]
    
    if (nrow(pts_proj) == 0) {
      stop("No valid records remaining after removing NA uncertainties")
    }
  }
  
  # ==========================================================================
  # SECTION 5: OPTIMIZED POINT FILTERING
  # ==========================================================================
  
  max_uncertainty <- max(pts_proj$uncertainty)
  site_search_area <- sf::st_buffer(site_proj_processing, dist = max_uncertainty)
  
  potential_pts_idx <- sf::st_intersects(pts_proj, site_search_area, sparse = FALSE)[, 1]
  pts_filtered <- pts_proj[potential_pts_idx, ]
  
  buffers_filtered_sf <- sf::st_buffer(pts_filtered, dist = pts_filtered$uncertainty)
  intersects_idx <- sf::st_intersects(buffers_filtered_sf, site_proj_processing, 
                                      sparse = FALSE)[, 1]
  
  pts_computed <- pts_filtered[intersects_idx, ]
  buffers_computed_sf <- buffers_filtered_sf[intersects_idx, ]
  
  if (nrow(pts_computed) == 0) {
    stop("No occurrence buffers intersect the study area")
  }
  
  if (has_exclusions) {
    buffers_computed_sf <- sf::st_make_valid(
      sf::st_difference(buffers_computed_sf, excl_proj)
    )
  }
  
  msg(paste("  Processing", nrow(pts_computed), "records..."))
  
  # ==========================================================================
  # SECTION 6: RASTER TEMPLATE
  # ==========================================================================
  
  if (is.null(template_raster)) {
    site_bbox <- sf::st_bbox(site_proj_processing)
    bbox_expanded <- site_bbox
    bbox_expanded["xmin"] <- site_bbox["xmin"] - cellsize
    bbox_expanded["xmax"] <- site_bbox["xmax"] + cellsize
    bbox_expanded["ymin"] <- site_bbox["ymin"] - cellsize
    bbox_expanded["ymax"] <- site_bbox["ymax"] + cellsize
    
    template_raster <- terra::rast(
      extent = terra::ext(bbox_expanded),
      resolution = cellsize,
      crs = crs_terra
    )
    terra::values(template_raster) <- NA
  }
  
  # ==========================================================================
  # SECTION 7: SPATIO-TEMPORAL SCORE CALCULATION
  # ==========================================================================
  
  v_buff <- terra::vect(buffers_computed_sf)
  cell_data <- terra::extract(template_raster, v_buff, cells = TRUE, ID = TRUE)
  counts_per_id <- table(cell_data$ID)
  
  spatial_count <- rep(1, nrow(pts_computed))
  spatial_count[as.integer(names(counts_per_id))] <- as.numeric(counts_per_id)
  
  pts_computed$spatial_score <- 1 / spatial_count
  pts_computed$time_score <- (1 - (tau / 100))^((year_study - pts_computed$year) / 100)
  pts_computed$st_ignorance <- pts_computed$spatial_score * pts_computed$time_score
  
  # ==========================================================================
  # SECTION 8: PER-TAXON RASTERIZATION
  # ==========================================================================
  
  taxa_list <- unique(pts_computed$Taxon)
  n_taxa <- length(taxa_list)
  
  msg(paste("  Calculating ignorance for", n_taxa, "taxa..."))
  
  # Create output storage
  if (keep_in_memory) {
    taxon_rasters <- list()
    taxon_files <- NULL
  } else {
    # Ensure temp directory exists
    if (!dir.exists(temp_dir)) {
      dir.create(temp_dir, recursive = TRUE)
    }
    taxon_files <- character(n_taxa)
    taxon_rasters <- NULL
  }
  
  pb <- utils::txtProgressBar(min = 0, max = n_taxa, style = 3)
  
  for (i in seq_along(taxa_list)) {
    tname <- taxa_list[i]
    taxon_idx <- pts_computed$Taxon == tname
    pts_taxon <- pts_computed[taxon_idx, ]
    bufs_taxon_sf <- buffers_computed_sf[taxon_idx, ]
    
    # Calculate taxon raster
    if (use_coverage_weighting) {
      # Coverage-weighted approach
      tax_r <- template_raster
      terra::values(tax_r) <- 0
      
      for (j in 1:nrow(pts_taxon)) {
        single_buf <- sf::st_sf(geometry = sf::st_geometry(bufs_taxon_sf[j, ]))
        coverage_r <- terra::rasterize(terra::vect(single_buf), 
                                       template_raster, cover = TRUE)
        coverage_r[is.na(coverage_r)] <- 0
        weighted_r <- coverage_r * pts_taxon$st_ignorance[j]
        tax_r <- max(tax_r, weighted_r, na.rm = TRUE)
      }
    } else {
      # Binary touch approach
      bufs_taxon_clean <- sf::st_sf(
        st_ignorance = pts_taxon$st_ignorance,
        geometry = sf::st_geometry(bufs_taxon_sf)
      )
      
      tax_r <- terra::rasterize(
        terra::vect(bufs_taxon_clean),
        template_raster,
        field = "st_ignorance",
        fun = "max",
        touches = TRUE
      )
      tax_r[is.na(tax_r)] <- 0
    }
    
    # Store result
    if (keep_in_memory) {
      taxon_rasters[[tname]] <- tax_r
    } else {
      # Save to temporary file
      temp_file <- file.path(temp_dir, paste0("taxon_", i, "_", 
                                              gsub("[^A-Za-z0-9]", "_", tname), 
                                              ".tif"))
      terra::writeRaster(tax_r, temp_file, overwrite = TRUE)
      taxon_files[i] <- temp_file
    }
    
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  msg("  Per-taxon calculation complete")
  
  # ==========================================================================
  # SECTION 9: RETURN RESULTS
  # ==========================================================================
  
  return(list(
    taxon_files = taxon_files,
    taxon_rasters = taxon_rasters,  # NULL if using files
    taxon_names = taxa_list,
    taxon_observations = pts_computed,
    baseline_year = year_study,
    template_raster = template_raster,
    parameters = list(
      tau = tau,
      cellsize = cellsize,
      CRS.new = CRS.new,
      use_coverage_weighting = use_coverage_weighting,
      site_buffer = site_buffer,
      buffer_width = buffer_width,
      mask_method = mask_method
    ),
    temp_dir = temp_dir,
    keep_in_memory = keep_in_memory
  ))
}