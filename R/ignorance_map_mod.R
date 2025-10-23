#' @title Ignorance map (Modernized & Optimized)
#'
#' @description Computes a Map of Relative Floristic Ignorance (MRFI) using
#' modern `sf` and `terra` packages for high performance.
#'
#' @param data_flor A data.frame with 5 columns: 'Taxon', 'Long', 'Lat',
#' 'uncertainty' (radius in meters), and 'year'.
#' @param site An `sf` object (or 'SpatialPolygonsDataFrame') representing the study
#' area. Assumed to be in a geographic CRS (e.g., EPSG:4326) if not set.
#' @param year_study The numeric year of the study (e.g., 2025). Defaults to
#' the current system year.
#' @param excl_areas An optional `sf` object (or 'SpatialPolygonsDataFrame')
#' to delimit unsuitable areas to be excluded.
#' @param CRS.new The numeric EPSG code for the projected CRS to use for all
#' calculations (must be in meters). Default = 3035 (ETRS89-LAEA Europe).
#' @param tau Percentual value of taxa loss in 100 years (0 <= tau < 100).
#' @param cellsize The resolution of the output map (in meters).
#' @param verbose Logical. If TRUE, prints progress messages.
#' @param check_overlap Logical. If TRUE, checks and plots point-site overlap.
#' @param output_dir Directory for output files. Defaults to working directory.
#' @param output_prefix Prefix for output filenames. Default = "Ignorance".
#'
#' @return A list with 4 objects:
#' \itemize{
#'  \item{`MRFI`}{ A `terra SpatRaster` of the Map of Relative Floristic Ignorance.}
#'  \item{`RICH`}{ A `terra SpatRaster` of species richness, computed without
#'  uncertainties.}
#'  \item{`Uncertainties`}{ A data.frame of the uncertainty and year values for
#'  all records used in the computation.}
#'  \item{`Statistics`}{ A data.frame summarizing the settings and results.}
#' }
#'
#' @importFrom sf st_as_sf st_make_valid st_crs st_transform st_buffer st_intersects st_intersection st_union st_bbox st_coordinates st_geometry st_difference
#' @importFrom terra rast values ncell cellFromXY mask crop crs vect rasterize extract ext global writeRaster app
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom ggplot2 ggplot aes geom_tile geom_sf coord_equal theme_classic labs theme xlab ylab scale_fill_distiller guide_legend ggtitle guides geom_histogram coord_cartesian scale_y_continuous after_stat
#' @importFrom grid unit grid.draw
#' @importFrom gridExtra grid.arrange tableGrob
#' @importFrom grDevices pdf dev.off
#' @export
#' @examples \dontrun{
#' # Note: Requires the 'ignobioR' data objects
#' data(floratus)
#' data(park)
#' data(unsuitablezone)
#'
#' # Coerce old sp data to sf for modern use
#' park_sf <- sf::st_as_sf(park)
#' unsuitablezone_sf <- sf::st_as_sf(unsuitablezone)
#'
#' # Short example
#' set.seed(123)
#' mrfi <- ignorance_map_mod(
#'   data_flor = floratus[sample(nrow(floratus), 2000), ],
#'   site = park_sf,
#'   tau = 80,
#'   cellsize = 2000
#' )
#'
#' # Plot the MRFI raster
#' terra::plot(mrfi$MRFI)
#'
#' # Extended example
#' mrfi_ext <- ignorance_map_mod(
#'   data_flor = floratus,
#'   excl_areas = unsuitablezone_sf,
#'   site = park_sf,
#'   tau = 20,
#'   cellsize = 2000
#' )
#'
#' terra::plot(mrfi_ext$MRFI)
#' terra::plot(mrfi_ext$RICH)
#' }

ignorance_map_mod <- function(data_flor, site, year_study = NULL, excl_areas = NULL,
                              CRS.new = 3035, tau, cellsize, verbose = TRUE,
                              check_overlap = TRUE, output_dir = getwd(),
                              output_prefix = "Ignorance") {
  
  # Helper function for conditional messages
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  # --- 1. Settings and Input Validation ---
  msg("Checking settings and inputs...")
  
  # Use current year if not specified
  if (is.null(year_study)) year_study <- as.numeric(format(Sys.Date(), "%Y"))
  
  # Validate year_study range
  if (year_study < 1800 || year_study > 2100) {
    warning("year_study seems unusual. Double-check the value.")
  }
  
  # Check for future dates in occurrence data
  if (max(data_flor$year) > year_study) {
    warning("Some occurrence dates are more recent than year_study")
  }
  
  # Validate tau parameter (percentage loss rate)
  if (tau < 0 || tau >= 100) stop("0 <= tau < 100 is required.")
  
  # Validate cellsize parameter
  if (cellsize <= 0) stop("cellsize must be positive")
  
  # Validate CRS parameter
  if (CRS.new <= 0 || !is.numeric(CRS.new)) {
    stop("CRS.new must be a valid positive EPSG code")
  }
  
  # Check required column names in data_flor
  req_cols <- c("Taxon", "Long", "Lat", "uncertainty", "year")
  if (!all(req_cols %in% names(data_flor))) {
    stop("data_flor must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  # Store initial record count for statistics
  total_initial_records <- nrow(data_flor)
  
  # Warn if uncertainties are too small relative to cell size
  # This threshold ensures records are large enough to be captured by the grid
  if (any(2 * data_flor$uncertainty < (cellsize / 20))) {
    stop("Some records have uncertainty too small vs. cellsize. They may not be captured by the raster grid.")
  }
  
  msg("Inputs validated.")
  
  # --- 2. Coerce to SF and Reproject ---
  msg(paste("Reprojecting inputs to EPSG:", CRS.new))
  
  # Define CRS in both formats (used throughout)
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
  site_proj <- sf::st_transform(sf::st_make_valid(site), crs_sf)
  
  # SMART FIX: Create positive buffer to ensure complete raster coverage
  # Buffer outward by half a cell so any cell partially inside the site gets included
  msg("Creating expanded site boundary for complete raster coverage...")
  site_proj_expanded <- sf::st_buffer(site_proj, dist = cellsize / 2)
  
  # Validate the buffered geometry (in case buffer created issues)
  site_proj_expanded <- sf::st_make_valid(site_proj_expanded)
  
  # Convert both original and expanded site to terra vectors
  site_vect_expanded <- terra::vect(site_proj_expanded)  # Used for rasterization
  site_vect_original <- terra::vect(site_proj)           # Used for final masking
  
  # Handle exclusion areas if provided
  excl_proj <- NULL
  if (!is.null(excl_areas)) {
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
  
  # Convert occurrence points to sf if not already
  if (!inherits(data_flor, "sf")) {
    data_flor <- sf::st_as_sf(data_flor, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)
  }
  
  # Reproject points to target CRS
  pts_proj <- sf::st_transform(data_flor, crs_sf)
  
  # Calculate points within site for later statistics (using original site boundary)
  pts_in_site <- sf::st_intersection(pts_proj, site_proj)
  
  # --- 3. Optional: Check and visualize overlap ---
  if (check_overlap) {
    # Test spatial intersection between points and site (using original boundary)
    overlap_check <- sf::st_intersects(pts_proj, site_proj, sparse = FALSE)
    n_overlap <- sum(overlap_check)
    msg(paste("Number of occurrence points overlapping the site:", n_overlap))
    
    if (n_overlap == 0) {
      warning("No points overlap the study area. Check projections or coordinates.")
    } else {
      # Create diagnostic plot showing both original and expanded boundaries
      plot(sf::st_geometry(site_proj), col = NA, border = "black", main = "Points vs Site", lwd = 2)
      plot(sf::st_geometry(site_proj_expanded), col = NA, border = "blue", add = TRUE, lty = 2)
      points(sf::st_coordinates(pts_proj), col = "red", pch = 20, cex = 0.5)
      legend("topright", legend = c("Original boundary", "Expanded boundary"), 
             col = c("black", "blue"), lty = c(1, 2), lwd = c(2, 1))
    }
  }
  
  # --- 4. OPTIMIZED: Pre-filter Points Before Buffering ---
  msg("Filtering records and creating buffers (optimized)...")
  
  # Get maximum uncertainty to create search buffer around site
  max_uncertainty <- max(pts_proj$uncertainty)
  
  # Create expanded search area (using expanded site + max buffer distance)
  site_search_area <- sf::st_buffer(site_proj_expanded, dist = max_uncertainty)
  
  # Pre-filter: keep only points whose locations could potentially intersect site
  # This is MUCH faster than creating all buffers first
  potential_pts_idx <- sf::st_intersects(pts_proj, site_search_area, sparse = FALSE)[, 1]
  pts_filtered <- pts_proj[potential_pts_idx, ]
  
  msg(paste("Pre-filtered to", nrow(pts_filtered), "potentially relevant points"))
  
  # NOW create buffers only for filtered points (huge memory/time savings)
  buffers_filtered_sf <- sf::st_buffer(pts_filtered, dist = pts_filtered$uncertainty)
  
  # Check which buffers actually intersect the expanded site
  intersects_idx <- sf::st_intersects(buffers_filtered_sf, site_proj_expanded, sparse = FALSE)[, 1]
  
  # Keep only points and buffers that intersect
  pts_computed <- pts_filtered[intersects_idx, ]
  buffers_computed_sf <- buffers_filtered_sf[intersects_idx, ]
  
  # Stop if no records remain
  if (nrow(pts_computed) == 0) stop("No occurrence buffers intersect the study area.")
  
  # Remove unsuitable areas from buffers if provided and validate geometries
  if (!is.null(excl_proj)) {
    buffers_computed_sf <- sf::st_make_valid(
      sf::st_difference(buffers_computed_sf, excl_proj)
    )
  }
  
  msg(paste("Retained", nrow(pts_computed), "records for computation."))
  
  # --- 5. Create Template Raster (Based on expanded site extent) ---
  msg("Creating template raster...")
  
  # Get bounding box from expanded site
  site_bbox <- sf::st_bbox(site_proj_expanded)
  
  # Create empty raster template with specified resolution
  r_template <- terra::rast(
    extent = terra::ext(site_bbox),
    resolution = cellsize,
    crs = crs_terra
  )
  terra::values(r_template) <- NA
  
  # --- 6. OPTIMIZED: Calculate Richness Map Using terra::rasterize ---
  msg("Calculating species richness map (RICH)...")
  
  # Convert points to terra vector format with taxon attribute
  pts_vect <- terra::vect(pts_computed)
  
  # Rasterize directly counting unique taxa per cell (no manual aggregation needed)
  r_rich <- terra::rasterize(
    pts_vect, 
    r_template, 
    field = "Taxon", 
    fun = function(x) length(unique(x))
  )
  
  # Set NA cells to 0 (cells with no observations)
  r_rich[is.na(r_rich)] <- 0
  
  # --- 7. OPTIMIZED: Calculate Per-Record Spatio-Temporal Scores ---
  msg("Calculating spatio-temporal scores...")
  
  # Convert buffers to terra vector format
  v_buff <- terra::vect(buffers_computed_sf)
  
  # Extract cells covered by each buffer (ID corresponds to buffer index)
  # No need to create r_cells raster - extract already returns what we need
  cell_data <- terra::extract(r_template, v_buff, cells = TRUE, ID = TRUE)
  
  # Count how many cells each buffer covers (spatial uncertainty)
  counts_per_id <- table(cell_data$ID)
  
  # Initialize spatial count vector
  spatial_count <- rep(1, nrow(pts_computed))
  spatial_count[as.integer(names(counts_per_id))] <- as.numeric(counts_per_id)
  
  # Calculate spatial score (inverse of cells covered)
  pts_computed$spatial_score <- 1 / spatial_count
  
  # Calculate temporal score using exponential decay based on tau
  pts_computed$time_score <- (1 - (tau / 100))^((year_study - pts_computed$year) / 100)
  
  # Combine spatial and temporal scores
  pts_computed$st_ignorance <- pts_computed$spatial_score * pts_computed$time_score
  
  # --- 8. Rasterize Per-Taxon Ignorance with Progress Bar ---
  msg("Drafting Map of Relative Floristic Ignorance (processing by taxon)...")
  
  # Get unique taxa list
  taxa_list <- unique(pts_computed$Taxon)
  
  # Initialize progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(taxa_list), style = 3)
  
  # Process each taxon separately to get maximum ignorance per cell
  # Suppress sf attribute warnings during rasterization
  raster_list <- suppressWarnings(
    lapply(seq_along(taxa_list), function(i) {
      tname <- taxa_list[i]
      
      # Filter points and buffers for this taxon
      taxon_idx <- pts_computed$Taxon == tname
      pts_taxon <- pts_computed[taxon_idx, ]
      bufs_taxon_sf <- buffers_computed_sf[taxon_idx, ]
      
      # Create a clean sf object with only geometry and the score column
      # This avoids the "attribute variables assumed constant" warning
      bufs_taxon_clean <- sf::st_sf(
        st_ignorance = pts_taxon$st_ignorance,
        geometry = sf::st_geometry(bufs_taxon_sf)
      )
      
      # Rasterize buffers, taking maximum ignorance value per cell
      tax_r <- terra::rasterize(
        terra::vect(bufs_taxon_clean),
        r_template,
        field = "st_ignorance",
        fun = "max",
        touches = TRUE  # Include cells touched by buffer edges
      )
      
      # Convert NA to 0 (cells with no observations of this taxon)
      # This is necessary for summing across taxa layers
      tax_r[is.na(tax_r)] <- 0
      
      # Update progress bar
      utils::setTxtProgressBar(pb, i)
      return(tax_r)
    })
  )
  close(pb)
  
  # --- 9. Finalize MRFI Raster ---
  msg("Finalizing rasters...")
  
  # Keep only valid SpatRaster objects
  raster_list_valid <- raster_list[sapply(raster_list, function(x) inherits(x, "SpatRaster"))]
  
  # Sum ignorance scores across all taxa
  if (length(raster_list_valid) > 0) {
    raster_stack <- terra::rast(raster_list_valid)
    raster_sum <- terra::app(raster_stack, fun = sum, na.rm = TRUE)
  } else {
    # No valid rasters (shouldn't happen)
    raster_sum <- r_template
    terra::values(raster_sum) <- 0
  }
  
  # Rescale to MRFI: maximum possible ignorance minus actual ignorance
  rmax <- max(terra::values(raster_sum), na.rm = TRUE)
  if (!is.finite(rmax)) rmax <- 0
  mrfi_r <- rmax - raster_sum
  
  # Mask rasters using a site mask with touches = TRUE (Option A) ---
  msg("Applying mask: keeping all cells touched by site polygon...")
  
  # Create a mask raster from the site polygon: cells touched by site polygon are set to 1
  site_mask_r <- terra::rasterize(
    site_vect_original,
    r_template,
    field = 1,
    touches = TRUE
  )
  
  # Convert zeros to NA so only intersecting cells are retained
  site_mask_r[site_mask_r == 0] <- NA
  
  # Apply mask (crop first for performance)
  mrfi_final <- terra::mask(terra::crop(mrfi_r, site_mask_r), site_mask_r)
  rich_final <- terra::mask(terra::crop(r_rich, site_mask_r), site_mask_r)
  
  # --- 10. Compile Statistics ---
  msg("Compiling statistics...")
  end_time <- Sys.time()
  
  # Calculate excluded records count
  total_used_records <- nrow(pts_computed)
  excluded_records_count <- total_initial_records - total_used_records
  
  # Create statistics summary table
  statistics_df <- data.frame(
    Statistic = c(
      "Started", 
      "Finished", 
      "Elapsed time (secs)", 
      "CRS (EPSG code)",
      "Cell size (m)", 
      "100 years % loss ratio (tau)",
      "Total initial records",
      "Total occurrences within site (points)",
      "Total occurrences computed (buffers)",
      "Records excluded from analysis",
      "Occ. uncertainty (median, m)",
      "Occ. dates (median, year)"
    ),
    Value = c(
      as.character(start_time),
      as.character(end_time),
      round(as.numeric(end_time - start_time, units = "secs")),
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
  
  # --- 11. Generate Plots ---
  msg("Generating plots and output files...")
  
  # Convert rasters to data frames for ggplot
  mrfi_df <- as.data.frame(mrfi_final, xy = TRUE)
  colnames(mrfi_df) <- c("x", "y", "value")
  
  rich_df <- as.data.frame(rich_final, xy = TRUE)
  colnames(rich_df) <- c("x", "y", "value")
  
  # Get max values for legend breaks
  mrfi_max_val <- terra::global(mrfi_final, "max", na.rm = TRUE)$max
  rich_max_val <- terra::global(rich_final, "max", na.rm = TRUE)$max
  
  # Create evenly spaced breaks for legends
  mrfi_breaks <- seq(0, mrfi_max_val, length.out = 12)
  rich_breaks <- seq(0, rich_max_val, length.out = 12)
  
  # Plot 1: Map of Relative Floristic Ignorance
  p1 <- ggplot2::ggplot(mrfi_df) +
    ggplot2::coord_equal() + 
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "IFI") +
    ggplot2::theme(
      legend.position = "right", 
      legend.direction = 'vertical', 
      legend.key.width = grid::unit(0.6, "cm")
    ) +
    ggplot2::xlab("Longitude") + 
    ggplot2::ylab("Latitude") +
    ggplot2::scale_fill_distiller(
      palette = "Spectral", 
      direction = -1, 
      guide = ggplot2::guide_legend(),
      limits = c(0, mrfi_max_val),
      breaks = mrfi_breaks,
      labels = round(mrfi_breaks, 0) 
    ) +
    ggplot2::geom_tile(
      mapping = ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value), 
      alpha = 0.8, 
      colour = "black", 
      linewidth = 0.1
    ) +
    ggplot2::geom_sf(
      data = site_proj,  # Use original boundary for display
      fill = NA, 
      color = "black", 
      size = 1, 
      inherit.aes = FALSE
    ) +
    ggplot2::ggtitle("Map of Relative Floristic Ignorance (MRFI)")
  
  # Plot 2: Species Richness Map
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
    ggplot2::geom_sf(
      data = site_proj,  # Use original boundary for display
      fill = NA, 
      color = "black", 
      size = 1, 
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_distiller(
      palette = "Spectral", 
      direction = +1, 
      guide = ggplot2::guide_legend(title = "Value"),
      limits = c(0, rich_max_val),
      breaks = rich_breaks,
      labels = round(rich_breaks, 0) 
    ) +
    ggplot2::ggtitle("Species richness map (without uncertainties)") +
    ggplot2::xlab("Longitude") + 
    ggplot2::ylab("Latitude") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Value"))
  
  # Plot 3: Temporal Uncertainty (Occurrence Year Distribution)
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
  
  # Plot 4: Spatial Uncertainty Distribution
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
  
  # --- 12. Save Output Files ---
  # Create file paths with custom prefix
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_output.pdf"))
  tif_path <- file.path(output_dir, paste0(output_prefix, "_MAP.tif"))
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
  
  # --- 13. Display Plots in Console ---
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
  
  # --- 14. Return Results ---
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