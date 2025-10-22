#' @title Ignorance map (Modernized)
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
#' @importFrom sf st_as_sf st_make_valid st_crs st_transform st_buffer st_intersects st_intersection st_union
#' @importFrom terra rast values ncell cellFromXY mask crop crs vect rasterize sprc extract
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom utils txtProgressBar setTxtProgressBar
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
#' 

ignorance_map_mod <- function(data_flor, site, year_study = NULL, excl_areas = NULL,
                              CRS.new = 3035, tau, cellsize, verbose = TRUE,
                              check_overlap = TRUE) {
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  # --- 1. Settings and Input Validation ---
  msg("Checking settings and inputs...")
  if (is.null(year_study)) year_study <- as.numeric(format(Sys.Date(), "%Y"))
  if (max(data_flor$year) > year_study) warning("Some occurrence dates are more recent than year_study")
  if (tau < 0 || tau >= 100) stop("0 <= tau < 100 is required.")
  
  req_cols <- c("Taxon", "Long", "Lat", "uncertainty", "year")
  if (!all(req_cols %in% names(data_flor))) {
    stop("data_flor must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  if (any(2 * data_flor$uncertainty < (cellsize / 20))) {
    stop("Some records have uncertainty too small vs. cellsize. They may be lost.")
  }
  
  msg("Inputs validated.")
  
  # --- 2. Coerce to SF and Reproject ---
  msg(paste("Reprojecting inputs to EPSG:", CRS.new))
  target_crs <- sf::st_crs(CRS.new)
  
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  if (is.na(sf::st_crs(site))) {
    msg("Input 'site' has no CRS. Assuming EPSG:4326.")
    sf::st_crs(site) <- 4326
  }
  site_proj <- sf::st_transform(sf::st_make_valid(site), target_crs)
  
  excl_proj <- NULL
  if (!is.null(excl_areas)) {
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    if (is.na(sf::st_crs(excl_areas))) {
      msg("Input 'excl_areas' has no CRS. Assuming EPSG:4326.")
      sf::st_crs(excl_areas) <- 4326
    }
    excl_proj <- sf::st_union(sf::st_transform(sf::st_make_valid(excl_areas), target_crs))
  }
  
  if (!inherits(data_flor, "sf")) {
    data_flor <- sf::st_as_sf(data_flor, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)
  }
  pts_proj <- sf::st_transform(data_flor, target_crs)
  
  # --- Optional: Check overlap ---
  if (check_overlap) {
    overlap_check <- sf::st_intersects(pts_proj, site_proj, sparse = FALSE)
    n_overlap <- sum(overlap_check)
    msg(paste("Number of occurrence points overlapping the site:", n_overlap))
    
    if (n_overlap == 0) {
      warning("No points overlap the study area. Check projections or coordinates.")
    } else {
      plot(sf::st_geometry(site_proj), col = NA, border = "black", main = "Points vs Site")
      points(sf::st_coordinates(pts_proj), col = "red", pch = 20, cex = 0.5)
    }
  }
  
  # --- 3. Filter Records ---
  msg("Creating buffers and filtering records...")
  all_buffers_sf <- sf::st_buffer(pts_proj, dist = pts_proj$uncertainty)
  intersects_idx <- sf::st_intersects(all_buffers_sf, site_proj, sparse = FALSE)
  pts_computed <- pts_proj[which(intersects_idx[, 1]), ]
  buffers_computed_sf <- all_buffers_sf[which(intersects_idx[, 1]), ]
  
  if (nrow(pts_computed) == 0) stop("No occurrence buffers intersect the study area.")
  if (!is.null(excl_proj)) buffers_computed_sf <- sf::st_difference(buffers_computed_sf, excl_proj)
  msg(paste("Retained", nrow(pts_computed), "records for computation."))
  
  # --- 4. Create Template Raster ---
  msg("Creating template raster...")
  
  combined_ext <- terra::ext(site_proj) + terra::ext(pts_computed)
  r_template <- terra::rast(
    extent = combined_ext,
    resolution = cellsize,
    crs = target_crs # The CRS derived earlier
  )
  terra::values(r_template) <- NA
  
  # --- 5. Calculate Richness Map ---
  msg("Calculating species richness map (RICH)...")
  r_rich <- terra::setValues(r_template, 0)
  
  pts_cells <- terra::cellFromXY(r_template, sf::st_coordinates(pts_computed))
  cell_taxa_df <- data.frame(cell = pts_cells, Taxon = pts_computed$Taxon)
  cell_taxa_df <- cell_taxa_df[!is.na(cell_taxa_df$cell), ]
  
  taxa_per_cell <- tapply(cell_taxa_df$Taxon, cell_taxa_df$cell, function(x) length(unique(x)))
  if (length(taxa_per_cell) > 0) r_rich[as.numeric(names(taxa_per_cell))] <- as.numeric(taxa_per_cell)
  
  # --- 6. Calculate Per-Record Scores ---
  msg("Calculating spatio-temporal scores...")
  r_cells <- r_template
  terra::values(r_cells) <- 1:terra::ncell(r_cells)
  
  v_buff <- terra::vect(buffers_computed_sf)
  cell_data <- terra::extract(r_cells, v_buff, cells = TRUE)
  counts_per_id <- table(cell_data$ID)
  spatial_count <- rep(1, nrow(pts_computed))
  if (length(counts_per_id) > 0) spatial_count[as.integer(names(counts_per_id))] <- counts_per_id
  pts_computed$spatial_score <- 1 / spatial_count
  pts_computed$time_score <- (1 - (tau / 100))^((year_study - pts_computed$year) / 100)
  pts_computed$st_ignorance <- pts_computed$spatial_score * pts_computed$time_score
  
  # --- 7. Rasterize Per-Taxon Ignorance with Progress Bar ---
  msg("Drafting Map of Relative Floristic Ignorance (safe with progress)...")
  taxa_list <- unique(pts_computed$Taxon)
  pb <- utils::txtProgressBar(min = 0, max = length(taxa_list), style = 3)
  
  raster_list <- lapply(seq_along(taxa_list), function(i) {
    tname <- taxa_list[i]
    pts_taxon <- pts_computed[pts_computed$Taxon == tname, ]
    bufs_taxon_sf <- buffers_computed_sf[pts_computed$Taxon == tname, ]
    bufs_taxon_sf$st_ignorance <- pts_taxon$st_ignorance
    
    tax_r <- terra::rasterize(
      terra::vect(bufs_taxon_sf),
      r_template,
      field = "st_ignorance",
      fun = "max",
      touches = TRUE
    )
    tax_r[is.na(tax_r)] <- 0
    
    utils::setTxtProgressBar(pb, i)
    return(tax_r)
  })
  close(pb)
  
  # --- 8. Finalize and Mask Rasters (UPDATED to Avoid NA Fill) ---
  msg("Finalizing rasters...")
  
  # 1. Keep only valid SpatRaster objects (as before)
  raster_list_valid <- raster_list[sapply(raster_list, function(x) inherits(x, "SpatRaster"))]
  
  if (length(raster_list_valid) > 0) {
    # Combine all rasters into one multilayer SpatRaster
    raster_stack <- terra::rast(raster_list_valid)
    
    # Sum across layers
    raster_sum <- terra::app(raster_stack, fun = sum, na.rm = TRUE)
  } else {
    raster_sum <- r_template
    terra::values(raster_sum) <- 0
  }
  
  # 2. Rescale to MRFI
  rmax <- max(terra::values(raster_sum), na.rm = TRUE)
  if (!is.finite(rmax)) rmax <- 0
  mrfi_r <- rmax - raster_sum # This is the un-masked MRFI map
  
  # 3. Create the site vector and a template for the site's bounding box
  site_vect <- terra::vect(site_proj)
  site_ext <- terra::ext(site_proj)
  site_box_template <- terra::rast(
    extent = site_ext,
    resolution = cellsize,
    crs = target_crs 
  )
  
  # 4. Mask to the actual study polygon (sets cells outside the polygon to NA)
  mrfi_masked <- terra::mask(mrfi_r, site_vect)
  rich_masked <- terra::mask(r_rich, site_vect)
  
  # 5. Fill cells *inside* the park boundary with zero for RICH if no data exists
  rich_masked <- terra::classify(rich_masked, cbind(NA, 0), others=TRUE)
  
  # 6. Final Step: Trim the Raster
  # The trim function cuts off rows and columns that are entirely NA.
  mrfi_final <- terra::mask(terra::crop(mrfi_r, site_vect), site_vect)
  rich_final <- terra::mask(terra::crop(r_rich, site_vect), site_vect)
  
  # --- 9. Compile Statistics ---
  msg("Compiling statistics...")
  end_time <- Sys.time()
  
  # 1. Ensure pts_in_site is calculated (needed for the original statistics)
  pts_in_site <- sf::st_intersection(pts_proj, site_proj)
  
  # 2. CALCULATE EXCLUDED RECORDS
  # total_initial_records must be available from the function start (e.g., nrow(data_flor))
  total_initial_records <- nrow(data_flor) 
  total_used_records <- nrow(pts_computed) # Final points used for MRFI calc
  excluded_records_count <- total_initial_records - total_used_records
  
  # 3. COMPILE FINAL DATA FRAME (retaining original structure)
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
  
  
  # --- 10. Plotting and Output Generation ---
  msg("Generating plots and output files...")
  
  # --- Prepare Data for ggplot ---
  # 1. MRFI Raster to DF
  test_spdf <- as.data.frame(mrfi_final, xy = TRUE)
  colnames(test_spdf) <- c("x", "y", "value")
  
  # 2. RICH Raster to DF
  test_spdf2 <- as.data.frame(rich_final, xy = TRUE)
  colnames(test_spdf2) <- c("x", "y", "value")
  
  # 3. Site Polygon for ggplot (sf to data.frame)
  # Uses the modern st_as_sf and ggplot2::fortify equivalent
  tip1 <- site_proj
  
  mrfi_max_val <- terra::global(mrfi_final, "max", na.rm = TRUE)$max
  rich_max_val <- terra::global(rich_final, "max", na.rm = TRUE)$max
  
  # Create a sequence of 5 breaks for cleaner legend display
  mrfi_breaks <- seq(0, mrfi_max_val, length.out = 12)
  rich_breaks <- seq(0, rich_max_val, length.out = 12)
  
  # --- Plot n° 1: MRFI ---
  p1 <- ggplot2::ggplot(test_spdf) +
    ggplot2::coord_equal() + ggplot2::theme_classic() +
    ggplot2::labs(fill = "IFI") +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical', legend.key.width = grid::unit(0.6, "cm")) +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") +
    ggplot2::scale_fill_distiller(
      palette = "Spectral", direction = -1, guide = ggplot2::guide_legend(),
      limits = c(0, mrfi_max_val),
      breaks = mrfi_breaks,
      labels = round(mrfi_breaks, 0) 
    ) +
    ggplot2::geom_tile(mapping = ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value), alpha = 0.8, colour = "black", linewidth = 0.1) +
    # Use sf::st_geometry to plot the outline cleanly
    ggplot2::geom_sf(data = tip1, fill = NA, color = "black", size = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("Map of Relative Floristic Ignorance (MRFI)")
  
  # --- Plot n° 2: Species Richness (FIXED) ---
  p2 <- ggplot2::ggplot(test_spdf2) + 
    ggplot2::coord_equal() + 
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical', legend.key.width = grid::unit(0.6, "cm")) +
    ggplot2::geom_tile(mapping = ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value), 
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::geom_sf(data = tip1, fill = NA, color = "black", size = 1, inherit.aes = FALSE) +
    ggplot2::scale_fill_distiller(
      palette = "Spectral", direction = +1, guide = ggplot2::guide_legend(title = "Value"),
      limits = c(0, rich_max_val),
      breaks = rich_breaks,
      labels = round(rich_breaks, 0) 
    ) +
    ggplot2::ggtitle("Species richness map (without uncertainties)") +
    ggplot2::xlab("Longitude") + 
    ggplot2::ylab("Latitude") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Value"))
  
  # --- Plot n° 3: Temporal Uncertainty (Year) ---
  p3 <- ggplot2::ggplot(pts_computed) +
    ggplot2::aes(x = .data$year, y = ggplot2::after_stat(density)) +
    ggplot2::geom_histogram(alpha = 0.6, fill = "#FF6666", binwidth = diff(range(pts_computed$year)) / 30) +
    ggplot2::coord_cartesian(xlim = c(min(pts_computed$year), year_study)) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    ggplot2::ggtitle("Occurrence date") +
    ggplot2::xlab("Year") + ggplot2::ylab("Frequency") +
    ggplot2::theme_classic()
  
  # --- Plot n° 4: Spatial Uncertainty (Uncertainty) ---
  p4 <- ggplot2::ggplot(pts_computed) +
    ggplot2::aes(x = .data$uncertainty, y = ggplot2::after_stat(density)) +
    ggplot2::geom_histogram(alpha = 0.6, fill = "#FF6666", binwidth = diff(range(pts_computed$uncertainty)) / 30) +
    ggplot2::coord_cartesian(xlim = c(min(pts_computed$uncertainty), max(pts_computed$uncertainty))) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    ggplot2::ggtitle("Occurrence spatial uncertainty") +
    ggplot2::xlab("Uncertainty (m)") + ggplot2::ylab("Frequency") +
    ggplot2::theme_classic()
  
  # --- PDF Output ---
  grDevices::pdf("Ignorance_output.pdf", onefile = TRUE)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  grid::grid.draw(gridExtra::grid.arrange(top = "Summary statistics", gridExtra::tableGrob(statistics_df)))
  grDevices::dev.off()
  
  # --- File Output ---
  terra::writeRaster(mrfi_final, filename = "MAPignorance.tif", overwrite = TRUE)
  utils::write.csv(taxa_list, row.names = FALSE, "Taxa considered to compute the Map of Relative Floristic Ignorance (MRFI).csv")
  msg(paste0("Done! The files have been saved here: ", getwd()))
  
  # --- Print Plots to Console ---
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  grid::grid.draw(gridExtra::grid.arrange(top = "Summary statistics", gridExtra::tableGrob(statistics_df)))
  
  # --- Return List (Section 11) ---
  
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
  
  msg("Done.")
}
