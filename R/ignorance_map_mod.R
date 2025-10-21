#' @title Ignorance map (modernized minimal working version)
#' @description Minimal working replacement of ignorance_map using sf + terra.
#'              Produces Map of Relative Floristic Ignorance (MRFI) and species richness (RICH).
#' @param data_flor data.frame with columns: Taxon, Long, Lat, uncertainty (m), year
#' @param site sf or sp polygon (study area). If sp, will be coerced; expected geographic initially (EPSG:4326) or any CRS.
#' @param year_study numeric year (defaults to current year)
#' @param excl_areas optional sf or sp polygon to exclude from suitable area (will be transformed)
#' @param CRS.new numeric EPSG code for projected CRS (default 3035)
#' @param tau percentual taxa loss in 100 years (0 <= tau < 100)
#' @param cellsize raster cell size in meters (numeric)
#' @param verbose logical
#' @return list with MRFI (terra SpatRaster), RICH (terra SpatRaster), Uncertainties (data.frame), Statistics (data.frame)
#' @export
ignorance_map_mod <- function(data_flor, site, year_study = NULL, excl_areas = NULL,
                              CRS.new = 3035, tau = 20, cellsize = 2000, verbose = TRUE) {
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  ## ----- basic checks -----
  if (is.null(year_study)) year_study <- as.numeric(format(Sys.Date(), "%Y"))
  if (!("Taxon" %in% names(data_flor) && "Long" %in% names(data_flor) &&
        "Lat" %in% names(data_flor) && "uncertainty" %in% names(data_flor) &&
        "year" %in% names(data_flor))) {
    stop("data_flor must contain columns: Taxon, Long, Lat, uncertainty, year")
  }
  if (any(data_flor$year > year_study)) {
    warning("Some occurrence dates are more recent than year_study")
  }
  if (!(tau >= 0 && tau < 100)) stop("0 <= tau < 100 is required")
  if (!is.numeric(cellsize) || cellsize <= 0) stop("cellsize must be positive (meters)")
  
  ## ----- coerce site and excl_areas to sf -----
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  if (!inherits(site, "sf")) stop("site must be an sf or sp (Spatial*) polygon object")
  site <- sf::st_make_valid(site)
  
  if (!is.null(excl_areas)) {
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    excl_areas <- sf::st_make_valid(excl_areas)
  }
  
  ## ----- transform inputs to projected CRS (meters) -----
  target_crs <- sf::st_crs(as.integer(CRS.new))
  if (is.na(target_crs)) stop("Invalid CRS.new EPSG code")
  site_proj <- sf::st_transform(site, target_crs)
  if (!is.null(excl_areas)) excl_proj <- sf::st_transform(excl_areas, target_crs) else excl_proj <- NULL
  
  ## ----- create points sf from data_flor (assume input long/lat decimal degrees if not specified) -----
  pts <- sf::st_as_sf(data_flor, coords = c("Long", "Lat"), crs = 4326, remove = FALSE) # assume WGS84 input
  pts_proj <- sf::st_transform(pts, target_crs)
  pts_proj$record_id <- seq_len(nrow(pts_proj))
  
  ## ----- filter records with very small uncertainty relative to cell size (like original) -----
  if (any(2 * pts_proj$uncertainty < (cellsize / 20))) {
    bad <- which(2 * pts_proj$uncertainty < (cellsize / 20))
    stop("There are ", length(bad), " occurrence records with uncertainty too small relative to cellsize. Inspect rows: ",
         paste(head(bad, 10), collapse = ", "))
  }
  
  ## ----- build base raster using terra -----
  # extent = bounding box of records that intersect site (we will crop later)
  # compute buffers first to know extent
  msg("Creating buffers and computing extent...")
  buff_sfc <- sf::st_buffer(pts_proj, dist = pts_proj$uncertainty) # vectorized; returns sfc_GEOMETRY
  # optionally remove exclusion areas from buffers (as original did for cont==1)
  if (!is.null(excl_proj)) {
    # subtract excl polygons from buffers
    buff_sfc <- sf::st_difference(buff_sfc, sf::st_union(excl_proj))
  }
  
  # combined extent (add margin equal to max uncertainty)
  bbox <- sf::st_bbox(sf::st_union(buff_sfc))
  margin <- max(pts_proj$uncertainty, na.rm = TRUE)
  xmin <- bbox["xmin"] - margin; xmax <- bbox["xmax"] + margin
  ymin <- bbox["ymin"] - margin; ymax <- bbox["ymax"] + margin
  
  # create terra raster with desired cellsize (in same CRS units: meters)
  r_template <- terra::rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                            resolution = cellsize, crs = sf::st_crs(target_crs)$wkt)
  # ensure NA initial
  terra::values(r_template) <- NA
  
  ## ----- compute richness raster (counts of unique taxa per cell) -----
  msg("Computing species richness (RICH) ...")
  # find cell indices for each point
  pts_xy <- sf::st_coordinates(pts_proj)
  cell_idx <- terra::cellFromXY(r_template, pts_xy)
  # combine taxon and cell -> count unique taxa per cell
  df_cells <- data.frame(cell = cell_idx, Taxon = pts_proj$Taxon, stringsAsFactors = FALSE)
  df_cells <- df_cells[!is.na(df_cells$cell), ]
  if (nrow(df_cells) == 0) stop("No points fall within raster extent - check data and cellsize")
  # count unique taxa per cell
  ux <- tapply(df_cells$Taxon, df_cells$cell, function(x) length(unique(x)))
  r_rich <- r_template
  terra::values(r_rich)[] <- 0
  idx <- as.integer(names(ux))
  terra::values(r_rich)[idx] <- as.numeric(ux)
  
  ## ----- per-record scores (spatial, temporal, st_ignorance) -----
  msg("Computing per-record spatial & temporal scores ...")
  # compute number of cells inside each buffer (spatial extent in raster cells)
  # use terra::cells() via terra::rasterize with getCover-like behavior:
  # but simpler: for each buffer, compute area in m2 and divide by cell area approx -> approximate number of cells
  cell_area <- (cellsize * cellsize)
  # safer: compute exact count using rasterize with mask
  # We'll rasterize each buffer to the template and count non-NA cells
  nrec <- nrow(pts_proj)
  spatial_count <- integer(nrec)
  for (i in seq_len(nrec)) {
    single_buffer <- buff_sfc[[i]]  # extract the i-th geometry as a single sfc object
    # skip empty geometries
    if (is.null(single_buffer) || sf::st_is_empty(single_buffer)) {
      spatial_count[i] <- 0
      next
    }
    tmpv <- terra::rasterize(terra::vect(sf::st_as_sf(single_buffer)), r_template, field = 1, touches = TRUE)
    spatial_count[i] <- sum(!is.na(terra::values(tmpv)))
  }
  if (any(spatial_count == 0)) {
    # keep them but avoid divide-by-zero: set to 1 cell (very tiny area)
    spatial_count[spatial_count == 0] <- 1
  }
  spatial_score <- 1 / spatial_count
  time_score <- (1 - tau/100) ^ ((year_study - pts_proj$year) / 100)
  st_ignorance_vals <- spatial_score * time_score
  
  pts_proj$spatial_score <- spatial_score
  pts_proj$time_score <- time_score
  pts_proj$st_ignorance <- st_ignorance_vals
  
  ## ----- For each taxon, rasterize buffers with per-record st_ignorance and take max across records of the taxon -----
  msg("Rasterizing per-taxon ignorance layers (this may take time)...")
  taxa <- unique(pts_proj$Taxon)
  n_taxa <- length(taxa)
  # initialize an accumulation raster (sum across taxa)
  raster_sum <- r_template
  terra::values(raster_sum)[] <- 0
  
  pb <- utils::txtProgressBar(min = 0, max = n_taxa, style = 3)
  for (ti in seq_along(taxa)) {
    tname <- taxa[ti]
    setTxtProgressBar(pb, ti)
    idx <- which(pts_proj$Taxon == tname)
    if (length(idx) == 0) next
    # build sf of buffers only for these records and add st_ignorance as attribute
    bufs <- sf::st_as_sf(data.frame(st_ignorance = pts_proj$st_ignorance[idx]), geometry = buff_sfc[idx])
    # union not wanted because we want per-record values and then max; but we will rasterize each record and take cell-wise max
    # rasterize each record: create a raster (initialized NA) with each record's value; then cell-wise max across records
    # Instead of stacking many rasters, we can create a raster of zeros then for each record set cells with its value where larger than current.
    tax_r <- r_template
    terra::values(tax_r)[] <- NA
    # loop records of the taxon
    for (j in seq_len(nrow(bufs))) {
      rec_geom <- bufs[j, , drop = FALSE]
      val <- as.numeric(bufs$st_ignorance[j])
      # rasterize for this single record (cells covered -> val)
      rec_r <- terra::rasterize(terra::vect(rec_geom), r_template, field = val, touches = TRUE)
      # ensure rec_r has NA where not covered; then update tax_r with pmax of existing and rec_r
      # terra stores NA as NA; use cell-based compare
      cur_vals <- terra::values(tax_r)
      rec_vals <- terra::values(rec_r)
      # if tax_r is entirely NA, simply assign rec_vals
      if (all(is.na(cur_vals))) {
        terra::values(tax_r) <- rec_vals
      } else {
        # pmax treating NAs appropriately: where both NA => NA, where one is NA => take other
        newvals <- cur_vals
        to_update <- !is.na(rec_vals) & (is.na(cur_vals) | rec_vals > cur_vals)
        newvals[to_update] <- rec_vals[to_update]
        terra::values(tax_r) <- newvals
      }
    }
    # replace NA with 0 (no coverage by that taxon)
    tax_vals <- terra::values(tax_r)
    tax_vals[is.na(tax_vals)] <- 0
    terra::values(tax_r) <- tax_vals
    # accumulate sum across taxa
    raster_sum_vals <- terra::values(raster_sum)
    raster_sum_vals[is.na(raster_sum_vals)] <- 0
    raster_sum_vals <- raster_sum_vals + terra::values(tax_r)
    terra::values(raster_sum) <- raster_sum_vals
  }
  close(pb)
  
  ## ----- compute MRFI: r.max - raster_sum (as original) -----
  msg("Computing MRFI (rescaling)...")
  rmax <- max(terra::values(raster_sum), na.rm = TRUE)
  mrfi_r <- raster_sum
  terra::values(mrfi_r) <- rmax - terra::values(raster_sum)
  # mask MRFI to study site (and remove excluded areas)
  v_site <- terra::vect(site_proj)
  mrfi_r_masked <- terra::mask(terra::crop(mrfi_r, v_site), v_site)
  # also compute richness masked similarly
  rich_masked <- terra::mask(terra::crop(r_rich, v_site), v_site)
  
  ## ----- produce outputs & stats -----
  end_time <- Sys.time()
  stats_names <- c("Started", "Finished", "Elapsed time", "CRS (EPSG code)",
                   "Cell size (km)", "100 years % loss ratio (tau)",
                   "Total occurrence within", "Total occurrences computed",
                   "Occurrence uncertainty (median value, m)", "Occurrence dates (median value, year)")
  values <- c(as.character(start_time), as.character(end_time), as.character(round(end_time - start_time, 2)),
              as.character(CRS.new), cellsize / 1000, tau,
              sum(!is.na(terra::values(rich_masked))), nrow(pts_proj),
              round(median(pts_proj$uncertainty, na.rm = TRUE)), round(median(pts_proj$year, na.rm = TRUE)))
  statistics <- data.frame(Name = stats_names, Value = values, stringsAsFactors = FALSE)
  
  msg("Done.")
  
  out <- list(
    MRFI = mrfi_r_masked,
    RICH = rich_masked,
    Uncertainties = data.frame(uncertainty = pts_proj$uncertainty, year = pts_proj$year, Taxon = pts_proj$Taxon),
    Statistics = statistics
  )
  return(out)
}
