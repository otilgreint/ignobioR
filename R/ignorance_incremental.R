#' @title Calculate Incremental MRFI Update (Internal)
#'
#' @description
#' Internal function that calculates updated MRFI using incremental approach.
#' Maintains fixed baseline (T1's max_ignorance) while applying temporal decay
#' to existing observations and incorporating new observations.
#'
#' @param data_flor_initial Complete initial dataset (for aging baseline taxa)
#' @param data_flor_new ONLY new observations (records not in initial dataset)
#' @param mrfi_baseline Output from ignorance_map() for initial survey
#' @param year_updated Year for updated survey
#' @param tau Temporal decay parameter (extracted from baseline if NULL)
#' @param cellsize Cell size (extracted from baseline if NULL)
#' @param CRS.new CRS code (extracted from baseline if NULL)
#' @param excl_areas Exclusion areas
#' @param site Study area polygon
#' @param use_coverage_weighting Use coverage-weighted rasterization
#' @param verbose Print progress messages
#' @param site_buffer Expand site boundary
#' @param buffer_width Buffer width if site_buffer = TRUE
#' @param mask_method Masking method
#'
#' @return List similar to ignorance_map() output:
#' \itemize{
#'   \item{MRFI}{ Updated MRFI raster (anchored to baseline)}
#'   \item{RICH}{ Updated richness raster}
#'   \item{baseline_max_ignorance}{ Maximum ignorance from T1 (reference)}
#'   \item{Statistics}{ Summary statistics}
#' }
#'
#' @importFrom terra rast app global writeRaster values rasterize vect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @keywords internal
#' @noRd

calculate_incremental_mrfi_internal <- function(
    data_flor_initial,
    data_flor_new,
    mrfi_baseline,
    year_updated,
    tau = NULL,
    cellsize = NULL,
    CRS.new = NULL,
    excl_areas = NULL,
    site = NULL,
    use_coverage_weighting = TRUE,
    verbose = TRUE,
    site_buffer = FALSE,
    buffer_width = NULL,
    mask_method = "touches",
    ...
) {
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  msg("  [INCREMENTAL MRFI CALCULATION]")
  
  # --- SAFETY CHECK: Remove incomplete records internally as a fallback ---
  data_flor_initial <- data_flor_initial[!is.na(data_flor_initial$year) & !is.na(data_flor_initial$Taxon), ]
  data_flor_new <- data_flor_new[!is.na(data_flor_new$year) & !is.na(data_flor_new$Taxon), ]
  # ----------------------------------------------------------------------
  
  # Extract parameters from baseline if not provided
  if (is.null(tau)) {
    tau <- as.numeric(mrfi_baseline$Statistics$Value[
      mrfi_baseline$Statistics$Statistic == "100 years % loss ratio (tau)"
    ])
  }
  
  if (is.null(cellsize)) {
    cellsize <- as.numeric(mrfi_baseline$Statistics$Value[
      mrfi_baseline$Statistics$Statistic == "Cell size (m)"
    ])
  }
  
  if (is.null(CRS.new)) {
    CRS.new <- as.numeric(mrfi_baseline$Statistics$Value[
      mrfi_baseline$Statistics$Statistic == "CRS (EPSG code)"
    ])
  }
  
  # Extract baseline year
  baseline_year_str <- mrfi_baseline$Statistics$Value[
    mrfi_baseline$Statistics$Statistic == "Started"
  ]
  # Handle case where "Started" might be a full date string or just a year
  year_initial <- tryCatch({
    as.numeric(substr(baseline_year_str, 1, 4))
  }, warning = function(w) {
    min(data_flor_initial$year, na.rm = TRUE) # Fallback
  })
  
  msg(paste0("  Baseline year: ", year_initial))
  msg(paste0("  Updated year: ", year_updated))
  msg(paste0("  Time elapsed: ", year_updated - year_initial, " years"))
  
  # Get baseline maximum ignorance (this is our fixed reference!)
  max_ignorance_baseline <- max(terra::values(mrfi_baseline$MRFI), na.rm = TRUE)
  msg(paste0("  Baseline max ignorance: ", round(max_ignorance_baseline, 2)))
  
  # Create temp directory for this calculation
  temp_dir_base <- tempdir()
  temp_dir_aged <- file.path(temp_dir_base, "aged_taxa")
  temp_dir_new <- file.path(temp_dir_base, "new_taxa")
  temp_dir_combined <- file.path(temp_dir_base, "combined_taxa")
  
  dir.create(temp_dir_aged, showWarnings = FALSE, recursive = TRUE)
  dir.create(temp_dir_new, showWarnings = FALSE, recursive = TRUE)
  dir.create(temp_dir_combined, showWarnings = FALSE, recursive = TRUE)
  
  # ==========================================================================
  # STEP 1: Calculate per-taxon ignorance for BASELINE (T1) taxa
  # ==========================================================================
  
  msg("  STEP 1: Processing baseline taxa...")
  
  baseline_per_taxon <- ignorance_per_taxon(
    data_flor = data_flor_initial,
    site = site,
    year_study = year_initial,  # Calculate at original time
    excl_areas = excl_areas,
    CRS.new = CRS.new,
    tau = tau,
    cellsize = cellsize,
    template_raster = mrfi_baseline$MRFI,
    use_coverage_weighting = use_coverage_weighting,
    temp_dir = temp_dir_base,
    keep_in_memory = FALSE,  # Use temp files
    verbose = FALSE,
    site_buffer = site_buffer,
    buffer_width = buffer_width,
    mask_method = mask_method
  )
  
  msg(paste0("  Processed ", length(baseline_per_taxon$taxon_names), " baseline taxa"))
  
  # ==========================================================================
  # STEP 2: Age baseline taxa to T2 time
  # ==========================================================================
  
  msg("  STEP 2: Applying temporal decay to baseline taxa...")
  
  aged_files <- character(length(baseline_per_taxon$taxon_names))
  valid_indices <- rep(TRUE, length(baseline_per_taxon$taxon_names)) # Track valid files
  
  pb <- utils::txtProgressBar(min = 0, max = length(baseline_per_taxon$taxon_names), style = 3)
  
  for (i in seq_along(baseline_per_taxon$taxon_names)) {
    taxon_name <- baseline_per_taxon$taxon_names[i]
    
    # Get observations for this taxon
    taxon_obs <- baseline_per_taxon$taxon_observations[
      baseline_per_taxon$taxon_observations$Taxon == taxon_name,
    ]
    
    # --- ROBUSTNESS FIX: Skip if observation data is empty or invalid ---
    if (nrow(taxon_obs) == 0) {
      valid_indices[i] <- FALSE
      utils::setTxtProgressBar(pb, i)
      next
    }
    
    # Age this taxon's raster
    aged_file <- file.path(temp_dir_aged, 
                           paste0("aged_", i, "_", 
                                  gsub("[^A-Za-z0-9]", "_", taxon_name), ".tif"))
    
    tryCatch({
      age_taxon_raster(
        taxon_file_or_raster = baseline_per_taxon$taxon_files[i],
        taxon_obs = taxon_obs,
        year_from = year_initial,
        year_to = year_updated,
        tau = tau,
        cellsize = cellsize,
        template_raster = mrfi_baseline$MRFI,
        use_coverage_weighting = use_coverage_weighting,
        save_to_file = TRUE,
        output_path = aged_file
      )
      aged_files[i] <- aged_file
    }, error = function(e) {
      # Log error but continue loop
      warning(paste("Failed to age taxon:", taxon_name, "-", e$message))
      valid_indices[i] <- FALSE
    })
    
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Filter out failed taxons
  aged_files <- aged_files[valid_indices]
  aged_names <- baseline_per_taxon$taxon_names[valid_indices]
  
  msg(paste0("  Aged ", length(aged_files), " baseline taxa successfully"))
  
  # ==========================================================================
  # STEP 3: Calculate per-taxon ignorance for NEW observations
  # ==========================================================================
  
  if (nrow(data_flor_new) > 0) {
    msg(paste0("  STEP 3: Processing ", nrow(data_flor_new), " new observations..."))
    
    new_per_taxon <- ignorance_per_taxon(
      data_flor = data_flor_new,
      site = site,
      year_study = year_updated,  # Calculate at updated time
      excl_areas = excl_areas,
      CRS.new = CRS.new,
      tau = tau,
      cellsize = cellsize,
      template_raster = mrfi_baseline$MRFI,
      use_coverage_weighting = use_coverage_weighting,
      temp_dir = temp_dir_new,
      keep_in_memory = FALSE,
      verbose = FALSE,
      site_buffer = site_buffer,
      buffer_width = buffer_width,
      mask_method = mask_method
    )
    
    msg(paste0("  Processed ", length(new_per_taxon$taxon_names), " taxa from new observations"))
    
    new_taxa <- setdiff(new_per_taxon$taxon_names, baseline_per_taxon$taxon_names)
    existing_taxa_with_new <- intersect(new_per_taxon$taxon_names, baseline_per_taxon$taxon_names)
    
    msg(paste0("  - ", length(new_taxa), " completely new taxa"))
    msg(paste0("  - ", length(existing_taxa_with_new), " existing taxa with new observations"))
    
  } else {
    msg("  STEP 3: No new observations to process")
    new_per_taxon <- list(taxon_files = character(0), taxon_names = character(0))
  }
  
  # ==========================================================================
  # STEP 4: Combine aged baseline + new observations
  # ==========================================================================
  
  msg("  STEP 4: Combining aged baseline with new observations...")
  
  combined_result <- combine_taxon_rasters(
    aged_files = aged_files,
    aged_names = aged_names, # Use the filtered names
    new_files = new_per_taxon$taxon_files,
    new_names = new_per_taxon$taxon_names,
    template_raster = mrfi_baseline$MRFI,
    output_dir = temp_dir_combined,
    verbose = verbose
  )
  
  # ==========================================================================
  # STEP 5: Calculate final IRFI with fixed baseline
  # ==========================================================================
  
  msg("  STEP 5: Calculating updated IRFI with fixed baseline...")
  msg(paste0("  Loading ", length(combined_result$combined_files), " combined rasters..."))
  
  # Load all combined rasters with progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(combined_result$combined_files), style = 3)
  combined_rasters <- lapply(seq_along(combined_result$combined_files), function(i) {
    r <- terra::rast(combined_result$combined_files[i])
    utils::setTxtProgressBar(pb, i)
    return(r)
  })
  close(pb)
  
  msg("  Stacking and summing ignorance values...")
  raster_stack <- terra::rast(combined_rasters)
  ignorance_sum_updated <- terra::app(raster_stack, fun = sum, na.rm = TRUE)
  
  # Apply IRFI formula with T1's max_ignorance (KEY: fixed baseline!)
  IRFI_updated <- max_ignorance_baseline - ignorance_sum_updated
  
  msg("  Updated IRFI calculated (can contain negative values)")
  
  # ==========================================================================
  # STEP 6: Calculate updated richness
  # ==========================================================================
  
  msg("  Calculating updated species richness...")
  
  # Combine initial + new data for richness calculation
  data_flor_complete <- rbind(data_flor_initial, data_flor_new)
  
  # Simple point-based richness (no uncertainty buffers)
  pts_sf <- sf::st_as_sf(data_flor_complete, coords = c("Long", "Lat"), 
                         crs = 4326, remove = FALSE)
  crs_sf <- sf::st_crs(CRS.new)
  pts_proj <- sf::st_transform(pts_sf, crs_sf)
  
  pts_vect <- terra::vect(pts_proj)
  r_rich_updated <- terra::rasterize(
    pts_vect,
    mrfi_baseline$MRFI,  # Use same template
    field = "Taxon",
    fun = function(x) length(unique(x))
  )
  r_rich_updated[is.na(r_rich_updated)] <- 0
  
  # ==========================================================================
  # STEP 7: Cleanup temporary files
  # ==========================================================================
  
  msg("  Cleaning up temporary files...")
  
  cleanup_temp_rasters(baseline_per_taxon$taxon_files, verbose = FALSE)
  cleanup_temp_rasters(aged_files, verbose = FALSE)
  cleanup_temp_rasters(new_per_taxon$taxon_files, verbose = FALSE)
  cleanup_temp_rasters(combined_result$combined_files, verbose = FALSE)
  
  # Remove temp directories
  unlink(temp_dir_aged, recursive = TRUE)
  unlink(temp_dir_new, recursive = TRUE)
  unlink(temp_dir_combined, recursive = TRUE)
  
  # ==========================================================================
  # STEP 8: Compile statistics
  # ==========================================================================
  
  end_time <- Sys.time()
  
  statistics_df <- data.frame(
    Statistic = c(
      "Calculation method",
      "Started", "Finished", "Elapsed time (secs)",
      "Baseline year", "Updated year", "Time elapsed (years)",
      "Baseline max ignorance (fixed reference)",
      "Total taxa in baseline", "Total new taxa discovered",
      "Total taxa in updated dataset",
      "New observations processed"
    ),
    Value = c(
      "Incremental (fixed baseline)",
      as.character(start_time), as.character(end_time),
      round(as.numeric(difftime(end_time, start_time, units = "secs")), 2),
      year_initial, year_updated, year_updated - year_initial,
      round(max_ignorance_baseline, 2),
      length(baseline_per_taxon$taxon_names),
      length(setdiff(new_per_taxon$taxon_names, baseline_per_taxon$taxon_names)),
      length(combined_result$combined_names),
      nrow(data_flor_new)
    )
  )
  
  msg("  Incremental MRFI calculation complete!")
  
  # ==========================================================================
  # STEP 9: Return results
  # ==========================================================================
  
  return(list(
    MRFI = IRFI_updated,
    RICH = r_rich_updated,
    baseline_max_ignorance = max_ignorance_baseline,
    Statistics = statistics_df
  ))
}