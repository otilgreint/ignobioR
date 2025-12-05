#' @title Age Taxon Rasters for Temporal Decay (Internal)
#'
#' @description
#' Internal function that applies temporal decay to taxon-specific ignorance
#' rasters by recalculating temporal scores from an old year to a new year.
#'
#' @param taxon_file_or_raster Either a file path to a raster or a SpatRaster object
#' @param taxon_obs Subset of observation data for this specific taxon (from pts_computed)
#' @param year_from Original calculation year
#' @param year_to New calculation year
#' @param tau Temporal decay parameter
#' @param cellsize Raster cell size in meters
#' @param template_raster SpatRaster template for consistency
#' @param use_coverage_weighting Logical. Use coverage-weighted rasterization
#' @param save_to_file Logical. If TRUE, save result to temp file
#' @param output_path Character. Output path if save_to_file = TRUE
#'
#' @return Either a SpatRaster (if save_to_file = FALSE) or file path (if TRUE)
#'
#' @importFrom terra rast vect rasterize
#' @importFrom sf st_buffer st_sf st_geometry
#'
#' @keywords internal
#' @noRd

age_taxon_raster <- function(
    taxon_file_or_raster,  # This parameter becomes unused
    taxon_obs,
    year_from,
    year_to,
    tau,
    cellsize,
    template_raster,
    use_coverage_weighting = TRUE,
    save_to_file = FALSE,
    output_path = NULL
) {
  
  # Safety: Empty observations
  if (nrow(taxon_obs) == 0) {
    aged_raster <- template_raster
    terra::values(aged_raster) <- 0
    if (save_to_file) {
      terra::writeRaster(aged_raster, output_path, overwrite = TRUE)
      return(output_path)
    }
    return(aged_raster)
  }
  
  # Safety: Remove NA uncertainty
  if (any(is.na(taxon_obs$uncertainty))) {
    n_removed <- sum(is.na(taxon_obs$uncertainty))
    warning(paste0("Removing ", n_removed, " observations with NA uncertainty"))
    taxon_obs <- taxon_obs[!is.na(taxon_obs$uncertainty), ]
  }
  
  # === KEY FIX: Recalculate st_ignorance with aged temporal scores ===
  taxon_obs$temporal_score_aged <- (1 - tau/100)^((year_to - taxon_obs$year) / 100)
  taxon_obs$st_ignorance <- taxon_obs$spatial_score * taxon_obs$temporal_score_aged
  
  # Recreate buffers
  buffers_sf <- sf::st_buffer(taxon_obs, dist = taxon_obs$uncertainty)
  
  # Rasterize with UPDATED st_ignorance
  if (use_coverage_weighting) {
    aged_raster <- template_raster
    terra::values(aged_raster) <- 0
    
    for (j in seq_len(nrow(taxon_obs))) {
      single_buf <- sf::st_sf(geometry = sf::st_geometry(buffers_sf[j, ]))
      coverage_r <- terra::rasterize(terra::vect(single_buf), 
                                     template_raster, cover = TRUE)
      coverage_r[is.na(coverage_r)] <- 0
      
      # NOW using updated st_ignorance!
      weighted_r <- coverage_r * taxon_obs$st_ignorance[j]
      aged_raster <- max(aged_raster, weighted_r, na.rm = TRUE)
    }
  } else {
    # Binary approach
    bufs_clean <- sf::st_sf(
      st_ignorance = taxon_obs$st_ignorance,  # Updated values
      geometry = sf::st_geometry(buffers_sf)
    )
    
    aged_raster <- terra::rasterize(
      terra::vect(bufs_clean),
      template_raster,
      field = "st_ignorance",
      fun = "max",
      touches = TRUE
    )
    aged_raster[is.na(aged_raster)] <- 0
  }
  
  # Save or return
  if (save_to_file) {
    terra::writeRaster(aged_raster, output_path, overwrite = TRUE)
    return(output_path)
  } else {
    return(aged_raster)
  }
}


#' @title Combine Aged and New Taxon Rasters (Internal)
#'
#' @description
#' Internal function that combines aged baseline rasters with new observation
#' rasters, taking the maximum ignorance value per taxon per cell.
#'
#' @param aged_files Character vector of paths to aged taxon rasters
#' @param aged_names Character vector of taxon names for aged rasters
#' @param new_files Character vector of paths to new observation rasters
#' @param new_names Character vector of taxon names for new rasters
#' @param template_raster SpatRaster template
#' @param output_dir Directory for combined output files
#' @param verbose Logical. Print progress messages
#'
#' @return List with:
#' \itemize{
#'   \item{combined_files}{ Character vector of paths to combined rasters}
#'   \item{combined_names}{ Character vector of taxon names}
#' }
#'
#' @importFrom terra rast writeRaster
#'
#' @keywords internal
#' @noRd

combine_taxon_rasters <- function(
    aged_files,
    aged_names,
    new_files,
    new_names,
    template_raster,
    output_dir,
    verbose = TRUE
) {
  
  msg <- function(...) if (verbose) message(...)
  
  # Input validation
  if(length(aged_files) != length(aged_names)) {
    stop("Mismatch between aged files and taxon names.")
  }
  
  # Identify taxa categories
  existing_taxa_with_new <- intersect(aged_names, new_names)
  only_aged <- setdiff(aged_names, new_names)
  only_new <- setdiff(new_names, aged_names)
  
  msg(paste0("  Combining: ", length(existing_taxa_with_new), 
             " taxa with both aged + new observations"))
  msg(paste0("  Keeping: ", length(only_aged), " aged-only taxa"))
  msg(paste0("  Adding: ", length(only_new), " new-only taxa"))
  
  combined_files <- character()
  combined_names <- character()
  
  # Calculate total operations for progress bar
  total_ops <- length(existing_taxa_with_new) + length(only_aged) + length(only_new)
  pb <- utils::txtProgressBar(min = 0, max = total_ops, style = 3)
  current_op <- 0
  
  # 1. Process taxa with BOTH aged and new observations
  if (length(existing_taxa_with_new) > 0) {
    for (taxon in existing_taxa_with_new) {
      aged_file <- aged_files[aged_names == taxon]
      new_file <- new_files[new_names == taxon]
      
      # Load, Max, Write
      tryCatch({
        aged_r <- terra::rast(aged_file)
        new_r <- terra::rast(new_file)
        combined_r <- max(aged_r, new_r, na.rm = TRUE)
        
        output_file <- file.path(output_dir, 
                                 paste0("combined_", gsub("[^A-Za-z0-9]", "_", taxon), ".tif"))
        terra::writeRaster(combined_r, output_file, overwrite = TRUE)
        
        combined_files <- c(combined_files, output_file)
        combined_names <- c(combined_names, taxon)
      }, error = function(e) {
        warning(paste("Error combining raster for taxon:", taxon, "-", e$message))
      })
      
      current_op <- current_op + 1
      utils::setTxtProgressBar(pb, current_op)
    }
  }
  
  # 2. Copy aged-only taxa
  if (length(only_aged) > 0) {
    for (taxon in only_aged) {
      aged_file <- aged_files[aged_names == taxon]
      output_file <- file.path(output_dir, 
                               paste0("combined_", gsub("[^A-Za-z0-9]", "_", taxon), ".tif"))
      
      if(file.exists(aged_file)) {
        file.copy(aged_file, output_file, overwrite = TRUE)
        combined_files <- c(combined_files, output_file)
        combined_names <- c(combined_names, taxon)
      }
      
      current_op <- current_op + 1
      utils::setTxtProgressBar(pb, current_op)
    }
  }
  
  # 3. Copy new-only taxa
  if (length(only_new) > 0) {
    for (taxon in only_new) {
      new_file <- new_files[new_names == taxon]
      output_file <- file.path(output_dir, 
                               paste0("combined_", gsub("[^A-Za-z0-9]", "_", taxon), ".tif"))
      
      if(file.exists(new_file)) {
        file.copy(new_file, output_file, overwrite = TRUE)
        combined_files <- c(combined_files, output_file)
        combined_names <- c(combined_names, taxon)
      }
      
      current_op <- current_op + 1
      utils::setTxtProgressBar(pb, current_op)
    }
  }
  
  close(pb)
  
  return(list(
    combined_files = combined_files,
    combined_names = combined_names
  ))
}


#' @title Cleanup Temporary Taxon Rasters (Internal)
#'
#' @description
#' Internal utility to delete temporary raster files after incremental
#' MRFI calculation is complete.
#'
#' @param file_paths Character vector of file paths to delete
#' @param verbose Logical. Print messages
#'
#' @return Logical indicating success
#'
#' @keywords internal
#' @noRd

cleanup_temp_rasters <- function(file_paths, verbose = TRUE) {
  
  msg <- function(...) if (verbose) message(...)
  
  if (is.null(file_paths) || length(file_paths) == 0) {
    return(TRUE)
  }
  
  n_files <- length(file_paths)
  n_deleted <- 0
  
  for (f in file_paths) {
    if (file.exists(f)) {
      tryCatch({
        file.remove(f)
        n_deleted <- n_deleted + 1
      }, error = function(e) {
        warning(paste("Could not delete", f, ":", e$message))
      })
    }
  }
  
  if (verbose && n_deleted > 0) {
    msg(paste0("  Cleaned up ", n_deleted, " temporary raster files"))
  }
  
  return(n_deleted == n_files)
}


#' @title Identify New Records Between Two Datasets (Internal)
#'
#' @description
#' Internal function that identifies records present in updated dataset but not
#' in initial dataset. Matches records by Taxon, coordinates, and year.
#'
#' @param data_initial Initial floristic dataset
#' @param data_updated Updated floristic dataset
#'
#' @return Data frame containing only new records (rows in updated but not in initial)
#'
#' @keywords internal
#' @noRd

identify_new_records <- function(data_initial, data_updated) {
  
  # Convert sf to data.frame if needed (keep only relevant columns)
  cols_to_match <- c("Taxon", "Long", "Lat", "uncertainty", "year")
  
  if (inherits(data_initial, "sf")) {
    initial_df <- as.data.frame(data_initial)[, cols_to_match]
  } else {
    initial_df <- data_initial[, cols_to_match]
  }
  
  if (inherits(data_updated, "sf")) {
    updated_df <- as.data.frame(data_updated)[, cols_to_match]
  } else {
    updated_df <- data_updated[, cols_to_match]
  }
  
  # Create unique identifiers for matching
  # Round coordinates to avoid floating point issues
  initial_df$match_id <- paste(
    initial_df$Taxon,
    round(initial_df$Long, 6),
    round(initial_df$Lat, 6),
    initial_df$uncertainty,
    initial_df$year,
    sep = "_"
  )
  
  updated_df$match_id <- paste(
    updated_df$Taxon,
    round(updated_df$Long, 6),
    round(updated_df$Lat, 6),
    updated_df$uncertainty,
    updated_df$year,
    sep = "_"
  )
  
  # Identify new records (in updated but not in initial)
  new_indices <- !(updated_df$match_id %in% initial_df$match_id)
  
  # Return new records with all original columns from data_updated
  if (inherits(data_updated, "sf")) {
    data_flor_new <- data_updated[new_indices, ]
  } else {
    data_flor_new <- data_updated[new_indices, ]
  }
  
  return(data_flor_new)
}