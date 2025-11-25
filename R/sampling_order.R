#' @title Sampling Order Optimization
#'
#' @description
#' Generates an optimized sampling order for field plots based on environmental 
#' diversity (NDVI, DEM) to maximize species accumulation efficiency. Can optionally
#' work within pre-defined Operational Geographic Units (OGUs) for field logistics.
#'
#' @param best_solution_sf An sf object from sampleboost() containing plot locations 
#'   and environmental data (must have columns: ndvi_mean, and optionally elevation)
#' @param ogus Optional sf polygon object defining Operational Geographic Units (OGUs).
#'   If provided, plots are assigned to OGUs via spatial intersection, and ordering
#'   is done separately within each OGU. OGUs are typically pre-defined study area
#'   zones (e.g., "North Valley", "South Ridge") based on ecological, logistical,
#'   or administrative criteria.
#' @param ogu_id_col Character. Name of column in `ogus` containing OGU names/identifiers
#'   (default = "OGU"). Ignored if `ogus` is NULL.
#' @param ndvi_weight Numeric. Weight for NDVI in environmental distance calculation (default = 1)
#' @param dem_weight Numeric. Weight for elevation in environmental distance calculation (default = 0.5)
#' @param start_plot Can be:
#'   \itemize{
#'     \item "auto" (default): Automatically selects most environmentally extreme plot in each OGU
#'     \item "centroid": Selects plot closest to environmental centroid (most "average" plot) - useful if starting from accessible center
#'     \item Integer: Single plot ID to start from (if no OGUs)
#'     \item Named vector: Starting plot ID for each OGU, e.g. c("North" = 5, "South" = 12)
#'   }
#' @param verbose Logical. Print progress messages (default = TRUE)
#' @param output_dir Character. Directory for output files (default = working directory)
#' @param output_prefix Character. Prefix for output filenames (default = "SamplingOrder")
#'
#' @return A list containing:
#' \itemize{
#'   \item{ordered_plots}{ sf object with sampling_order and OGU columns added}
#'   \item{field_ready}{ Data frame with complete field information}
#'   \item{env_distance_matrix}{ Environmental distance matrix between all plots}
#'   \item{summary}{ Data frame with ordering statistics}
#' }
#' 
#' Additionally, two files are automatically saved to output_dir:
#' \itemize{
#'   \item{[prefix]_field-ready.csv}{ Complete field sheet with Plot_ID, Sampling_Order, 
#'   OGU (if applicable), coordinates (lat/lon), NDVI, elevation - everything in one file}
#'   \item{[prefix]_summary.txt}{ Summary statistics and field instructions}
#' }
#'
#' @section OGU Format:
#' The `ogus` parameter should be an sf object with POLYGON or MULTIPOLYGON geometry containing:
#' \itemize{
#'   \item **Geometry**: Polygon boundaries defining each OGU zone
#'   \item **ID column**: Character or factor column with OGU names (e.g., "North_Valley", "Lake_Shore")
#'   \item **CRS**: Must have a valid CRS (will be matched to plot CRS automatically)
#' }
#' 
#' **Example OGU creation:**
#' \preformatted{
#' # From shapefile
#' ogus <- st_read("study_zones.shp")
#' 
#' # Or create manually
#' library(sf)
#' north <- st_polygon(list(matrix(c(x_coords, y_coords), ncol=2)))
#' south <- st_polygon(list(matrix(c(x_coords2, y_coords2), ncol=2)))
#' ogus <- st_sf(
#'   OGU = c("North", "South"),
#'   geometry = st_sfc(north, south),
#'   crs = 3035
#' )
#' }
#'
#' @details
#' **Ordering Algorithm (Maximin Diversity)**:
#' Iteratively selects plots that maximize environmental diversity by choosing the 
#' plot that is most different from all previously selected plots. This ensures 
#' rapid saturation of the species accumulation curve.
#' 
#' **Environmental Distance**: Euclidean distance in standardized environmental space:
#' \deqn{d_{ij} = \sqrt{w_{ndvi} \times (NDVI_i - NDVI_j)^2 + w_{dem} \times (DEM_i - DEM_j)^2}}
#' 
#' **With OGUs**:
#' If OGUs are provided, maximin diversity ordering is applied separately within each
#' OGU. The user decides which OGU to sample first based on field logistics (accessibility,
#' weather, permits, etc.). Each OGU gets its own sampling sequence starting from 1.
#' 
#' **Without OGUs**:
#' Pure diversity-based ordering across all plots globally.
#'
#' @importFrom sf st_transform st_crs st_coordinates st_as_sf st_drop_geometry st_geometry st_sf st_intersects st_join
#' @importFrom stats dist
#' @importFrom utils write.csv
#' @export
#'
#' @examples
#' \dontrun{
#' # Run sampleboost first
#' result <- sampleboost(ndvi, ignorance, site, nplot = 50, plot_radius = 20, perm = 100)
#' 
#' # Simple diversity ordering (no OGUs) - auto selects most extreme plot
#' order1 <- sampling_order(result$best_solution_sf)
#' 
#' # Start from most "average" plot (closest to centroid) - useful for accessible center
#' order2 <- sampling_order(result$best_solution_sf, start_plot = "centroid")
#' 
#' # With user-defined OGUs - auto selects most extreme in each OGU
#' my_zones <- st_read("study_zones.shp")
#' order3 <- sampling_order(
#'   result$best_solution_sf, 
#'   ogus = my_zones,
#'   ogu_id_col = "zone_name"
#' )
#' 
#' # Specify starting plot for each OGU (based on access points!)
#' # E.g., Plot 5 is near north parking, Plot 23 is near south trail
#' order4 <- sampling_order(
#'   result$best_solution_sf, 
#'   ogus = my_zones,
#'   ogu_id_col = "zone_name",
#'   start_plot = c("North_Valley" = 5, "South_Ridge" = 23, "Lake_Shore" = 12)
#' )
#' 
#' # Custom output location
#' final <- sampling_order(
#'   result$best_solution_sf, 
#'   ogus = my_zones,
#'   ogu_id_col = "zone_name",
#'   start_plot = c("North_Valley" = 5, "South_Ridge" = 23),
#'   output_dir = "~/fieldwork",
#'   output_prefix = "LakeVico_2025"
#' )
#' }

sampling_order <- function(best_solution_sf, 
                           ogus = NULL,
                           ogu_id_col = "OGU",
                           ndvi_weight = 1, 
                           dem_weight = 0.5,
                           start_plot = "auto",
                           verbose = TRUE,
                           output_dir = file.path(getwd(), "output"),
                           output_prefix = "SamplingOrder") {
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  msg <- function(...) if (verbose) message(...)
  
  # Maximin ordering: iteratively select plot most different from already selected
  maximin_order <- function(dist_matrix, start_idx = 1) {
    n <- nrow(dist_matrix)
    selected <- integer(n)
    selected[1] <- start_idx
    
    for (i in 2:n) {
      remaining <- setdiff(1:n, selected[1:(i-1)])
      
      # For each remaining plot, find minimum distance to any selected plot
      min_dists <- sapply(remaining, function(r) {
        min(dist_matrix[r, selected[1:(i-1)]])
      })
      
      # Select plot with maximum of these minimum distances
      next_plot <- remaining[which.max(min_dists)]
      selected[i] <- next_plot
    }
    
    return(selected)
  }
  
  # Standardize to [0,1]
  standardize <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  # ============================================================================
  # SECTION 1: INPUT VALIDATION AND DATA EXTRACTION
  # ============================================================================
  
  msg("Validating inputs and extracting environmental data...")
  
  if (!inherits(best_solution_sf, "sf")) {
    stop("best_solution_sf must be an sf object (output from sampleboost).\n",
         "If loading from CSV, convert it first:\n",
         "  library(sf)\n",
         "  df <- read.csv('SampleBoost_best-solution.csv')\n",
         "  best_solution_sf <- st_as_sf(df, coords = c('x', 'y'), crs = 3035)")
  }
  
  # Check that geometry column exists
  geom_col <- attr(best_solution_sf, "sf_column")
  if (is.null(geom_col) || !geom_col %in% names(best_solution_sf)) {
    stop("best_solution_sf does not contain valid geometry.\n",
         "If loading from CSV, make sure to convert using st_as_sf() with coords parameter.")
  }
  
  # Check for required columns
  if (!"ndvi_mean" %in% names(best_solution_sf)) {
    stop("best_solution_sf must contain 'ndvi_mean' column")
  }
  
  has_dem <- "elevation" %in% names(best_solution_sf) && !all(is.na(best_solution_sf$elevation))
  
  n_plots <- nrow(best_solution_sf)
  
  # Create Plot_ID if it doesn't exist
  if (!"Plot_ID" %in% names(best_solution_sf)) {
    best_solution_sf$Plot_ID <- 1:n_plots
    msg("  Created Plot_ID column (1 to ", n_plots, ")")
  }
  
  # Validate OGUs if provided
  if (!is.null(ogus)) {
    if (!inherits(ogus, "sf")) {
      stop("ogus must be an sf object with polygon geometry")
    }
    
    if (!ogu_id_col %in% names(ogus)) {
      stop("Column '", ogu_id_col, "' not found in ogus object.\n",
           "Available columns: ", paste(names(ogus)[names(ogus) != attr(ogus, "sf_column")], collapse = ", "))
    }
    
    # Ensure CRS match
    if (!identical(sf::st_crs(ogus), sf::st_crs(best_solution_sf))) {
      msg("  Transforming OGUs to match plot CRS...")
      ogus <- sf::st_transform(ogus, sf::st_crs(best_solution_sf))
    }
    
    msg("  OGUs provided: ", paste(unique(ogus[[ogu_id_col]]), collapse = ", "))
  }
  
  # Extract environmental data
  ndvi_data <- best_solution_sf$ndvi_mean
  dem_data <- if (has_dem) best_solution_sf$elevation else rep(0, n_plots)
  coords <- sf::st_coordinates(best_solution_sf)
  
  msg(paste0("  ", n_plots, " plots detected"))
  msg(paste0("  DEM data available: ", ifelse(has_dem, "Yes", "No")))
  
  # ============================================================================
  # SECTION 2: CALCULATE ENVIRONMENTAL DISTANCE MATRIX
  # ============================================================================
  
  msg("Calculating environmental distance matrix...")
  
  # Standardize to [0,1]
  ndvi_std <- standardize(ndvi_data)
  dem_std <- if (has_dem) standardize(dem_data) else rep(0, n_plots)
  
  # Calculate environmental distance matrix
  env_matrix <- matrix(0, n_plots, n_plots)
  for (i in 1:n_plots) {
    for (j in 1:n_plots) {
      if (i != j) {
        env_matrix[i, j] <- sqrt(
          ndvi_weight * (ndvi_std[i] - ndvi_std[j])^2 +
            dem_weight * (dem_std[i] - dem_std[j])^2
        )
      }
    }
  }
  
  msg(paste0("  Mean environmental distance: ", round(mean(env_matrix[env_matrix > 0]), 3)))
  
  # ============================================================================
  # SECTION 3: ASSIGN PLOTS TO OGUs (IF PROVIDED)
  # ============================================================================
  
  ogu_assignments <- rep("All_Plots", n_plots)  # Default: single group
  
  if (!is.null(ogus)) {
    msg("Assigning plots to OGUs via spatial intersection...")
    
    # Spatial join to assign plots to OGUs
    plots_with_ogus <- sf::st_join(best_solution_sf, ogus[ogu_id_col], 
                                   join = sf::st_intersects, left = TRUE)
    
    ogu_assignments <- plots_with_ogus[[ogu_id_col]]
    
    # Handle plots not in any OGU
    if (any(is.na(ogu_assignments))) {
      n_unassigned <- sum(is.na(ogu_assignments))
      msg(paste0("  WARNING: ", n_unassigned, " plots not within any OGU (will be assigned to 'Unassigned')"))
      ogu_assignments[is.na(ogu_assignments)] <- "Unassigned"
    }
    
    # Convert to character for consistency
    ogu_assignments <- as.character(ogu_assignments)
    
    # Report OGU sizes
    ogu_sizes <- table(ogu_assignments)
    msg("  Plots per OGU:")
    for (ogu_name in names(ogu_sizes)) {
      msg(paste0("    ", ogu_name, ": ", ogu_sizes[ogu_name], " plots"))
    }
  }
  
  # ============================================================================
  # SECTION 4: APPLY MAXIMIN ORDERING (WITHIN OGUs IF PROVIDED)
  # ============================================================================
  
  msg("Applying maximin diversity ordering...")
  
  # Parse start_plot parameter
  start_plot_is_named_vector <- is.numeric(start_plot) && !is.null(names(start_plot))
  
  # Initialize
  sampling_order_within_ogu <- integer(n_plots)
  
  # Get unique OGUs
  unique_ogus <- unique(ogu_assignments)
  
  for (ogu_name in unique_ogus) {
    ogu_plots <- which(ogu_assignments == ogu_name)
    n_ogu_plots <- length(ogu_plots)
    
    if (!is.null(ogus) && length(unique_ogus) > 1) {
      msg(paste0("  Ordering OGU '", ogu_name, "' (", n_ogu_plots, " plots)..."))
    }
    
    # Extract submatrix for this OGU
    ogu_dist_matrix <- env_matrix[ogu_plots, ogu_plots]
    
    # Determine start plot for this OGU
    if (is.character(start_plot) && start_plot == "auto") {
      # Most environmentally extreme within OGU
      ogu_ndvi <- ndvi_std[ogu_plots]
      ogu_dem <- dem_std[ogu_plots]
      ogu_centroid <- c(mean(ogu_ndvi), mean(ogu_dem))
      distances_from_center <- sqrt((ogu_ndvi - ogu_centroid[1])^2 + 
                                      (ogu_dem - ogu_centroid[2])^2)
      start_idx_local <- which.max(distances_from_center)
      if (length(unique_ogus) == 1 || verbose) {
        msg(paste0("    Starting from Plot #", best_solution_sf$Plot_ID[ogu_plots[start_idx_local]], 
                   " (most extreme)"))
      }
      
    } else if (is.character(start_plot) && start_plot == "centroid") {
      # Closest to environmental centroid (most "average" plot)
      ogu_ndvi <- ndvi_std[ogu_plots]
      ogu_dem <- dem_std[ogu_plots]
      ogu_centroid <- c(mean(ogu_ndvi), mean(ogu_dem))
      distances_from_center <- sqrt((ogu_ndvi - ogu_centroid[1])^2 + 
                                      (ogu_dem - ogu_centroid[2])^2)
      start_idx_local <- which.min(distances_from_center)  # Min instead of max!
      if (length(unique_ogus) == 1 || verbose) {
        msg(paste0("    Starting from Plot #", best_solution_sf$Plot_ID[ogu_plots[start_idx_local]], 
                   " (closest to centroid)"))
      }
      
    } else if (start_plot_is_named_vector) {
      # Named vector: user specified start plot per OGU
      if (ogu_name %in% names(start_plot)) {
        specified_plot <- start_plot[ogu_name]
        if (specified_plot %in% best_solution_sf$Plot_ID[ogu_plots]) {
          start_idx_local <- which(best_solution_sf$Plot_ID[ogu_plots] == specified_plot)
          msg(paste0("    Starting from Plot #", specified_plot, " (user-specified)"))
        } else {
          warning("Plot ", specified_plot, " not found in OGU '", ogu_name, "'. Using first plot.")
          start_idx_local <- 1
        }
      } else {
        warning("No start plot specified for OGU '", ogu_name, "'. Using auto selection.")
        # Fall back to auto
        ogu_ndvi <- ndvi_std[ogu_plots]
        ogu_dem <- dem_std[ogu_plots]
        ogu_centroid <- c(mean(ogu_ndvi), mean(ogu_dem))
        distances_from_center <- sqrt((ogu_ndvi - ogu_centroid[1])^2 + 
                                        (ogu_dem - ogu_centroid[2])^2)
        start_idx_local <- which.max(distances_from_center)
      }
      
    } else if (is.numeric(start_plot)) {
      # Single integer: use if plot is in this OGU
      if (start_plot %in% best_solution_sf$Plot_ID[ogu_plots]) {
        start_idx_local <- which(best_solution_sf$Plot_ID[ogu_plots] == start_plot)
        msg(paste0("    Starting from Plot #", start_plot, " (user-specified)"))
      } else {
        if (length(unique_ogus) > 1) {
          # Multi-OGU: use auto for OGUs that don't contain the specified plot
          ogu_ndvi <- ndvi_std[ogu_plots]
          ogu_dem <- dem_std[ogu_plots]
          ogu_centroid <- c(mean(ogu_ndvi), mean(ogu_dem))
          distances_from_center <- sqrt((ogu_ndvi - ogu_centroid[1])^2 + 
                                          (ogu_dem - ogu_centroid[2])^2)
          start_idx_local <- which.max(distances_from_center)
        } else {
          warning("Plot ", start_plot, " not found. Using first plot.")
          start_idx_local <- 1
        }
      }
    } else {
      # Default to first plot
      start_idx_local <- 1
    }
    
    # Get ordering within OGU (local indices)
    local_order <- maximin_order(ogu_dist_matrix, start_idx_local)
    
    # Assign order numbers (1, 2, 3, ... within each OGU)
    sampling_order_within_ogu[ogu_plots[local_order]] <- 1:n_ogu_plots
  }
  
  # ============================================================================
  # SECTION 5: CREATE OUTPUT
  # ============================================================================
  
  msg("Creating output...")
  
  # Add OGU and sampling order to data
  best_solution_sf$OGU <- ogu_assignments
  best_solution_sf$sampling_order <- sampling_order_within_ogu
  
  # Create ordered sf object (order by OGU, then by sampling_order)
  ordered_sf <- best_solution_sf[order(best_solution_sf$OGU, best_solution_sf$sampling_order), ]
  
  # Summary statistics
  summary_df <- data.frame(
    n_plots = n_plots,
    n_ogus = length(unique_ogus),
    has_ogus = !is.null(ogus),
    has_dem = has_dem,
    ndvi_weight = ndvi_weight,
    dem_weight = if (has_dem) dem_weight else NA,
    mean_env_distance = round(mean(env_matrix[env_matrix > 0]), 3),
    stringsAsFactors = FALSE
  )
  
  msg("Done!")
  
  # ============================================================================
  # SECTION 6: SAVE OUTPUTS
  # ============================================================================
  
  msg("Saving outputs...")
  
  # Define output paths
  csv_field <- file.path(output_dir, paste0(output_prefix, "_field-ready.csv"))
  txt_summary <- file.path(output_dir, paste0(output_prefix, "_summary.txt"))
  
  # Transform to WGS84 for lat/lon
  ordered_wgs84 <- sf::st_transform(ordered_sf, 4326)
  coords_wgs84 <- sf::st_coordinates(ordered_wgs84)
  
  # Create field-ready CSV with everything in one file
  field_ready <- data.frame(
    Plot_ID = ordered_sf$Plot_ID,
    OGU = ordered_sf$OGU,
    Sampling_Order = ordered_sf$sampling_order,
    Latitude = round(coords_wgs84[, 2], 6),
    Longitude = round(coords_wgs84[, 1], 6),
    NDVI = round(ordered_sf$ndvi_mean, 3),
    Elevation_m = if(has_dem) round(ordered_sf$elevation, 1) else NA,
    stringsAsFactors = FALSE
  )
  
  # Sort by Plot_ID for easy lookup in field
  field_ready <- field_ready[order(field_ready$Plot_ID), ]
  
  utils::write.csv(field_ready, csv_field, row.names = FALSE)
  msg(paste0("  Saved: ", csv_field))
  
  # Save summary as text file
  sink(txt_summary)
  cat("SAMPLING ORDER OPTIMIZATION SUMMARY\n")
  cat("=====================================\n\n")
  cat(paste0("Number of plots: ", n_plots, "\n"))
  cat(paste0("Number of OGUs: ", length(unique_ogus), "\n"))
  if (!is.null(ogus)) {
    cat(paste0("OGU names: ", paste(unique_ogus, collapse = ", "), "\n"))
  }
  cat(paste0("DEM data used: ", ifelse(has_dem, "Yes", "No"), "\n"))
  cat(paste0("NDVI weight: ", ndvi_weight, "\n"))
  if (has_dem) cat(paste0("DEM weight: ", dem_weight, "\n"))
  cat(paste0("Mean environmental distance: ", round(mean(env_matrix[env_matrix > 0]), 3), "\n\n"))
  
  if (!is.null(ogus) && length(unique_ogus) > 1) {
    cat("OGU INFORMATION:\n")
    cat("-------------------\n")
    for (ogu_name in unique_ogus) {
      n_in_ogu <- sum(ogu_assignments == ogu_name)
      cat(paste0(ogu_name, ": ", n_in_ogu, " plots\n"))
    }
    cat("\n")
  }
  
  cat("FIELD-READY CSV CONTENTS:\n")
  cat("------------------------\n")
  cat("Plot_ID: Original plot identifier\n")
  if (!is.null(ogus)) {
    cat("OGU: Operational Geographic Unit (pre-defined study zone)\n")
  }
  cat("Sampling_Order: Optimized sequence within each OGU (1 = sample first)\n")
  cat("Latitude/Longitude: WGS84 coordinates for GPS\n")
  cat("NDVI: Vegetation index value\n")
  if (has_dem) cat("Elevation_m: Elevation in meters\n")
  cat("\n")
  
  cat("FIELD INSTRUCTIONS:\n")
  cat("------------------\n")
  if (is.null(ogus) || length(unique_ogus) == 1) {
    cat("1. Open the field-ready CSV on your device\n")
    cat("2. Sort by Sampling_Order column\n")
    cat("3. Navigate to the Lat/Lon of the first plot\n")
    cat("4. Sample that plot, then move to the next in order\n")
    cat("\nThis order maximizes environmental diversity capture.\n")
    cat("If you run out of time/budget, you've still sampled across the full gradient.\n")
  } else {
    cat("1. Open the field-ready CSV on your device\n")
    cat("2. Choose which OGU to sample first (based on accessibility, weather, permits)\n")
    cat("3. Filter to that OGU and sort by Sampling_Order\n")
    cat("4. Navigate to plots using Lat/Lon in order (1, 2, 3, ...)\n")
    cat("5. When OGU complete, move to next OGU\n")
    cat("\nEach OGU is ordered independently to maximize diversity within that zone.\n")
    cat("OGUs are typically defined based on ecological zones, accessibility,\n")
    cat("management units, or other a priori knowledge of the study area.\n")
  }
  
  cat("\nOUTPUT FILES:\n")
  cat("-------------\n")
  cat(paste0(basename(csv_field), " - Complete field sheet (TAKE THIS TO FIELD!)\n"))
  cat(paste0(basename(txt_summary), " - This summary file\n"))
  sink()
  msg(paste0("  Saved: ", txt_summary))
  
  msg(paste0("All files saved to: ", output_dir))
  
  # ============================================================================
  # SECTION 7: RETURN
  # ============================================================================
  
  return(list(
    ordered_plots = ordered_sf,
    field_ready = field_ready,
    env_distance_matrix = env_matrix,
    summary = summary_df
  ))
}