#' @title Sampling Order Optimization
#'
#' @description
#' Generates an optimized sampling order for field plots based on environmental 
#' diversity (NDVI, DEM) to maximize species accumulation efficiency. Can operate 
#' with or without spatial clustering for logistical efficiency in large study areas.
#'
#' @param best_solution_sf An sf object from sampleboost() containing plot locations 
#'   and environmental data (must have columns: ndvi_mean, and optionally elevation)
#' @param n_clusters Integer or "auto". NULL = no clustering (pure diversity ordering).
#'   Integer = exact number of clusters. "auto" = suggest 2, 3, and 4 cluster scenarios.
#' @param method Character. Clustering method: "spatial" (k-means on coordinates), 
#'   "elevation" (elevation-based zones, requires DEM), or "auto" (try both if DEM available)
#' @param ndvi_weight Numeric. Weight for NDVI in environmental distance calculation (default = 1)
#' @param dem_weight Numeric. Weight for elevation in environmental distance calculation (default = 0.5)
#' @param start_plot Integer or "auto". Plot ID to start from (default = 1), or "auto" 
#'   to automatically select most environmentally extreme plot
#' @param verbose Logical. Print progress messages (default = TRUE)
#' @param output_dir Character. Directory for output files (default = working directory)
#' @param output_prefix Character. Prefix for output filenames (default = "SamplingOrder")
#'
#' @return A list containing:
#' \itemize{
#'   \item{ordered_plots}{ sf object with sampling_order and cluster columns added}
#'   \item{env_distance_matrix}{ Environmental distance matrix between all plots}
#'   \item{summary}{ Data frame with ordering statistics}
#' }
#' 
#' Additionally, two files are automatically saved to output_dir:
#' \itemize{
#'   \item{[prefix]_field-ready.csv}{ Complete field sheet with Plot_ID, Sampling_Order, 
#'   coordinates (lat/lon), NDVI, elevation, and cluster info - everything in one file}
#'   \item{[prefix]_summary.txt}{ Summary statistics and field instructions}
#' }
#' 
#' If n_clusters = "auto", returns a list of scenarios instead of a single result,
#' and no files are saved.
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
#' **Clustering Options**:
#' \itemize{
#'   \item{No clustering}{ Pure diversity-based ordering across all plots}
#'   \item{Spatial clustering}{ k-means on X,Y coordinates to minimize travel distance}
#'   \item{Elevation clustering}{ Natural breaks in elevation for accessibility}
#' }
#' 
#' Within each cluster, maximin diversity ordering is applied. User selects which 
#' cluster to sample first based on field logistics.
#'
#' @importFrom sf st_transform st_crs st_coordinates st_as_sf st_drop_geometry st_geometry st_sf
#' @importFrom stats kmeans dist quantile
#' @importFrom utils write.csv
#' @export
#'
#' @examples
#' \dontrun{
#' # Run sampleboost first
#' result <- sampleboost(ndvi, ignorance, site, nplot = 50, plot_radius = 20, perm = 100)
#' 
#' # Simple diversity ordering - creates field-ready CSV
#' order1 <- sampling_order(result$best_solution_sf)
#' 
#' # With spatial clustering (3 clusters)
#' order2 <- sampling_order(result$best_solution_sf, n_clusters = 3, method = "spatial")
#' 
#' # Suggest multiple scenarios (no files saved - for comparison only)
#' scenarios <- sampling_order(result$best_solution_sf, n_clusters = "auto")
#' 
#' # Pick best and save final outputs
#' final <- sampling_order(
#'   result$best_solution_sf, 
#'   n_clusters = 3,
#'   output_dir = "~/fieldwork",
#'   output_prefix = "Final_Order"
#' )
#' 
#' # With elevation-based clustering (requires DEM in original sampleboost)
#' order3 <- sampling_order(result$best_solution_sf, n_clusters = 3, method = "elevation")
#' }

sampling_order <- function(best_solution_sf, 
                           n_clusters = NULL, 
                           method = "auto",
                           ndvi_weight = 1, 
                           dem_weight = 0.5,
                           start_plot = 1,
                           verbose = TRUE,
                           output_dir = getwd(),
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
  
  if (method == "elevation" && !has_dem) {
    stop("method='elevation' requires DEM data (elevation column in best_solution_sf)")
  }
  
  n_plots <- nrow(best_solution_sf)
  
  # Create Plot_ID if it doesn't exist
  if (!"Plot_ID" %in% names(best_solution_sf)) {
    best_solution_sf$Plot_ID <- 1:n_plots
    msg("  Created Plot_ID column (1 to ", n_plots, ")")
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
  # SECTION 3: HANDLE AUTO-SCENARIOS
  # ============================================================================
  
  if (!is.null(n_clusters) && n_clusters == "auto") {
    msg("Generating multiple clustering scenarios...")
    
    scenarios <- list()
    
    # No clustering
    scenarios$no_clusters <- sampling_order(
      best_solution_sf, n_clusters = NULL, method = method,
      ndvi_weight = ndvi_weight, dem_weight = dem_weight,
      start_plot = start_plot, verbose = FALSE,
      output_dir = tempdir(), output_prefix = "temp"
    )
    
    # 2, 3, 4 spatial clusters
    for (k in 2:4) {
      scenario_name <- paste0("clusters_", k)
      scenarios[[scenario_name]] <- sampling_order(
        best_solution_sf, n_clusters = k, method = "spatial",
        ndvi_weight = ndvi_weight, dem_weight = dem_weight,
        start_plot = start_plot, verbose = FALSE,
        output_dir = tempdir(), output_prefix = "temp"
      )
    }
    
    # Elevation-based if DEM available
    if (has_dem) {
      scenarios$elevation_zones <- sampling_order(
        best_solution_sf, n_clusters = 3, method = "elevation",
        ndvi_weight = ndvi_weight, dem_weight = dem_weight,
        start_plot = start_plot, verbose = FALSE,
        output_dir = tempdir(), output_prefix = "temp"
      )
    }
    
    msg(paste0("Generated ", length(scenarios), " scenarios: ", paste(names(scenarios), collapse = ", ")))
    msg("Note: Review scenarios and re-run with chosen n_clusters to save outputs")
    return(scenarios)
  }
  
  # ============================================================================
  # SECTION 4: CLUSTERING (IF REQUESTED)
  # ============================================================================
  
  cluster_assignments <- rep(1, n_plots)  # Default: all in one cluster
  
  if (!is.null(n_clusters) && n_clusters > 1) {
    msg(paste0("Creating ", n_clusters, " clusters using method: ", method, "..."))
    
    if (method == "auto") {
      method <- if (has_dem) "spatial" else "spatial"  # Default to spatial
    }
    
    if (method == "spatial") {
      # K-means clustering on coordinates
      kmeans_result <- stats::kmeans(coords, centers = n_clusters, nstart = 25)
      cluster_assignments <- kmeans_result$cluster
      msg(paste0("  Spatial k-means clustering complete"))
      
    } else if (method == "elevation") {
      # Elevation-based clustering using quantiles
      elevation_breaks <- stats::quantile(dem_data, probs = seq(0, 1, length.out = n_clusters + 1))
      cluster_assignments <- cut(dem_data, breaks = elevation_breaks, 
                                 labels = FALSE, include.lowest = TRUE)
      msg(paste0("  Elevation-based clustering complete"))
      msg(paste0("  Elevation ranges: ", 
                 paste(round(elevation_breaks), collapse = " -> "), " m"))
    }
    
    # Report cluster sizes
    cluster_sizes <- table(cluster_assignments)
    msg(paste0("  Cluster sizes: ", paste(cluster_sizes, collapse = ", ")))
  }
  
  # ============================================================================
  # SECTION 5: APPLY MAXIMIN ORDERING
  # ============================================================================
  
  msg("Applying maximin diversity ordering...")
  
  if (is.null(n_clusters) || n_clusters == 1) {
    # No clustering - global ordering
    
    # Determine start plot
    if (start_plot == "auto") {
      # Select most environmentally extreme plot (farthest from centroid)
      env_centroid <- c(mean(ndvi_std), mean(dem_std))
      distances_from_center <- sqrt((ndvi_std - env_centroid[1])^2 + 
                                      (dem_std - env_centroid[2])^2)
      start_idx <- which.max(distances_from_center)
      msg(paste0("  Auto-selected start plot: #", start_idx, " (most environmentally extreme)"))
    } else {
      start_idx <- start_plot
      msg(paste0("  Starting from plot #", start_idx))
    }
    
    sampling_order_seq <- maximin_order(env_matrix, start_idx)
    
  } else {
    # Clustered ordering
    sampling_order_seq <- integer(n_plots)
    current_position <- 1
    
    for (cluster_id in sort(unique(cluster_assignments))) {
      cluster_plots <- which(cluster_assignments == cluster_id)
      n_cluster_plots <- length(cluster_plots)
      
      msg(paste0("  Ordering cluster ", cluster_id, " (", n_cluster_plots, " plots)..."))
      
      # Extract submatrix for this cluster
      cluster_dist_matrix <- env_matrix[cluster_plots, cluster_plots]
      
      # Determine start plot within cluster
      if (start_plot == "auto") {
        # Most extreme within cluster
        cluster_ndvi <- ndvi_std[cluster_plots]
        cluster_dem <- dem_std[cluster_plots]
        cluster_centroid <- c(mean(cluster_ndvi), mean(cluster_dem))
        distances_from_center <- sqrt((cluster_ndvi - cluster_centroid[1])^2 + 
                                        (cluster_dem - cluster_centroid[2])^2)
        start_idx_local <- which.max(distances_from_center)
      } else {
        start_idx_local <- 1  # Start from first plot in cluster
      }
      
      # Get ordering within cluster (local indices)
      local_order <- maximin_order(cluster_dist_matrix, start_idx_local)
      
      # Convert to global indices and assign sequential order numbers
      global_indices <- cluster_plots[local_order]
      sampling_order_seq[current_position:(current_position + n_cluster_plots - 1)] <- global_indices
      current_position <- current_position + n_cluster_plots
    }
  }
  
  # ============================================================================
  # SECTION 6: CREATE OUTPUT
  # ============================================================================
  
  msg("Creating output...")
  
  # Preserve geometry before manipulating data
  geom_backup <- sf::st_geometry(best_solution_sf)
  
  # Create ordered data frame (drop geometry temporarily)
  ordered_df <- sf::st_drop_geometry(best_solution_sf)
  ordered_df$sampling_order <- match(1:n_plots, sampling_order_seq)
  ordered_df$cluster <- as.factor(cluster_assignments)
  
  # Get the order indices
  order_indices <- order(ordered_df$sampling_order)
  
  # Reorder both data and geometry consistently
  ordered_df <- ordered_df[order_indices, ]
  geom_ordered <- geom_backup[order_indices]
  
  # Reconstruct sf object with ordered geometry
  ordered_sf <- sf::st_sf(ordered_df, geometry = geom_ordered, crs = sf::st_crs(best_solution_sf))
  
  # Summary statistics
  summary_df <- data.frame(
    n_plots = n_plots,
    n_clusters = if (is.null(n_clusters)) 1 else n_clusters,
    method = method,
    has_dem = has_dem,
    ndvi_weight = ndvi_weight,
    dem_weight = if (has_dem) dem_weight else NA,
    mean_env_distance = round(mean(env_matrix[env_matrix > 0]), 3),
    stringsAsFactors = FALSE
  )
  
  msg("Done!")
  
  # ============================================================================
  # SECTION 7: SAVE OUTPUTS
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
    Plot_ID = ordered_df$Plot_ID,
    Sampling_Order = ordered_df$sampling_order,
    Cluster = ordered_df$cluster,
    Latitude = round(coords_wgs84[, 2], 6),
    Longitude = round(coords_wgs84[, 1], 6),
    NDVI = round(ordered_df$ndvi_mean, 3),
    Elevation_m = if(has_dem) round(ordered_df$elevation, 1) else NA,
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
  cat(paste0("Number of clusters: ", if (is.null(n_clusters)) 1 else n_clusters, "\n"))
  cat(paste0("Clustering method: ", method, "\n"))
  cat(paste0("DEM data used: ", ifelse(has_dem, "Yes", "No"), "\n"))
  cat(paste0("NDVI weight: ", ndvi_weight, "\n"))
  if (has_dem) cat(paste0("DEM weight: ", dem_weight, "\n"))
  cat(paste0("Mean environmental distance: ", round(mean(env_matrix[env_matrix > 0]), 3), "\n\n"))
  
  if (!is.null(n_clusters) && n_clusters > 1) {
    cat("CLUSTER INFORMATION:\n")
    cat("-------------------\n")
    for (cl in sort(unique(cluster_assignments))) {
      n_in_cluster <- sum(cluster_assignments == cl)
      cat(paste0("Cluster ", cl, ": ", n_in_cluster, " plots\n"))
    }
    cat("\n")
  }
  
  cat("FIELD-READY CSV CONTENTS:\n")
  cat("------------------------\n")
  cat("Plot_ID: Original plot identifier\n")
  cat("Sampling_Order: Optimized sequence (1 = sample first)\n")
  if (!is.null(n_clusters) && n_clusters > 1) {
    cat("Cluster: Spatial group (choose which to sample first)\n")
  }
  cat("Latitude/Longitude: WGS84 coordinates for GPS\n")
  cat("NDVI: Vegetation index value\n")
  if (has_dem) cat("Elevation_m: Elevation in meters\n")
  cat("\n")
  
  cat("FIELD INSTRUCTIONS:\n")
  cat("------------------\n")
  if (is.null(n_clusters) || n_clusters == 1) {
    cat("1. Open the field-ready CSV on your device\n")
    cat("2. Sort by Sampling_Order column\n")
    cat("3. Navigate to the Lat/Lon of the first plot\n")
    cat("4. Sample that plot, then move to the next in order\n")
    cat("\nThis order maximizes environmental diversity capture.\n")
    cat("If you run out of time/budget, you've still sampled across the full gradient.\n")
  } else {
    cat("1. Open the field-ready CSV on your device\n")
    cat("2. Choose which Cluster to sample first (based on accessibility)\n")
    cat("3. Filter to that cluster and sort by Sampling_Order\n")
    cat("4. Navigate to plots using Lat/Lon in order\n")
    cat("5. When cluster complete, move to next cluster\n")
    cat("\nEach cluster is ordered to maximize diversity within that spatial group.\n")
  }
  
  cat("\nOUTPUT FILES:\n")
  cat("-------------\n")
  cat(paste0(basename(csv_field), " - Complete field sheet (TAKE THIS TO FIELD!)\n"))
  cat(paste0(basename(txt_summary), " - This summary file\n"))
  sink()
  msg(paste0("  Saved: ", txt_summary))
  
  msg(paste0("All files saved to: ", output_dir))
  
  # ============================================================================
  # SECTION 8: RETURN
  # ============================================================================
  
  return(list(
    ordered_plots = ordered_sf,
    field_ready = field_ready,
    env_distance_matrix = env_matrix,
    summary = summary_df
  ))
}