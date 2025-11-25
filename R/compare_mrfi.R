#' @title Compare Two Maps of Relative Floristic Ignorance (Incremental Approach)
#'
#' @description
#' Compares initial and updated floristic surveys using incremental MRFI calculation.
#' The updated MRFI is calculated as an UPDATE of the baseline (maintaining fixed
#' maximum ignorance) rather than independent recalculation, ensuring comparable
#' scales even when new taxa are discovered.
#' 
#' Visualization focuses on KNOWLEDGE GAINED:
#' \itemize{
#'   \item{Page 1: Initial MRFI + Knowledge Gained (IRFI improvement)}
#'   \item{Page 2: Initial Richness + Richness Gain}
#'   \item{Page 3: Richness gain by IRFI category + Improvement vs Initial IRFI}
#'   \item{Page 4: Summary statistics}
#' }
#' 
#' This function implements the evaluation phase (Q3) of the Next Generation
#' Floristics workflow (D'Antraccoli et al. 2022).
#'
#' @param data_flor_initial Data frame with initial floristic data. Must have
#'   columns: Taxon, Long, Lat, uncertainty, year.
#' @param data_flor_updated Data frame with updated/complete floristic data.
#'   Must have columns: Taxon, Long, Lat, uncertainty, year.
#' @param site An `sf` object representing the study area.
#' @param tau Numeric. Temporal decay parameter (% taxa loss per 100 years).
#' @param cellsize Numeric. Raster cell size in meters.
#' @param year_initial Numeric. Year of initial survey (default: extracted from data).
#' @param year_updated Numeric. Year of updated survey (default: extracted from data).
#' @param excl_areas Optional `sf` object for exclusion areas.
#' @param CRS.new Numeric EPSG code for projected CRS (default = 3035).
#' @param target_quantile Numeric. Quantile for IRFI categorization (default = 0.75).
#' @param output_dir Directory for output files (default = working directory + "/output").
#' @param output_prefix Filename prefix (default = "MRFI_comparison").
#' @param verbose Logical. Print progress messages (default = TRUE).
#' @param site_buffer Logical. Expand site boundary for analysis (default = FALSE).
#' @param buffer_width Numeric. Buffer distance in meters if site_buffer = TRUE.
#' @param mask_method Character. Masking method: "touches" or "none" (default = "touches").
#' @param use_coverage_weighting Logical. Use coverage-weighted rasterization (default = TRUE).
#' @param ... Additional parameters passed to \code{ignorance_map()}.
#'
#' @return A list with class "mrfi_comparison":
#' \itemize{
#'   \item{\code{mrfi_initial}}{Initial MRFI output (ignorance_map object)}
#'   \item{\code{mrfi_updated}}{Updated MRFI output (incremental calculation)}
#'   \item{\code{improvement_raster}}{SpatRaster of IRFI improvements}
#'   \item{\code{richness_gain_raster}}{SpatRaster of richness gains}
#'   \item{\code{summary}}{Data frame with comparison statistics}
#'   \item{\code{richness_comparison}}{Data frame with richness statistics}
#'   \item{\code{targeting_efficiency}}{Data frame with targeting metrics}
#' }
#'
#' @details
#' **Incremental MRFI Calculation:**
#' 
#' The updated MRFI is calculated as:
#' \code{IRFI_updated = max_ignorance_T1 - Σ(aged_baseline + new_observations)}
#' 
#' Where:
#' \itemize{
#'   \item{aged_baseline: T1 observations with temporal scores updated to T2}
#'   \item{new_observations: New observations between T1 and T2}
#'   \item{max_ignorance_T1: Fixed reference from initial survey (constant)}
#' }
#' 
#' This maintains a constant baseline (max_ignorance from initial survey),
#' preventing the paradox where discovering new taxa appears to increase ignorance.
#' 
#' **Negative IRFI values:** Can occur in well-sampled cells where new observations
#' continue to reduce ignorance below the initial reference. These are shown
#' accurately in the improvement map.
#'
#' **PDF Output (4 pages, A4 landscape):**
#' \itemize{
#'   \item{Page 1: Initial MRFI (left) + Knowledge Gained map (right)}
#'   \item{Page 2: Initial Richness (left) + Richness Gain map (right)}
#'   \item{Page 3: Richness by IRFI boxplot (left) + Improvement scatterplot (right)}
#'   \item{Page 4: Summary statistics table}
#' }
#'
#' @importFrom terra rast values compareGeom writeRaster global res ext as.polygons crs rasterize vect
#' @importFrom sf st_as_sf st_crs st_transform st_buffer st_bbox st_coordinates st_geometry st_make_valid st_intersection st_intersects st_difference st_union
#' @importFrom stats cor quantile median
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom ggplot2 ggplot aes geom_tile geom_histogram geom_vline geom_point geom_smooth geom_boxplot scale_fill_gradientn scale_fill_gradient2 scale_fill_gradient coord_equal theme_classic labs ggtitle element_text theme geom_sf
#' @importFrom gridExtra grid.arrange tableGrob
#' @importFrom grDevices pdf dev.off
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid textGrob gpar unit
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare two survey periods
#' comparison <- compare_mrfi(
#'   data_flor_initial = floratus_2020,
#'   data_flor_updated = floratus_2023,
#'   site = park,
#'   tau = 30,
#'   cellsize = 2000
#' )
#' 
#' # View results
#' print(comparison$summary)
#' terra::plot(comparison$improvement_raster)
#' 
#' # Evaluate if sampling was effective (Q3 from NGF workflow)
#' mean_improvement_pct <- comparison$summary$Value[
#'   comparison$summary$Statistic == "Mean IRFI improvement (%)"
#' ]
#' 
#' if (as.numeric(mean_improvement_pct) > 20) {
#'   message("Substantial improvement! Inventory may be sufficient.")
#' } else {
#'   message("Limited improvement. Additional sampling recommended.")
#' }
#' }

compare_mrfi <- function(
    data_flor_initial,
    data_flor_updated,
    site,
    tau,
    cellsize,
    year_initial = NULL,
    year_updated = NULL,
    excl_areas = NULL,
    CRS.new = 3035,
    target_quantile = 0.75,
    output_dir = file.path(getwd(), "output"),
    output_prefix = "MRFI_comparison",
    verbose = TRUE,
    site_buffer = FALSE,
    buffer_width = NULL,
    mask_method = "touches",
    use_coverage_weighting = TRUE,
    ...
) {
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  # ============================================================================
  # SECTION 1: INPUT VALIDATION
  # ============================================================================
  
  msg("Validating inputs...")
  
  if (missing(data_flor_initial) || is.null(data_flor_initial)) {
    stop("data_flor_initial is required")
  }
  
  if (missing(data_flor_updated) || is.null(data_flor_updated)) {
    stop("data_flor_updated is required")
  }
  
  if (missing(site) || is.null(site)) {
    stop("site is required")
  }
  
  if (missing(tau) || is.null(tau)) {
    stop("tau is required")
  }
  
  if (missing(cellsize) || is.null(cellsize)) {
    stop("cellsize is required")
  }
  
  if (!is.numeric(target_quantile) || target_quantile < 0 || target_quantile > 1) {
    stop("target_quantile must be between 0 and 1")
  }
  
  # Set default years if not provided
  if (is.null(year_initial)) {
    year_initial <- min(data_flor_initial$year, na.rm = TRUE)
    msg(paste0("year_initial not specified, using earliest year from data: ", year_initial))
  }
  
  if (is.null(year_updated)) {
    year_updated <- max(data_flor_updated$year, na.rm = TRUE)
    msg(paste0("year_updated not specified, using latest year from data: ", year_updated))
  }
  
  # Validate year order
  if (year_updated <= year_initial) {
    stop("year_updated must be greater than year_initial")
  }
  
  # ============================================================================
  # SECTION 2: CALCULATE BASELINE MRFI
  # ============================================================================
  
  msg("Calculating baseline MRFI (initial survey)...")
  
  mrfi_initial <- ignorance_map(
    data_flor = data_flor_initial,
    site = site,
    year_study = year_initial,
    excl_areas = excl_areas,
    CRS.new = CRS.new,
    tau = tau,
    cellsize = cellsize,
    verbose = verbose,
    check_overlap = FALSE,
    output_dir = output_dir,
    output_prefix = paste0(output_prefix, "_initial"),
    site_buffer = site_buffer,
    buffer_width = buffer_width,
    mask_method = mask_method,
    use_coverage_weighting = use_coverage_weighting,
    ...
  )
  
  msg("  Baseline MRFI calculated successfully.")
  
  # ============================================================================
  # SECTION 3: IDENTIFY NEW RECORDS
  # ============================================================================
  
  msg("Identifying new records...")
  
  data_flor_new <- identify_new_records(
    data_initial = data_flor_initial,
    data_updated = data_flor_updated
  )
  
  msg(paste0("  Found ", nrow(data_flor_new), " new observations"))
  msg(paste0("  (", nrow(data_flor_updated) - nrow(data_flor_initial), 
             " additional records in updated dataset)"))
  
  # ============================================================================
  # SECTION 4: CALCULATE INCREMENTAL UPDATED MRFI
  # ============================================================================
  
  msg("Calculating updated MRFI (incremental approach)...")
  msg("  Updated MRFI uses baseline as reference (fixed max_ignorance)...")
  
  mrfi_updated <- calculate_incremental_mrfi_internal(
    data_flor_initial = data_flor_initial,  # Need this for aging baseline taxa
    data_flor_new = data_flor_new,          # Only new observations
    mrfi_baseline = mrfi_initial,
    year_updated = year_updated,
    site = site,
    excl_areas = excl_areas,
    CRS.new = CRS.new,
    tau = tau,
    cellsize = cellsize,
    verbose = verbose,
    site_buffer = site_buffer,
    buffer_width = buffer_width,
    mask_method = mask_method,
    use_coverage_weighting = use_coverage_weighting,
    ...
  )
  
  msg("  Incremental MRFI calculated successfully.")
  
  # ============================================================================
  # SECTION 5: EXTRACT RASTERS AND CALCULATE IMPROVEMENTS
  # ============================================================================
  
  mrfi_initial_raster <- mrfi_initial$MRFI
  mrfi_updated_raster <- mrfi_updated$MRFI
  rich_initial_raster <- mrfi_initial$RICH
  rich_updated_raster <- mrfi_updated$RICH
  
  # Check compatibility
  msg("Checking raster compatibility...")
  if (!terra::compareGeom(mrfi_initial_raster, mrfi_updated_raster, stopOnError = FALSE)) {
    stop("Initial and updated MRFI rasters have incompatible geometries.")
  }
  
  # Calculate improvements
  msg("Calculating improvements...")
  improvement_raster <- mrfi_initial_raster - mrfi_updated_raster
  richness_gain_raster <- rich_updated_raster - rich_initial_raster
  
  # Extract values
  initial_vals <- terra::values(mrfi_initial_raster, na.rm = FALSE)
  updated_vals <- terra::values(mrfi_updated_raster, na.rm = FALSE)
  improvement_vals <- terra::values(improvement_raster, na.rm = FALSE)
  
  rich_initial_vals <- terra::values(rich_initial_raster, na.rm = FALSE)
  rich_updated_vals <- terra::values(rich_updated_raster, na.rm = FALSE)
  richness_gain_vals <- terra::values(richness_gain_raster, na.rm = FALSE)
  
  # Create cell-level data frame
  cell_data <- data.frame(
    cell_id = seq_along(initial_vals),
    irfi_initial = as.vector(initial_vals),
    irfi_updated = as.vector(updated_vals),
    irfi_improvement = as.vector(improvement_vals),
    rich_initial = as.vector(rich_initial_vals),
    rich_updated = as.vector(rich_updated_vals),
    richness_gain = as.vector(richness_gain_vals),
    stringsAsFactors = FALSE
  )
  
  cell_data_valid <- cell_data[complete.cases(cell_data), ]
  
  if (nrow(cell_data_valid) == 0) {
    stop("No valid cells found for comparison.")
  }
  
  msg(paste0("Analyzing ", nrow(cell_data_valid), " valid cells..."))
  
  cell_data_valid$irfi_improvement_pct <- (cell_data_valid$irfi_improvement / 
                                             cell_data_valid$irfi_initial) * 100
  cell_data_valid$irfi_improvement_pct[is.infinite(cell_data_valid$irfi_improvement_pct)] <- NA
  
  # ============================================================================
  # SECTION 6: SUMMARY STATISTICS
  # ============================================================================
  
  msg("Computing summary statistics...")
  
  # IRFI statistics
  mean_initial <- mean(cell_data_valid$irfi_initial, na.rm = TRUE)
  mean_updated <- mean(cell_data_valid$irfi_updated, na.rm = TRUE)
  mean_improvement <- mean(cell_data_valid$irfi_improvement, na.rm = TRUE)
  mean_improvement_pct <- (mean_improvement / mean_initial) * 100
  
  median_initial <- median(cell_data_valid$irfi_initial, na.rm = TRUE)
  median_updated <- median(cell_data_valid$irfi_updated, na.rm = TRUE)
  median_improvement <- median(cell_data_valid$irfi_improvement, na.rm = TRUE)
  median_improvement_pct <- (median_improvement / median_initial) * 100
  
  max_improvement <- max(cell_data_valid$irfi_improvement, na.rm = TRUE)
  max_worsening <- min(cell_data_valid$irfi_improvement, na.rm = TRUE)
  
  n_improved <- sum(cell_data_valid$irfi_improvement > 0, na.rm = TRUE)
  n_worsened <- sum(cell_data_valid$irfi_improvement < 0, na.rm = TRUE)
  n_unchanged <- sum(cell_data_valid$irfi_improvement == 0, na.rm = TRUE)
  
  pct_improved <- (n_improved / nrow(cell_data_valid)) * 100
  pct_worsened <- (n_worsened / nrow(cell_data_valid)) * 100
  
  # Richness statistics
  mean_rich_initial <- mean(cell_data_valid$rich_initial, na.rm = TRUE)
  mean_rich_updated <- mean(cell_data_valid$rich_updated, na.rm = TRUE)
  mean_richness_gain <- mean(cell_data_valid$richness_gain, na.rm = TRUE)
  mean_richness_gain_pct <- (mean_richness_gain / mean_rich_initial) * 100
  
  total_rich_initial <- max(cell_data_valid$rich_initial, na.rm = TRUE)
  total_rich_updated <- max(cell_data_valid$rich_updated, na.rm = TRUE)
  total_taxa_gain <- total_rich_updated - total_rich_initial
  
  # Richness comparison table
  richness_stats <- data.frame(
    Statistic = c(
      "Mean richness (initial)",
      "Mean richness (updated)",
      "Mean richness gain",
      "Mean richness gain (%)",
      "Total taxa discovered"
    ),
    Value = c(
      round(mean_rich_initial, 1),
      round(mean_rich_updated, 1),
      round(mean_richness_gain, 1),
      round(mean_richness_gain_pct, 1),
      round(total_taxa_gain, 0)
    ),
    stringsAsFactors = FALSE
  )
  
  # ============================================================================
  # SECTION 7: TARGETING EFFICIENCY
  # ============================================================================
  
  msg("Evaluating targeting efficiency...")
  
  cor_ignorance_improvement <- cor(
    cell_data_valid$irfi_initial,
    cell_data_valid$irfi_improvement,
    use = "complete.obs"
  )
  
  ignorance_threshold <- quantile(
    cell_data_valid$irfi_initial,
    probs = target_quantile,
    na.rm = TRUE
  )
  
  most_ignorant_cells <- cell_data_valid$irfi_initial >= ignorance_threshold
  
  improvement_in_target <- mean(
    cell_data_valid$irfi_improvement[most_ignorant_cells],
    na.rm = TRUE
  )
  
  improvement_overall <- mean(cell_data_valid$irfi_improvement, na.rm = TRUE)
  targeting_ratio <- improvement_in_target / improvement_overall
  targeting_score <- min(targeting_ratio * 50, 100)
  
  targeting_stats <- data.frame(
    Statistic = c(
      "Correlation (initial IRFI vs improvement)",
      paste0("Improvement in top ", (1-target_quantile)*100, "% ignorant cells"),
      "Targeting efficiency score (0-100)"
    ),
    Value = c(
      round(cor_ignorance_improvement, 3),
      round(improvement_in_target, 2),
      round(targeting_score, 1)
    ),
    stringsAsFactors = FALSE
  )
  
  # ============================================================================
  # SECTION 8: BUILD SUMMARY TABLE
  # ============================================================================
  
  end_time <- Sys.time()
  
  summary_df <- data.frame(
    Statistic = c(
      "Analysis timestamp",
      "Elapsed time (secs)",
      "Comparison method",
      "Valid cells analyzed",
      "Mean IRFI (initial)",
      "Mean IRFI (updated)",
      "Mean IRFI improvement",
      "Mean IRFI improvement (%)",
      "Median IRFI improvement",
      "Cells improved (%)",
      "Cells worsened (%)",
      "Mean richness (initial)",
      "Mean richness (updated)",
      "Mean richness gain",
      "Total taxa discovered",
      "Correlation (IRFI vs improvement)",
      "Targeting efficiency score"
    ),
    Value = c(
      as.character(end_time),
      round(as.numeric(difftime(end_time, start_time, units = "secs")), 2),
      "Incremental (fixed baseline)",
      nrow(cell_data_valid),
      round(mean_initial, 2),
      round(mean_updated, 2),
      round(mean_improvement, 2),
      round(mean_improvement_pct, 1),
      round(median_improvement, 2),
      round(pct_improved, 1),
      round(pct_worsened, 1),
      round(mean_rich_initial, 1),
      round(mean_rich_updated, 1),
      round(mean_richness_gain, 1),
      round(total_taxa_gain, 0),
      round(cor_ignorance_improvement, 3),
      round(targeting_score, 1)
    )
  )
  
  # ============================================================================
  # SECTION 9: PREPARE SITE BOUNDARIES FOR PLOTTING
  # ============================================================================
  
  msg("Preparing visualization components...")
  
  # Extract site geometry from baseline (already projected)
  site_proj_original <- sf::st_as_sf(site)
  crs_sf <- sf::st_crs(CRS.new)
  
  if (is.na(sf::st_crs(site_proj_original))) {
    sf::st_crs(site_proj_original) <- 4326
  }
  
  site_proj_original <- sf::st_transform(sf::st_make_valid(site_proj_original), crs_sf)
  
  # ============================================================================
  # SECTION 10: GENERATE PLOTS
  # ============================================================================
  
  msg("Generating comparison plots...")
  
  # Convert rasters to data frames
  initial_mrfi_df <- as.data.frame(mrfi_initial_raster, xy = TRUE)
  names(initial_mrfi_df) <- c("x", "y", "value")
  
  improvement_df <- as.data.frame(improvement_raster, xy = TRUE)
  names(improvement_df) <- c("x", "y", "value")
  
  initial_rich_df <- as.data.frame(rich_initial_raster, xy = TRUE)
  names(initial_rich_df) <- c("x", "y", "value")
  
  richness_gain_df <- as.data.frame(richness_gain_raster, xy = TRUE)
  names(richness_gain_df) <- c("x", "y", "value")
  
  # Calculate scales
  global_max_irfi <- terra::global(mrfi_initial_raster, "max", na.rm = TRUE)$max
  irfi_breaks <- seq(0, global_max_irfi, length.out = 9)
  
  # Use 95th percentile of absolute values (robust, excludes extreme outliers)
  improvement_vals <- terra::values(improvement_raster, na.rm = TRUE)
  improvement_limit <- quantile(abs(improvement_vals), probs = 0.95, na.rm = TRUE)
  improvement_breaks <- seq(-improvement_limit, improvement_limit, length.out = 9)
  
  richness_max_initial <- terra::global(rich_initial_raster, "max", na.rm = TRUE)$max
  richness_breaks <- seq(0, richness_max_initial, length.out = 9)
  
  richness_gain_max <- terra::global(richness_gain_raster, "max", na.rm = TRUE)$max
  richness_gain_breaks <- seq(0, richness_gain_max, length.out = 9)
  
  # ---------------------------------------------------------------------
  # PAGE 1 PLOTS: INITIAL MRFI + KNOWLEDGE GAINED
  # ---------------------------------------------------------------------
  
  p1_initial_mrfi <- ggplot2::ggplot(initial_mrfi_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical',
                   legend.key.width = grid::unit(0.6, "cm"),
                   plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
      limits = c(0, global_max_irfi),
      breaks = irfi_breaks,
      labels = round(irfi_breaks, 1),
      na.value = "transparent",
      guide = ggplot2::guide_legend(title = "IRFI", keyheight = grid::unit(1.2, "lines"),
                                    keywidth = grid::unit(1.2, "lines"),
                                    label.theme = ggplot2::element_text(size = 9),
                                    title.theme = ggplot2::element_text(size = 10, face = "bold"))
    ) +
    ggplot2::geom_sf(data = site_proj_original, fill = NA, color = "black",
                     linewidth = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("Initial Survey - MRFI") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  p1_improvement <- ggplot2::ggplot(improvement_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical',
                   legend.key.width = grid::unit(0.6, "cm"),
                   plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradient2(
      low = "#D73027",
      mid = "white",
      high = "#1A9850",
      midpoint = 0,
      limits = c(-improvement_limit, improvement_limit),
      breaks = improvement_breaks,
      labels = round(improvement_breaks, 1),
      name = "ΔIRFI",
      na.value = "transparent",
      guide = ggplot2::guide_legend(title = "IRFI Change", keyheight = grid::unit(1.2, "lines"),
                                    keywidth = grid::unit(1.2, "lines"),
                                    label.theme = ggplot2::element_text(size = 9),
                                    title.theme = ggplot2::element_text(size = 10, face = "bold"))
    ) +
    ggplot2::geom_sf(data = site_proj_original, fill = NA, color = "black",
                     linewidth = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("Knowledge Gained (IRFI Improvement)") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  # ---------------------------------------------------------------------
  # PAGE 2 PLOTS: INITIAL RICHNESS + RICHNESS GAIN
  # ---------------------------------------------------------------------
  
  p2_initial_richness <- ggplot2::ggplot(initial_rich_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical',
                   legend.key.width = grid::unit(0.6, "cm"),
                   plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(9, "Spectral"),
      limits = c(0, richness_max_initial),
      breaks = richness_breaks,
      labels = round(richness_breaks, 1),
      na.value = "transparent",
      guide = ggplot2::guide_legend(title = "N Taxa", keyheight = grid::unit(1.2, "lines"),
                                    keywidth = grid::unit(1.2, "lines"),
                                    label.theme = ggplot2::element_text(size = 9),
                                    title.theme = ggplot2::element_text(size = 10, face = "bold"))
    ) +
    ggplot2::geom_sf(data = site_proj_original, fill = NA, color = "black",
                     linewidth = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("Initial Survey - Species Richness") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  p2_richness_gain <- ggplot2::ggplot(richness_gain_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical',
                   legend.key.width = grid::unit(0.6, "cm"),
                   plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#1A9850",
      limits = c(0, richness_gain_max),
      breaks = richness_gain_breaks,
      labels = round(richness_gain_breaks, 1),
      name = "Taxa Gained",
      na.value = "transparent",
      guide = ggplot2::guide_legend(title = "Taxa Gained", keyheight = grid::unit(1.2, "lines"),
                                    keywidth = grid::unit(1.2, "lines"),
                                    label.theme = ggplot2::element_text(size = 9),
                                    title.theme = ggplot2::element_text(size = 10, face = "bold"))
    ) +
    ggplot2::geom_sf(data = site_proj_original, fill = NA, color = "black",
                     linewidth = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("Richness Gain (Taxa Added)") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  # ---------------------------------------------------------------------
  # PAGE 3 PLOTS: RICHNESS BY IRFI + IMPROVEMENT SCATTERPLOT
  # ---------------------------------------------------------------------
  
  # Categorize cells by IRFI quartile
  cell_data_valid$irfi_category <- cut(
    cell_data_valid$irfi_initial,
    breaks = quantile(cell_data_valid$irfi_initial, probs = c(0, 0.25, 0.5, 0.75, 1)),
    labels = c("Low IRFI\n(Q1)", "Medium-Low\n(Q2)", "Medium-High\n(Q3)", "High IRFI\n(Q4)"),
    include.lowest = TRUE
  )
  
  p3_boxplot <- ggplot2::ggplot(cell_data_valid, ggplot2::aes(x = .data$irfi_category, 
                                                              y = .data$richness_gain)) +
    ggplot2::geom_boxplot(fill = "#1A9850", alpha = 0.6, outlier.color = "red") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 9)) +
    ggplot2::ggtitle("Richness Gain by Initial IRFI Category") +
    ggplot2::xlab("Initial IRFI Quartile") +
    ggplot2::ylab("Taxa Gained") +
    ggplot2::labs(subtitle = "Did high-ignorance cells gain more taxa?")
  
  p3_scatter <- ggplot2::ggplot(cell_data_valid, 
                                ggplot2::aes(x = .data$irfi_initial, 
                                             y = .data$irfi_improvement)) +
    ggplot2::geom_point(alpha = 0.4, color = "darkblue", size = 1.5) +
    ggplot2::geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::ggtitle(paste0("Targeting Efficiency (r = ", round(cor_ignorance_improvement, 2), ")")) +
    ggplot2::xlab("Initial IRFI") +
    ggplot2::ylab("IRFI Improvement") +
    ggplot2::labs(subtitle = "Positive correlation = effective targeting")
  
  # ============================================================================
  # SECTION 11: SAVE PDF OUTPUT
  # ============================================================================
  
  msg("Saving PDF report...")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_report.pdf"))
  tif_improvement_path <- file.path(output_dir, paste0(output_prefix, "_improvement.tif"))
  tif_richness_gain_path <- file.path(output_dir, paste0(output_prefix, "_richness_gain.tif"))
  
  # Save rasters
  terra::writeRaster(improvement_raster, filename = tif_improvement_path, overwrite = TRUE)
  terra::writeRaster(richness_gain_raster, filename = tif_richness_gain_path, overwrite = TRUE)
  msg(paste0("  Improvement raster saved: ", tif_improvement_path))
  msg(paste0("  Richness gain raster saved: ", tif_richness_gain_path))
  
  # Create PDF (4 pages)
  grDevices::pdf(pdf_path, width = 11.69, height = 8.27, onefile = TRUE)
  
  # Page 1
  gridExtra::grid.arrange(
    p1_initial_mrfi, p1_improvement, ncol = 2,
    top = grid::textGrob("Page 1: Initial State vs Knowledge Gained", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  
  # Page 2
  gridExtra::grid.arrange(
    p2_initial_richness, p2_richness_gain, ncol = 2,
    top = grid::textGrob("Page 2: Initial Richness vs Richness Gain", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  
  # Page 3
  gridExtra::grid.arrange(
    p3_boxplot, p3_scatter, ncol = 2,
    top = grid::textGrob("Page 3: Knowledge Gain Analysis", 
                         gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
  
  # Page 4
  grid::grid.draw(
    gridExtra::grid.arrange(
      top = grid::textGrob("Page 4: Summary Statistics", 
                           gp = grid::gpar(fontsize = 14, fontface = "bold")),
      gridExtra::tableGrob(summary_df, rows = NULL)
    )
  )
  
  grDevices::dev.off()
  
  msg(paste0("  PDF report saved: ", pdf_path))
  msg("Done! MRFI comparison complete.")
  
  # ============================================================================
  # SECTION 12: RETURN RESULTS
  # ============================================================================
  
  results <- list(
    mrfi_initial = mrfi_initial,
    mrfi_updated = mrfi_updated,
    improvement_raster = improvement_raster,
    richness_gain_raster = richness_gain_raster,
    summary = summary_df,
    richness_comparison = richness_stats,
    targeting_efficiency = targeting_stats
  )
  
  class(results) <- c("mrfi_comparison", "list")
  return(results)
}