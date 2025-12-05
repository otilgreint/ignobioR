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
#'   \item{Page 1: Initial MRFI + IRFI Change (negative = knowledge gained)}
#'   \item{Page 2: Initial Richness + Richness Gain}
#'   \item{Page 3: Richness gain by IRFI category + Quadrant analysis}
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
#'   \item{\code{delta_irfi_raster}}{SpatRaster of IRFI changes (negative = reduction)}
#'   \item{\code{richness_gain_raster}}{SpatRaster of richness gains}
#'   \item{\code{summary}}{Data frame with comparison statistics}
#'   \item{\code{richness_comparison}}{Data frame with richness statistics}
#'   \item{\code{targeting_efficiency}}{Data frame with targeting metrics}
#'   \item{\code{quadrant_analysis}}{Data frame with quadrant breakdown}
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
#' **Sign Convention (dIRFI):**
#' \itemize{
#'   \item{dIRFI = IRFI_updated - IRFI_initial}
#'   \item{Negative values indicate knowledge gain (ignorance reduction)}
#'   \item{Positive values indicate knowledge loss (ignorance increase from temporal decay)}
#' }
#'
#' **PDF Output (4 pages, A4 landscape):**
#' \itemize{
#'   \item{Page 1: Initial MRFI (left) + IRFI Change map (right, negative = gain)}
#'   \item{Page 2: Initial Richness (left) + Richness Gain map (right)}
#'   \item{Page 3: Richness by IRFI boxplot (left) + Quadrant analysis (right)}
#'   \item{Page 4: Summary statistics table}
#' }
#'
#' @importFrom terra rast values compareGeom writeRaster global res ext as.polygons crs rasterize vect
#' @importFrom sf st_as_sf st_crs st_transform st_buffer st_bbox st_coordinates st_geometry st_make_valid st_intersection st_intersects st_difference st_union
#' @importFrom stats cor quantile median kruskal.test cor.test
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom ggplot2 ggplot aes geom_tile geom_histogram geom_vline geom_point geom_smooth geom_boxplot scale_fill_gradientn scale_fill_gradient2 scale_fill_gradient coord_equal theme_classic labs ggtitle element_text theme geom_sf annotate
#' @importFrom gridExtra grid.arrange tableGrob
#' @importFrom grDevices pdf dev.off
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid textGrob gpar unit
#'
#' @export

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
    data_flor_initial = data_flor_initial,
    data_flor_new = data_flor_new,
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
  # SECTION 5: EXTRACT RASTERS AND CALCULATE CHANGES
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
  
  # Calculate changes (NEGATIVE = knowledge gain)
  msg("Calculating IRFI changes...")
  delta_irfi_raster <- mrfi_updated_raster - mrfi_initial_raster  # Negative when ignorance decreases
  richness_gain_raster <- rich_updated_raster - rich_initial_raster
  
  # Extract values
  initial_vals <- terra::values(mrfi_initial_raster, na.rm = FALSE)
  updated_vals <- terra::values(mrfi_updated_raster, na.rm = FALSE)
  delta_irfi_vals <- terra::values(delta_irfi_raster, na.rm = FALSE)
  
  rich_initial_vals <- terra::values(rich_initial_raster, na.rm = FALSE)
  rich_updated_vals <- terra::values(rich_updated_raster, na.rm = FALSE)
  richness_gain_vals <- terra::values(richness_gain_raster, na.rm = FALSE)
  
  # Create cell-level data frame
  cell_data <- data.frame(
    cell_id = seq_along(initial_vals),
    irfi_initial = as.vector(initial_vals),
    irfi_updated = as.vector(updated_vals),
    delta_irfi = as.vector(delta_irfi_vals),  # Negative = improvement
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
  
  # Calculate percentage (as reduction, so flip sign for display)
  cell_data_valid$irfi_reduction_pct <- (-cell_data_valid$delta_irfi / 
                                           cell_data_valid$irfi_initial) * 100
  cell_data_valid$irfi_reduction_pct[is.infinite(cell_data_valid$irfi_reduction_pct)] <- NA
  
  # ============================================================================
  # SECTION 6: SUMMARY STATISTICS
  # ============================================================================
  
  msg("Computing summary statistics...")
  
  # IRFI statistics (note: delta_irfi is negative for improvement)
  mean_initial <- mean(cell_data_valid$irfi_initial, na.rm = TRUE)
  mean_updated <- mean(cell_data_valid$irfi_updated, na.rm = TRUE)
  mean_delta <- mean(cell_data_valid$delta_irfi, na.rm = TRUE)
  mean_reduction_pct <- (-mean_delta / mean_initial) * 100  # Report as positive %
  
  median_initial <- median(cell_data_valid$irfi_initial, na.rm = TRUE)
  median_updated <- median(cell_data_valid$irfi_updated, na.rm = TRUE)
  median_delta <- median(cell_data_valid$delta_irfi, na.rm = TRUE)
  
  max_reduction <- min(cell_data_valid$delta_irfi, na.rm = TRUE)  # Most negative = best
  max_increase <- max(cell_data_valid$delta_irfi, na.rm = TRUE)   # Most positive = worst
  
  n_improved <- sum(cell_data_valid$delta_irfi < 0, na.rm = TRUE)  # Negative = improved
  n_worsened <- sum(cell_data_valid$delta_irfi > 0, na.rm = TRUE)  # Positive = worsened
  n_unchanged <- sum(cell_data_valid$delta_irfi == 0, na.rm = TRUE)
  
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
  
  # Correlation: higher initial IRFI should correlate with MORE NEGATIVE delta (better reduction)
  cor_ignorance_reduction <- cor(
    cell_data_valid$irfi_initial,
    cell_data_valid$delta_irfi,
    use = "complete.obs"
  )
  
  ignorance_threshold <- quantile(
    cell_data_valid$irfi_initial,
    probs = target_quantile,
    na.rm = TRUE
  )
  
  most_ignorant_cells <- cell_data_valid$irfi_initial >= ignorance_threshold
  
  # Mean delta in high-IRFI cells (more negative = better)
  delta_in_target <- mean(
    cell_data_valid$delta_irfi[most_ignorant_cells],
    na.rm = TRUE
  )
  
  delta_overall <- mean(cell_data_valid$delta_irfi, na.rm = TRUE)
  
  # Targeting ratio: ratio of reductions (flip signs for interpretation)
  targeting_ratio <- abs(delta_in_target) / abs(delta_overall)
  targeting_score <- min(targeting_ratio * 50, 100)
  
  targeting_stats <- data.frame(
    Statistic = c(
      "Correlation (initial IRFI vs delta IRFI)",
      paste0("Mean delta in top ", (1-target_quantile)*100, "% ignorant cells"),
      "Targeting efficiency score (0-100)"
    ),
    Value = c(
      round(cor_ignorance_reduction, 3),
      round(delta_in_target, 2),
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
      "Mean IRFI reduction",
      "Mean IRFI reduction (%)",
      "Median delta IRFI",
      "Cells improved (%)",
      "Cells worsened (%)",
      "Mean richness (initial)",
      "Mean richness (updated)",
      "Mean richness gain",
      "Total taxa discovered",
      "Correlation (IRFI vs delta)",
      "Targeting efficiency score"
    ),
    Value = c(
      as.character(end_time),
      round(as.numeric(difftime(end_time, start_time, units = "secs")), 2),
      "Incremental (fixed baseline)",
      nrow(cell_data_valid),
      round(mean_initial, 2),
      round(mean_updated, 2),
      round(abs(mean_delta), 2),  # Report as positive (reduction magnitude)
      round(mean_reduction_pct, 1),
      round(median_delta, 2),     # Keep negative sign
      round(pct_improved, 1),
      round(pct_worsened, 1),
      round(mean_rich_initial, 1),
      round(mean_rich_updated, 1),
      round(mean_richness_gain, 1),
      round(total_taxa_gain, 0),
      round(cor_ignorance_reduction, 3),
      round(targeting_score, 1)
    )
  )
  
  # ============================================================================
  # SECTION 9: PREPARE SITE BOUNDARIES FOR PLOTTING
  # ============================================================================
  
  msg("Preparing visualization components...")
  
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
  
  delta_irfi_df <- as.data.frame(delta_irfi_raster, xy = TRUE)
  names(delta_irfi_df) <- c("x", "y", "value")
  
  initial_rich_df <- as.data.frame(rich_initial_raster, xy = TRUE)
  names(initial_rich_df) <- c("x", "y", "value")
  
  richness_gain_df <- as.data.frame(richness_gain_raster, xy = TRUE)
  names(richness_gain_df) <- c("x", "y", "value")
  
  # Calculate scales
  global_max_irfi <- terra::global(mrfi_initial_raster, "max", na.rm = TRUE)$max
  irfi_breaks <- seq(0, global_max_irfi, length.out = 9)
  
  # Scale for delta IRFI (negative = green, positive = red)
  delta_vals <- terra::values(delta_irfi_raster, na.rm = TRUE)
  delta_min <- min(delta_vals, na.rm = TRUE)  # Most negative (best reduction)
  delta_max <- max(delta_vals, na.rm = TRUE)  # Most positive (worst increase)
  delta_breaks <- seq(delta_min, delta_max, length.out = 9)
  
  richness_max_initial <- terra::global(rich_initial_raster, "max", na.rm = TRUE)$max
  richness_breaks <- seq(0, richness_max_initial, length.out = 9)
  
  richness_gain_max <- terra::global(richness_gain_raster, "max", na.rm = TRUE)$max
  richness_gain_breaks <- seq(0, richness_gain_max, length.out = 9)
  
  # ---------------------------------------------------------------------
  # PAGE 1 PLOTS: INITIAL MRFI + IRFI CHANGE
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
  
  p1_delta <- ggplot2::ggplot(delta_irfi_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right", legend.direction = 'vertical',
                   legend.key.width = grid::unit(0.6, "cm"),
                   plot.title = ggplot2::element_text(size = 11, face = "bold")) +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradient2(
      low = "#1A9850",      # Green for negative (knowledge gain)
      mid = "white",
      high = "#D73027",     # Red for positive (knowledge loss)
      midpoint = 0,
      limits = c(delta_min, delta_max),
      breaks = delta_breaks,
      labels = round(delta_breaks, 1),
      name = "dIRFI",
      na.value = "transparent",
      guide = ggplot2::guide_legend(
        title = "dIRFI",
        reverse = TRUE,        # Negatives (Green) at top
        keyheight = grid::unit(1.2, "lines"),
        keywidth = grid::unit(1.2, "lines"),
        label.theme = ggplot2::element_text(size = 9),
        title.theme = ggplot2::element_text(size = 10, face = "bold")
      )
    ) +
    ggplot2::geom_sf(data = site_proj_original, fill = NA, color = "black",
                     linewidth = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("IRFI Change (negative = knowledge gained)") +
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
  # PAGE 3 PLOTS: RICHNESS BY IRFI + QUADRANT ANALYSIS
  # ---------------------------------------------------------------------
  
  # Categorize cells by IRFI quartile (for boxplot)
  cell_data_valid$irfi_category <- cut(
    cell_data_valid$irfi_initial,
    breaks = quantile(cell_data_valid$irfi_initial, probs = c(0, 0.25, 0.5, 0.75, 1)),
    labels = c("Low IRFI\n(Q1)", "Medium-Low\n(Q2)", "Medium-High\n(Q3)", "High IRFI\n(Q4)"),
    include.lowest = TRUE
  )
  
  # Statistical test for boxplot: Kruskal-Wallis
  kw_test <- kruskal.test(richness_gain ~ irfi_category, data = cell_data_valid)
  kw_p_value <- kw_test$p.value
  
  kw_p_text <- if (kw_p_value < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", format(round(kw_p_value, 3), nsmall = 3))
  }
  
  p3_boxplot <- ggplot2::ggplot(cell_data_valid, ggplot2::aes(x = .data$irfi_category, 
                                                              y = .data$richness_gain)) +
    ggplot2::geom_boxplot(fill = "#1A9850", alpha = 0.6, outlier.color = "red") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9),
                   axis.text.x = ggplot2::element_text(size = 9)) +
    ggplot2::ggtitle("Richness Gain by Initial IRFI Category") +
    ggplot2::xlab("Initial IRFI Quartile") +
    ggplot2::ylab("Taxa Gained") +
    ggplot2::labs(subtitle = paste0("Did high-ignorance cells gain more taxa? (Kruskal-Wallis: ", kw_p_text, ")"))
  
  # Statistical test for quadrant plot: Pearson correlation
  cor_test <- cor.test(cell_data_valid$irfi_initial, 
                       cell_data_valid$delta_irfi,
                       method = "pearson")
  cor_p_value <- cor_test$p.value
  
  cor_p_text <- if (cor_p_value < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", format(round(cor_p_value, 3), nsmall = 3))
  }
  
  # Calculate thresholds (RENAMED FOR CLARITY)
  q4_threshold <- quantile(cell_data_valid$irfi_initial, probs = 0.75, na.rm = TRUE)  # 75th percentile = Q4 boundary
  zero_threshold <- 0
  
  # Log threshold decision
  if (verbose) {
    msg(paste0("Quadrant thresholds:"))
    msg(paste0("  IRFI: 75th percentile (Q4 boundary) = ", round(q4_threshold, 2)))
    msg(paste0("  Delta: Zero (improved vs worsened)"))
  }
  
  # Assign quadrants (FIXED - consistent naming)
  cell_data_valid$quadrant <- NA
  cell_data_valid$quadrant[cell_data_valid$irfi_initial < q4_threshold & 
                             cell_data_valid$delta_irfi < zero_threshold] <- "Q1-Q3\nImproved"
  cell_data_valid$quadrant[cell_data_valid$irfi_initial >= q4_threshold & 
                             cell_data_valid$delta_irfi < zero_threshold] <- "Q4 (High IRFI)\nImproved"
  cell_data_valid$quadrant[cell_data_valid$irfi_initial < q4_threshold & 
                             cell_data_valid$delta_irfi >= zero_threshold] <- "Q1-Q3\nWorsened/Unchanged"
  cell_data_valid$quadrant[cell_data_valid$irfi_initial >= q4_threshold & 
                             cell_data_valid$delta_irfi >= zero_threshold] <- "Q4 (High IRFI)\nWorsened/Unchanged"
  
  # Convert to factor with specific order
  cell_data_valid$quadrant <- factor(
    cell_data_valid$quadrant,
    levels = c(
      "Q1-Q3\nImproved",
      "Q4 (High IRFI)\nImproved",
      "Q1-Q3\nWorsened/Unchanged",
      "Q4 (High IRFI)\nWorsened/Unchanged"
    )
  )
  
  # Count cells in each quadrant
  quadrant_counts <- table(cell_data_valid$quadrant)
  quadrant_pcts <- round((quadrant_counts / nrow(cell_data_valid)) * 100, 1)
  
  # Create labels with counts and percentages
  quadrant_labels <- paste0(
    c("Q1-Q3\nImproved",
      "Q4 (Most Ignorant)\nImproved ✓",
      "Q1-Q3\nWorsened/Unchanged",
      "Q4 (Most Ignorant)\nWorsened/Unchanged ✗"),
    "\n",
    quadrant_counts,
    " cells (", 
    quadrant_pcts, 
    "%)"
  )
  names(quadrant_labels) <- levels(cell_data_valid$quadrant)
  
  # Define colors for quadrants
  quadrant_colors <- c(
    "Q1-Q3\nImproved" = "#BEBEBE",                      # Gray - routine
    "Q4 (High IRFI)\nImproved" = "#1A9850",             # Green - success!
    "Q1-Q3\nWorsened/Unchanged" = "#FDB462",            # Orange - not priority
    "Q4 (High IRFI)\nWorsened/Unchanged" = "#D73027"    # Red - critical!
  )
  
  max_irfi <- max(cell_data_valid$irfi_initial, na.rm = TRUE)
  irfi_breaks <- unique(c(0, pretty(c(0, max_irfi), n = 10)))
  irfi_breaks <- irfi_breaks[irfi_breaks <= max_irfi]
  
  # Create quadrant plot with pseudo-log scale
  p3_scatter <- ggplot2::ggplot(cell_data_valid, 
                                ggplot2::aes(x = .data$irfi_initial, 
                                             y = .data$delta_irfi,
                                             color = .data$quadrant)) +
    # Quadrant dividing lines
    ggplot2::geom_vline(xintercept = q4_threshold, linetype = "dashed", 
                        color = "black", linewidth = 1, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = zero_threshold, linetype = "dashed", 
                        color = "black", linewidth = 1, alpha = 0.8) +
    
    # Points colored by quadrant
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    
    # Color scale
    ggplot2::scale_color_manual(
      values = quadrant_colors,
      labels = quadrant_labels,
      name = "Sampling Efficiency"
    ) +
    
    # PSEUDO-LOG SCALE for x-axis (handles IRFI = 0 without -Inf)
    ggplot2::scale_x_continuous(
      trans = ggforce::power_trans(2),
      breaks = irfi_breaks, # Dynamic breaks are used here
      labels = irfi_breaks  # Labels match the breaks
    )
    
    # Theme and styling (FIXED - single theme() call)
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "gray40"),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8, lineheight = 1.3),
      legend.title = ggplot2::element_text(size = 9, face = "bold"),
      legend.key.size = grid::unit(1.2, "lines"),
      legend.spacing.y = grid::unit(0.3, "cm"),
      legend.key.spacing.y = grid::unit(0.4, "cm"),
      panel.grid.major = ggplot2::element_blank(), 
      panel.grid.minor = ggplot2::element_blank(), 
      panel.background = ggplot2::element_blank()
    ) +
    
    # Labels
    ggplot2::ggtitle("Sampling Efficiency Matrix") +
    ggplot2::xlab("Initial IRFI (pseudo-log scale)") +
    ggplot2::ylab("dIRFI") +
    ggplot2::labs(subtitle = paste0("Q4 focus; r = ", 
                                    round(cor_ignorance_reduction, 2), 
                                    ", ", cor_p_text))
  
  # Store quadrant summary for results
  quadrant_summary <- data.frame(
    Quadrant = c("Q1-Q3, Improved", 
                 "Q4 (High IRFI), Improved (SUCCESS)", 
                 "Q1-Q3, Worsened/Unchanged", 
                 "Q4 (High IRFI), Worsened/Unchanged (CRITICAL)"),
    IRFI_Group = c("Q1-Q3", "Q4", "Q1-Q3", "Q4"),
    Change = c("Improved", "Improved", "Worsened/Unchanged", "Worsened/Unchanged"),
    Count = as.vector(quadrant_counts),
    Percentage = as.vector(quadrant_pcts),
    Interpretation = c(
      "Routine areas showing knowledge improvement - expected pattern",
      "Highest-ignorance quartile successfully improved - targeting effective",
      "Routine areas showing decline - not a conservation priority",
      "Highest-ignorance quartile failed to improve - critical gaps requiring immediate attention"
    ),
    stringsAsFactors = FALSE
  )
  
  # ============================================================================
  # SECTION 11: SAVE PDF OUTPUT
  # ============================================================================
  
  msg("Saving PDF report...")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_report.pdf"))
  tif_delta_path <- file.path(output_dir, paste0(output_prefix, "_delta_irfi.tif"))
  tif_richness_gain_path <- file.path(output_dir, paste0(output_prefix, "_richness_gain.tif"))
  
  # Save rasters
  terra::writeRaster(delta_irfi_raster, filename = tif_delta_path, overwrite = TRUE)
  terra::writeRaster(richness_gain_raster, filename = tif_richness_gain_path, overwrite = TRUE)
  msg(paste0("  Delta IRFI raster saved: ", tif_delta_path))
  msg(paste0("  Richness gain raster saved: ", tif_richness_gain_path))
  
  # Create PDF (4 pages)
  grDevices::pdf(pdf_path, width = 11.69, height = 8.27, onefile = TRUE)
  
  # Page 1
  gridExtra::grid.arrange(
    p1_initial_mrfi, p1_delta, ncol = 2,
    top = grid::textGrob("Page 1: Initial State vs IRFI Change", 
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
    delta_irfi_raster = delta_irfi_raster,
    richness_gain_raster = richness_gain_raster,
    summary = summary_df,
    richness_comparison = richness_stats,
    targeting_efficiency = targeting_stats,
    quadrant_analysis = quadrant_summary
  )
  
  class(results) <- c("mrfi_comparison", "list")
  return(results)
}