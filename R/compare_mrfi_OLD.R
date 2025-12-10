#' @title Compare Two Maps of Relative Floristic Ignorance (Incremental Approach)
#'
#' @description
#' Compares initial and updated floristic surveys using incremental MRFI calculation.
#' Focuses on absolute change (dIRFI) with clear visualization using discrete quantile scales.
#'
#' @param data_flor_initial Data frame with initial floristic data.
#' @param data_flor_updated Data frame with updated floristic data.
#' @param site An sf object representing the study area.
#' @param tau Numeric. Temporal decay parameter.
#' @param cellsize Numeric. Raster cell size in meters.
#' @param year_initial Numeric. Year of initial survey.
#' @param year_updated Numeric. Year of updated survey.
#' @param excl_areas Optional sf object for exclusion areas.
#' @param CRS.new Numeric EPSG code (default = 3035).
#' @param target_quantile Numeric. Quantile threshold (default = 0.75).
#' @param output_dir Character. Output directory.
#' @param output_prefix Character. Filename prefix.
#' @param verbose Logical. Print messages.
#' @param site_buffer Logical. Expand site boundary.
#' @param buffer_width Numeric. Buffer width in meters.
#' @param mask_method Character. Masking method.
#' @param use_coverage_weighting Logical. Use coverage weighting.
#' @param ... Additional parameters.
#'
#' @return List with comparison results.
#'
#' @importFrom terra rast values compareGeom writeRaster global
#' @importFrom sf st_as_sf st_crs st_transform st_bbox st_make_valid st_intersection
#' @importFrom stats cor quantile median kruskal.test cor.test
#' @importFrom ggplot2 ggplot aes geom_tile geom_boxplot geom_point scale_fill_gradientn scale_fill_gradient2 scale_fill_gradient coord_equal theme_classic labs ggtitle theme geom_sf element_text guide_legend
#' @importFrom gridExtra grid.arrange tableGrob
#' @importFrom grDevices pdf dev.off
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid textGrob gpar unit
#' @importFrom scales rescale
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
  
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  # Input validation
  msg("Validating inputs...")
  if (missing(data_flor_initial) || is.null(data_flor_initial)) stop("data_flor_initial is required")
  if (missing(data_flor_updated) || is.null(data_flor_updated)) stop("data_flor_updated is required")
  if (missing(site) || is.null(site)) stop("site is required")
  if (missing(tau) || is.null(tau)) stop("tau is required")
  if (missing(cellsize) || is.null(cellsize)) stop("cellsize is required")
  
  if (is.null(year_initial)) {
    year_initial <- min(data_flor_initial$year, na.rm = TRUE)
    msg(paste0("year_initial: ", year_initial))
  }
  if (is.null(year_updated)) {
    year_updated <- max(data_flor_updated$year, na.rm = TRUE)
    msg(paste0("year_updated: ", year_updated))
  }
  if (year_updated <= year_initial) stop("year_updated must be greater than year_initial")
  
  # Calculate baseline
  msg("Calculating baseline MRFI...")
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
  
  # Identify new records
  msg("Identifying new records...")
  data_flor_new <- identify_new_records(
    data_initial = data_flor_initial,
    data_updated = data_flor_updated
  )
  msg(paste0("  Found ", nrow(data_flor_new), " new observations"))
  
  # Calculate updated MRFI
  msg("Calculating updated MRFI...")
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
  
  # Extract rasters
  mrfi_initial_raster <- mrfi_initial$MRFI
  mrfi_updated_raster <- mrfi_updated$MRFI
  rich_initial_raster <- mrfi_initial$RICH
  rich_updated_raster <- mrfi_updated$RICH
  
  msg("Checking compatibility...")
  if (!terra::compareGeom(mrfi_initial_raster, mrfi_updated_raster, stopOnError = FALSE)) {
    stop("Rasters have incompatible geometries")
  }
  
  # Calculate changes
  msg("Calculating changes...")
  delta_irfi_raster <- mrfi_updated_raster - mrfi_initial_raster
  richness_gain_raster <- rich_updated_raster - rich_initial_raster
  
  # Extract values
  cell_data <- data.frame(
    irfi_initial = as.vector(terra::values(mrfi_initial_raster)),
    irfi_updated = as.vector(terra::values(mrfi_updated_raster)),
    delta_irfi = as.vector(terra::values(delta_irfi_raster)),
    rich_initial = as.vector(terra::values(rich_initial_raster)),
    rich_updated = as.vector(terra::values(rich_updated_raster)),
    richness_gain = as.vector(terra::values(richness_gain_raster))
  )
  
  cell_data_valid <- cell_data[complete.cases(cell_data), ]
  msg(paste0("Analyzing ", nrow(cell_data_valid), " valid cells"))
  
  # Statistics
  msg("Computing statistics...")
  mean_initial <- mean(cell_data_valid$irfi_initial, na.rm = TRUE)
  mean_updated <- mean(cell_data_valid$irfi_updated, na.rm = TRUE)
  mean_delta <- mean(cell_data_valid$delta_irfi, na.rm = TRUE)
  mean_reduction_pct <- (-mean_delta / mean_initial) * 100
  
  n_improved <- sum(cell_data_valid$delta_irfi < 0, na.rm = TRUE)
  n_worsened <- sum(cell_data_valid$delta_irfi > 0, na.rm = TRUE)
  pct_improved <- (n_improved / nrow(cell_data_valid)) * 100
  pct_worsened <- (n_worsened / nrow(cell_data_valid)) * 100
  
  mean_rich_gain <- mean(cell_data_valid$richness_gain, na.rm = TRUE)
  total_taxa_gain <- max(cell_data_valid$rich_updated) - max(cell_data_valid$rich_initial)
  
  cor_ignorance_delta <- cor(cell_data_valid$irfi_initial, 
                             cell_data_valid$delta_irfi, 
                             use = "complete.obs")
  
  # Targeting efficiency
  ignorance_threshold <- quantile(cell_data_valid$irfi_initial, probs = target_quantile, na.rm = TRUE)
  most_ignorant <- cell_data_valid$irfi_initial >= ignorance_threshold
  delta_in_target <- mean(cell_data_valid$delta_irfi[most_ignorant], na.rm = TRUE)
  delta_overall <- mean(cell_data_valid$delta_irfi, na.rm = TRUE)
  targeting_ratio <- abs(delta_in_target) / abs(delta_overall)
  targeting_score <- min(targeting_ratio * 50, 100)
  
  # Summary table
  end_time <- Sys.time()
  summary_df <- data.frame(
    Statistic = c(
      "Analysis timestamp",
      "Elapsed time (secs)",
      "Valid cells analyzed",
      "Mean IRFI (initial)",
      "Mean IRFI (updated)",
      "Mean IRFI reduction",
      "Mean IRFI reduction (%)",
      "Median delta IRFI",
      "Cells improved (%)",
      "Cells worsened (%)",
      "Mean richness gain",
      "Total taxa discovered",
      "Correlation (IRFI vs delta)",
      "Targeting efficiency score"
    ),
    Value = c(
      as.character(end_time),
      round(as.numeric(difftime(end_time, start_time, units = "secs")), 2),
      nrow(cell_data_valid),
      round(mean_initial, 2),
      round(mean_updated, 2),
      round(abs(mean_delta), 2),
      round(mean_reduction_pct, 1),
      round(median(cell_data_valid$delta_irfi, na.rm = TRUE), 2),
      round(pct_improved, 1),
      round(pct_worsened, 1),
      round(mean_rich_gain, 1),
      round(total_taxa_gain, 0),
      round(cor_ignorance_delta, 3),
      round(targeting_score, 1)
    )
  )
  
  # Prepare visualization
  msg("Preparing visualizations...")
  crs_sf <- sf::st_crs(CRS.new)
  site_plot <- sf::st_as_sf(site)
  if (is.na(sf::st_crs(site_plot))) sf::st_crs(site_plot) <- 4326
  site_plot <- sf::st_transform(sf::st_make_valid(site_plot), crs_sf)
  
  excl_plot <- NULL
  if (!is.null(excl_areas)) {
    excl_temp <- sf::st_as_sf(excl_areas)
    if (is.na(sf::st_crs(excl_temp))) sf::st_crs(excl_temp) <- 4326
    excl_temp <- sf::st_transform(sf::st_make_valid(excl_temp), crs_sf)
    excl_plot <- sf::st_intersection(excl_temp, site_plot)
    if (length(excl_plot) == 0 || all(sf::st_is_empty(excl_plot))) excl_plot <- NULL
  }
  
  # Convert to data frames
  initial_mrfi_df <- as.data.frame(mrfi_initial_raster, xy = TRUE)
  names(initial_mrfi_df) <- c("x", "y", "value")
  
  delta_irfi_df <- as.data.frame(delta_irfi_raster, xy = TRUE)
  names(delta_irfi_df) <- c("x", "y", "value")
  
  initial_rich_df <- as.data.frame(rich_initial_raster, xy = TRUE)
  names(initial_rich_df) <- c("x", "y", "value")
  
  richness_gain_df <- as.data.frame(richness_gain_raster, xy = TRUE)
  names(richness_gain_df) <- c("x", "y", "value")
  
  # Calculate quantile breaks
  global_max_irfi <- terra::global(mrfi_initial_raster, "max", na.rm = TRUE)$max
  irfi_breaks <- seq(0, global_max_irfi, length.out = 9)
  
  delta_vals <- terra::values(delta_irfi_raster, na.rm = TRUE)
  delta_breaks <- stats::quantile(delta_vals, probs = seq(0, 1, length.out = 9), na.rm = TRUE)
  delta_breaks <- unique(delta_breaks)
  
  rich_max <- terra::global(rich_initial_raster, "max", na.rm = TRUE)$max
  rich_breaks <- seq(0, rich_max, length.out = 9)
  
  rich_gain_max <- terra::global(richness_gain_raster, "max", na.rm = TRUE)$max
  rich_gain_breaks <- seq(0, rich_gain_max, length.out = 9)
  
  # Create plots with DISCRETE legends
  p1_initial <- ggplot2::ggplot(initial_mrfi_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
      values = scales::rescale(irfi_breaks),
      breaks = irfi_breaks,
      labels = round(irfi_breaks, 1),
      limits = c(min(irfi_breaks), max(irfi_breaks)),
      na.value = "transparent",
      guide = ggplot2::guide_legend(
        title = "IRFI",
        keyheight = grid::unit(1.2, "lines"),
        keywidth = grid::unit(1.2, "lines")
      )
    ) +
    ggplot2::geom_sf(data = site_plot, fill = NA, color = "black", linewidth = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("Initial Survey - MRFI") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  if (!is.null(excl_plot)) {
    p1_initial <- p1_initial +
      ggplot2::geom_sf(data = excl_plot, fill = "gray", color = "black",
                       linewidth = 0.5, linetype = "dotted", inherit.aes = FALSE, alpha = 0.3)
  }
  
  p1_delta <- ggplot2::ggplot(delta_irfi_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradient2(
      low = "#1A9850",
      mid = "white",
      high = "#D73027",
      midpoint = 0,
      breaks = delta_breaks,
      labels = round(delta_breaks, 1),
      limits = c(min(delta_breaks), max(delta_breaks)),
      na.value = "transparent",
      guide = ggplot2::guide_legend(
        title = "dIRFI",
        reverse = TRUE,
        keyheight = grid::unit(1.2, "lines"),
        keywidth = grid::unit(1.2, "lines")
      )
    ) +
    ggplot2::geom_sf(data = site_plot, fill = NA, color = "black", linewidth = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("IRFI Change (negative = knowledge gained)") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  if (!is.null(excl_plot)) {
    p1_delta <- p1_delta +
      ggplot2::geom_sf(data = excl_plot, fill = "gray", color = "black",
                       linewidth = 0.5, linetype = "dotted", inherit.aes = FALSE, alpha = 0.3)
  }
  
  p2_initial_rich <- ggplot2::ggplot(initial_rich_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = RColorBrewer::brewer.pal(9, "Spectral"),
      breaks = rich_breaks,
      labels = round(rich_breaks, 1),
      limits = c(0, rich_max),
      na.value = "transparent",
      guide = ggplot2::guide_legend(
        title = "N Taxa",
        keyheight = grid::unit(1.2, "lines"),
        keywidth = grid::unit(1.2, "lines")
      )
    ) +
    ggplot2::geom_sf(data = site_plot, fill = NA, color = "black", linewidth = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("Initial Survey - Species Richness") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  if (!is.null(excl_plot)) {
    p2_initial_rich <- p2_initial_rich +
      ggplot2::geom_sf(data = excl_plot, fill = "gray", color = "black",
                       linewidth = 0.5, linetype = "dotted", inherit.aes = FALSE, alpha = 0.3)
  }
  
  p2_rich_gain <- ggplot2::ggplot(richness_gain_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#1A9850",
      breaks = rich_gain_breaks,
      labels = round(rich_gain_breaks, 1),
      limits = c(0, rich_gain_max),
      na.value = "transparent",
      guide = ggplot2::guide_legend(
        title = "Taxa Gained",
        keyheight = grid::unit(1.2, "lines"),
        keywidth = grid::unit(1.2, "lines")
      )
    ) +
    ggplot2::geom_sf(data = site_plot, fill = NA, color = "black", linewidth = 1, inherit.aes = FALSE) +
    ggplot2::ggtitle("Richness Gain (Taxa Added)") +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  if (!is.null(excl_plot)) {
    p2_rich_gain <- p2_rich_gain +
      ggplot2::geom_sf(data = excl_plot, fill = "gray", color = "black",
                       linewidth = 0.5, linetype = "dotted", inherit.aes = FALSE, alpha = 0.3)
  }
  
  # Boxplot and scatterplot
  cell_data_valid$irfi_category <- cut(
    cell_data_valid$irfi_initial,
    breaks = quantile(cell_data_valid$irfi_initial, probs = c(0, 0.25, 0.5, 0.75, 1)),
    labels = c("Low IRFI\n(Q1)", "Medium-Low\n(Q2)", "Medium-High\n(Q3)", "High IRFI\n(Q4)"),
    include.lowest = TRUE
  )
  
  kw_test <- kruskal.test(richness_gain ~ irfi_category, data = cell_data_valid)
  kw_p_text <- if (kw_test$p.value < 0.001) "p < 0.001" else paste0("p = ", round(kw_test$p.value, 3))
  
  p3_boxplot <- ggplot2::ggplot(cell_data_valid, ggplot2::aes(x = irfi_category, y = richness_gain)) +
    ggplot2::geom_boxplot(fill = "#1A9850", alpha = 0.6, outlier.color = "red") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9)) +
    ggplot2::ggtitle("Richness Gain by Initial IRFI Category") +
    ggplot2::xlab("Initial IRFI Quartile") +
    ggplot2::ylab("Taxa Gained") +
    ggplot2::labs(subtitle = paste0("Kruskal-Wallis: ", kw_p_text))
  
  q4_threshold <- quantile(cell_data_valid$irfi_initial, probs = 0.75, na.rm = TRUE)
  
  cell_data_valid$quadrant <- NA
  cell_data_valid$quadrant[cell_data_valid$irfi_initial < q4_threshold & 
                             cell_data_valid$delta_irfi < 0] <- "Q1-Q3\nImproved"
  cell_data_valid$quadrant[cell_data_valid$irfi_initial >= q4_threshold & 
                             cell_data_valid$delta_irfi < 0] <- "Q4\nImproved"
  cell_data_valid$quadrant[cell_data_valid$irfi_initial < q4_threshold & 
                             cell_data_valid$delta_irfi >= 0] <- "Q1-Q3\nWorsened"
  cell_data_valid$quadrant[cell_data_valid$irfi_initial >= q4_threshold & 
                             cell_data_valid$delta_irfi >= 0] <- "Q4\nWorsened"
  
  cell_data_valid$quadrant <- factor(
    cell_data_valid$quadrant,
    levels = c("Q1-Q3\nImproved", "Q4\nImproved", "Q1-Q3\nWorsened", "Q4\nWorsened")
  )
  
  quad_counts <- table(cell_data_valid$quadrant)
  quad_pcts <- round((quad_counts / nrow(cell_data_valid)) * 100, 1)
  quad_labels <- paste0(c("Q1-Q3 Improved", "Q4 Improved ✓", "Q1-Q3 Worsened", "Q4 Worsened ✗"),
                        "\n", quad_counts, " (", quad_pcts, "%)")
  names(quad_labels) <- levels(cell_data_valid$quadrant)
  
  quad_colors <- c("Q1-Q3\nImproved" = "#BEBEBE", "Q4\nImproved" = "#1A9850",
                   "Q1-Q3\nWorsened" = "#FDB462", "Q4\nWorsened" = "#D73027")
  
  cor_test <- cor.test(cell_data_valid$irfi_initial, cell_data_valid$delta_irfi, method = "pearson")
  cor_p_text <- if (cor_test$p.value < 0.001) "p < 0.001" else paste0("p = ", round(cor_test$p.value, 3))
  
  max_irfi <- max(cell_data_valid$irfi_initial, na.rm = TRUE)
  irfi_breaks_scatter <- unique(c(0, pretty(c(0, max_irfi), n = 10)))
  irfi_breaks_scatter <- irfi_breaks_scatter[irfi_breaks_scatter <= max_irfi]
  
  p3_scatter <- ggplot2::ggplot(cell_data_valid, ggplot2::aes(x = irfi_initial, y = delta_irfi, color = quadrant)) +
    ggplot2::geom_vline(xintercept = q4_threshold, linetype = "dashed", color = "black", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
    ggplot2::geom_point(alpha = 0.6, size = 2) +
    ggplot2::scale_color_manual(values = quad_colors, labels = quad_labels, name = "Sampling Efficiency") +
    ggplot2::scale_x_continuous(trans = ggforce::power_trans(2), breaks = irfi_breaks_scatter, labels = irfi_breaks_scatter) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 8.5)) +
    ggplot2::ggtitle("Sampling Efficiency Matrix") +
    ggplot2::xlab("Initial IRFI (pseudo-log)") +
    ggplot2::ylab("dIRFI") +
    ggplot2::labs(subtitle = paste0("r = ", round(cor_ignorance_delta, 2), ", ", cor_p_text))
  
  quadrant_summary <- data.frame(
    Quadrant = c("Q1-Q3 Improved", "Q4 Improved (SUCCESS)", "Q1-Q3 Worsened", "Q4 Worsened (CRITICAL)"),
    Count = as.vector(quad_counts),
    Percentage = as.vector(quad_pcts)
  )
  
  # Save PDF
  msg("Saving PDF...")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_report.pdf"))
  terra::writeRaster(delta_irfi_raster, file.path(output_dir, paste0(output_prefix, "_delta_irfi.tif")), overwrite = TRUE)
  terra::writeRaster(richness_gain_raster, file.path(output_dir, paste0(output_prefix, "_richness_gain.tif")), overwrite = TRUE)
  
  site_bbox <- sf::st_bbox(site_plot)
  aspect_ratio <- (site_bbox$xmax - site_bbox$xmin) / (site_bbox$ymax - site_bbox$ymin)
  layout_ncol <- if (aspect_ratio > 1.5) 1 else 2
  
  grDevices::pdf(pdf_path, width = 11.69, height = 8.27, onefile = TRUE)
  
  gridExtra::grid.arrange(p1_initial, p1_delta, ncol = layout_ncol,
                          top = grid::textGrob("Page 1: Initial State vs IRFI Change", 
                                               gp = grid::gpar(fontsize = 14, fontface = "bold")))
  
  gridExtra::grid.arrange(p2_initial_rich, p2_rich_gain, ncol = layout_ncol,
                          top = grid::textGrob("Page 2: Initial Richness vs Richness Gain", 
                                               gp = grid::gpar(fontsize = 14, fontface = "bold")))
  
  gridExtra::grid.arrange(p3_boxplot, p3_scatter, ncol = 2,
                          top = grid::textGrob("Page 3: Knowledge Gain Analysis", 
                                               gp = grid::gpar(fontsize = 14, fontface = "bold")))
  
  grid::grid.draw(gridExtra::grid.arrange(
    top = grid::textGrob("Page 4: Summary Statistics", gp = grid::gpar(fontsize = 14, fontface = "bold")),
    gridExtra::tableGrob(summary_df, rows = NULL)
  ))
  
  grDevices::dev.off()
  
  msg(paste0("Done! PDF saved: ", pdf_path))
  
  # Return results
  results <- list(
    mrfi_initial = mrfi_initial,
    mrfi_updated = mrfi_updated,
    delta_irfi_raster = delta_irfi_raster,
    richness_gain_raster = richness_gain_raster,
    summary = summary_df,
    quadrant_analysis = quadrant_summary
  )
  
  class(results) <- c("mrfi_comparison", "list")
  return(results)
}