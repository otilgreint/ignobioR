#' @title Compare Two Maps of Relative Floristic Ignorance
#'
#' @description
#' Compares initial and updated floristic surveys focusing on targeting efficiency
#' and distribution of knowledge improvement. Uses 4-class efficiency system
#' based on quartile status and sampling presence.
#'
#' @param data_flor_initial Data frame with initial floristic data.
#' @param data_flor_updated Data frame with updated floristic data.
#' @param site An sf object representing the study area.
#' @param tau Numeric. Temporal decay parameter.
#' @param cellsize Numeric. Raster cell size in meters.
#' @param year_initial Numeric. Year of initial survey (default = min year in data).
#' @param year_updated Numeric. Year of updated survey (default = max year in data).
#' @param excl_areas Optional sf object for exclusion areas.
#' @param CRS.new Numeric EPSG code (default = 3035).
#' @param target_quantile Numeric. Quantile threshold for "high ignorance" (default = 0.75).
#' @param output_dir Character. Output directory (default = "./output").
#' @param output_prefix Character. Filename prefix (default = "MRFI_comparison").
#' @param verbose Logical. Print messages (default = TRUE).
#' @param site_buffer Logical. Expand site boundary (default = FALSE).
#' @param buffer_width Numeric. Buffer width in meters if site_buffer = TRUE.
#' @param mask_method Character. Masking method: "touches" or "covers" (default = "touches").
#' @param use_coverage_weighting Logical. Use coverage weighting in ignorance_map (default = TRUE).
#' @param ... Additional parameters passed to ignorance_map().
#'
#' @return List with comparison results including efficiency maps, statistics, and distribution metrics.
#'
#' @importFrom terra rast values compareGeom writeRaster global
#' @importFrom sf st_as_sf st_crs st_transform st_bbox st_make_valid st_intersection
#' @importFrom stats quantile median IQR sd
#' @importFrom ggplot2 ggplot aes geom_tile geom_histogram geom_point scale_fill_gradientn scale_fill_manual scale_size_continuous scale_color_manual coord_equal theme_classic labs ggtitle theme geom_sf element_text guide_legend
#' @importFrom gridExtra grid.arrange tableGrob ttheme_default
#' @importFrom grDevices pdf dev.off
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid textGrob gpar unit
#' @importFrom scales rescale
#' @importFrom ineq Gini
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
  
  # ============================================================================
  # SECTION 1: INPUT VALIDATION
  # ============================================================================
  
  msg("Validating inputs...")
  if (missing(data_flor_initial) || is.null(data_flor_initial)) stop("data_flor_initial is required")
  if (missing(data_flor_updated) || is.null(data_flor_updated)) stop("data_flor_updated is required")
  if (missing(site) || is.null(site)) stop("site is required")
  if (missing(tau) || is.null(tau)) stop("tau is required")
  if (missing(cellsize) || is.null(cellsize)) stop("cellsize is required")
  
  if (is.null(year_initial)) {
    year_initial <- min(data_flor_initial$year, na.rm = TRUE)
    msg(paste0("  year_initial set to: ", year_initial))
  }
  if (is.null(year_updated)) {
    year_updated <- max(data_flor_updated$year, na.rm = TRUE)
    msg(paste0("  year_updated set to: ", year_updated))
  }
  if (year_updated <= year_initial) stop("year_updated must be greater than year_initial")
  
  # ============================================================================
  # SECTION 2: CALCULATE INITIAL MRFI
  # ============================================================================
  
  msg("Calculating baseline MRFI (T1)...")
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
    output_prefix = paste0(output_prefix, "_T1"),
    site_buffer = site_buffer,
    buffer_width = buffer_width,
    mask_method = mask_method,
    use_coverage_weighting = use_coverage_weighting,
    ...
  )
  
  # ============================================================================
  # SECTION 3: IDENTIFY NEW RECORDS
  # ============================================================================
  
  msg("Identifying new records (T1 -> T2)...")
  
  data_flor_initial$record_id <- paste(
    data_flor_initial$Taxon,
    data_flor_initial$Long,
    data_flor_initial$Lat,
    data_flor_initial$year,
    sep = "_"
  )
  
  data_flor_updated$record_id <- paste(
    data_flor_updated$Taxon,
    data_flor_updated$Long,
    data_flor_updated$Lat,
    data_flor_updated$year,
    sep = "_"
  )
  
  new_records_idx <- !data_flor_updated$record_id %in% data_flor_initial$record_id
  data_flor_new <- data_flor_updated[new_records_idx, ]
  
  msg(paste0("  Found ", nrow(data_flor_new), " new occurrences"))
  
  # ============================================================================
  # SECTION 4: CALCULATE UPDATED MRFI
  # ============================================================================
  
  msg("Calculating updated MRFI (T2)...")
  mrfi_updated <- ignorance_map(
    data_flor = data_flor_updated,
    site = site,
    year_study = year_updated,
    excl_areas = excl_areas,
    CRS.new = CRS.new,
    tau = tau,
    cellsize = cellsize,
    verbose = verbose,
    check_overlap = FALSE,
    output_dir = output_dir,
    output_prefix = paste0(output_prefix, "_T2"),
    site_buffer = site_buffer,
    buffer_width = buffer_width,
    mask_method = mask_method,
    use_coverage_weighting = use_coverage_weighting,
    ...
  )
  
  # ============================================================================
  # SECTION 5: EXTRACT AND VALIDATE RASTERS
  # ============================================================================
  
  msg("Extracting rasters and validating compatibility...")
  mrfi_T1 <- mrfi_initial$MRFI
  mrfi_T2 <- mrfi_updated$MRFI
  rich_T1 <- mrfi_initial$RICH
  rich_T2 <- mrfi_updated$RICH
  
  if (!terra::compareGeom(mrfi_T1, mrfi_T2, stopOnError = FALSE)) {
    stop("Rasters have incompatible geometries")
  }
  
  # ============================================================================
  # SECTION 6: CREATE NEW OCCURRENCES RASTER
  # ============================================================================
  
  msg("Creating new occurrences map...")
  
  crs_sf <- sf::st_crs(CRS.new)
  
  if (nrow(data_flor_new) > 0) {
    new_records_sf <- sf::st_as_sf(
      data_flor_new,
      coords = c("Long", "Lat"),
      crs = 4326,
      remove = FALSE
    )
    new_records_sf <- sf::st_transform(new_records_sf, crs_sf)
    
    new_occ_raster <- terra::rasterize(
      terra::vect(new_records_sf),
      mrfi_T1,
      fun = "length"
    )
    new_occ_raster[is.na(new_occ_raster)] <- 0
    
  } else {
    msg("  Warning: No new records found!")
    new_occ_raster <- mrfi_T1
    terra::values(new_occ_raster) <- 0
    new_records_sf <- NULL
  }
  
  # ============================================================================
  # SECTION 7: CREATE 4-CLASS EFFICIENCY MAP
  # ============================================================================
  
  msg("Creating targeting efficiency classification (4-class system)...")
  
  # Extract clean values
  irfi_T1_values <- terra::values(mrfi_T1)
  irfi_T1_values <- irfi_T1_values[!is.na(irfi_T1_values)]
  
  # Calculate Q75 threshold
  q75_T1 <- quantile(irfi_T1_values, target_quantile, na.rm = TRUE)
  msg(paste0("  High-ignorance threshold (Q", target_quantile * 100, "): ", round(q75_T1, 2)))
  
  efficiency_raster <- mrfi_T1
  
  # Extract cellwise values for classification
  irfi_T1_vals <- terra::values(mrfi_T1)
  sampled_vals <- terra::values(new_occ_raster)
  
  efficiency_class <- rep(NA_integer_, length(irfi_T1_vals))
  
  # Class 1: Optimal (Q4 + sampled)
  efficiency_class[!is.na(irfi_T1_vals) &
                     irfi_T1_vals >= q75_T1 &
                     sampled_vals > 0] <- 1
  
  # Class 2: Missed (Q4 + not sampled)
  efficiency_class[!is.na(irfi_T1_vals) &
                     irfi_T1_vals >= q75_T1 &
                     sampled_vals == 0] <- 2
  
  # Class 3: Stable (Q1-Q3 + not sampled)
  efficiency_class[!is.na(irfi_T1_vals) &
                     irfi_T1_vals < q75_T1 &
                     sampled_vals == 0] <- 3
  
  # Class 4: Additional (Q1-Q3 + sampled)
  efficiency_class[!is.na(irfi_T1_vals) &
                     irfi_T1_vals < q75_T1 &
                     sampled_vals > 0] <- 4
  
  terra::values(efficiency_raster) <- efficiency_class
  
  # Count cells in each class
  eff_counts <- table(efficiency_class)
  eff_total <- sum(!is.na(efficiency_class))
  eff_pct <- (eff_counts / eff_total) * 100
  
  msg(paste0("  Optimal (Q4 + sampled): ",
             eff_counts["1"], " (", round(eff_pct["1"], 1), "%)"))
  msg(paste0("  Missed (Q4 + not sampled): ",
             eff_counts["2"], " (", round(eff_pct["2"], 1), "%) PRIORITY"))
  msg(paste0("  Stable (Q1-Q3 + not sampled): ",
             eff_counts["3"], " (", round(eff_pct["3"], 1), "%)"))
  msg(paste0("  Additional (Q1-Q3 + sampled): ",
             eff_counts["4"], " (", round(eff_pct["4"], 1), "%)"))
  
  # ============================================================================
  # SECTION 8: CALCULATE DISTRIBUTION METRICS
  # ============================================================================
  
  msg("Computing distribution metrics...")
  
  irfi_T2_values <- terra::values(mrfi_T2, na.rm = TRUE)
  
  gini_T1 <- ineq::Gini(irfi_T1_values)
  gini_T2 <- ineq::Gini(irfi_T2_values)
  gini_change_pct <- ((gini_T2 - gini_T1) / gini_T1) * 100
  
  cv_T1 <- sd(irfi_T1_values) / mean(irfi_T1_values)
  cv_T2 <- sd(irfi_T2_values) / mean(irfi_T2_values)
  cv_change_pct <- ((cv_T2 - cv_T1) / cv_T1) * 100
  
  iqr_T1 <- IQR(irfi_T1_values)
  iqr_T2 <- IQR(irfi_T2_values)
  iqr_change_pct <- ((iqr_T2 - iqr_T1) / iqr_T1) * 100
  
  q75_T1_dist <- quantile(irfi_T1_values, 0.75)
  q75_T2_dist <- quantile(irfi_T2_values, 0.75)
  
  n_q4_T1 <- sum(irfi_T1_values >= q75_T1_dist)
  n_q4_T2 <- sum(irfi_T2_values >= q75_T2_dist)
  n_q4_change <- n_q4_T2 - n_q4_T1
  
  distribution_metrics <- data.frame(
    Metric = c("Gini coefficient", "Coefficient of variation", 
               "IQR (interquartile range)", "Number of Q4 cells"),
    T1 = c(round(gini_T1, 3), round(cv_T1, 3), 
           round(iqr_T1, 2), n_q4_T1),
    T2 = c(round(gini_T2, 3), round(cv_T2, 3), 
           round(iqr_T2, 2), n_q4_T2),
    Change = c(
      paste0(round(gini_change_pct, 1), "%"),
      paste0(round(cv_change_pct, 1), "%"),
      paste0(round(iqr_change_pct, 1), "%"),
      paste0(ifelse(n_q4_change >= 0, "+", ""), n_q4_change)
    ),
    Interpretation = c(
      ifelse(gini_change_pct < 0, "More even", "Less even"),
      ifelse(cv_change_pct < 0, "Less variable", "More variable"),
      ifelse(iqr_change_pct < 0, "Compressed", "Expanded"),
      ifelse(n_q4_change < 0, "Fewer high-IRFI", "More high-IRFI")
    )
  )
  
  # ============================================================================
  # SECTION 9: CALCULATE KNOWLEDGE METRICS
  # ============================================================================
  
  msg("Computing knowledge metrics...")
  
  mean_irfi_T1 <- mean(irfi_T1_values)
  mean_irfi_T2 <- mean(irfi_T2_values)
  mean_irfi_change_pct <- ((mean_irfi_T2 - mean_irfi_T1) / mean_irfi_T1) * 100
  
  total_taxa_T1 <- length(unique(data_flor_initial$Taxon))
  total_taxa_T2 <- length(unique(data_flor_updated$Taxon))
  new_taxa <- total_taxa_T2 - total_taxa_T1
  
  total_occ_T1 <- nrow(data_flor_initial)
  total_occ_T2 <- nrow(data_flor_updated)
  new_occ <- total_occ_T2 - total_occ_T1
  
  mean_rich_T1 <- terra::global(rich_T1, "mean", na.rm = TRUE)$mean
  mean_rich_T2 <- terra::global(rich_T2, "mean", na.rm = TRUE)$mean
  mean_rich_change_pct <- ((mean_rich_T2 - mean_rich_T1) / mean_rich_T1) * 100
  
  # ============================================================================
  # SECTION 10: PREPARE VISUALIZATION ELEMENTS
  # ============================================================================
  
  msg("Preparing visualizations...")
  
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
  
  add_boundaries <- function(p) {
    p <- p + ggplot2::geom_sf(data = site_plot, fill = NA, color = "black", 
                              linewidth = 1, inherit.aes = FALSE)
    if (!is.null(excl_plot)) {
      p <- p + ggplot2::geom_sf(data = excl_plot, fill = "gray", color = "black",
                                linewidth = 0.5, linetype = "dotted", 
                                inherit.aes = FALSE, alpha = 0.3)
    }
    return(p)
  }
  
  # ============================================================================
  # SECTION 11: PAGE 1 - SEPARATE SCALES FOR EACH MAP
  # ============================================================================
  
  msg("Creating Page 1: Comparative MRFI maps (separate scales)...")
  
  min_T1 <- terra::global(mrfi_T1, "min", na.rm = TRUE)$min
  max_T1 <- terra::global(mrfi_T1, "max", na.rm = TRUE)$max
  breaks_T1 <- seq(min_T1, max_T1, length.out = 9)
  
  min_T2 <- terra::global(mrfi_T2, "min", na.rm = TRUE)$min
  max_T2 <- terra::global(mrfi_T2, "max", na.rm = TRUE)$max
  breaks_T2 <- seq(min_T2, max_T2, length.out = 9)
  
  mrfi_T1_df <- as.data.frame(mrfi_T1, xy = TRUE)
  names(mrfi_T1_df) <- c("x", "y", "value")
  
  mrfi_T2_df <- as.data.frame(mrfi_T2, xy = TRUE)
  names(mrfi_T2_df) <- c("x", "y", "value")
  
  # T1 map with box legend
  p1_T1 <- ggplot2::ggplot(mrfi_T1_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
      values = scales::rescale(breaks_T1),
      breaks = breaks_T1,
      labels = round(breaks_T1, 1),
      limits = c(min_T1, max_T1),
      na.value = "transparent",
      guide = ggplot2::guide_legend(
        title = "IRFI",
        keyheight = grid::unit(1.2, "lines"),
        keywidth = grid::unit(1.2, "lines"),
        nrow = 9
      )
    ) +
    ggplot2::ggtitle(paste0("Initial Knowledge Gaps (T1: ", year_initial, ")")) +
    ggplot2::labs(subtitle = paste0("Range: ", round(min_T1, 1), " - ", round(max_T1, 1))) +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  p1_T1 <- add_boundaries(p1_T1)
  
  # T2 map with box legend
  p1_T2 <- ggplot2::ggplot(mrfi_T2_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = value),
                       alpha = 0.8, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(9, "Spectral")),
      values = scales::rescale(breaks_T2),
      breaks = breaks_T2,
      labels = round(breaks_T2, 1),
      limits = c(min_T2, max_T2),
      na.value = "transparent",
      guide = ggplot2::guide_legend(
        title = "IRFI",
        keyheight = grid::unit(1.2, "lines"),
        keywidth = grid::unit(1.2, "lines"),
        nrow = 9
      )
    ) +
    ggplot2::ggtitle(paste0("Current Knowledge Gaps (T2: ", year_updated, ")")) +
    ggplot2::labs(subtitle = paste0("Range: ", round(min_T2, 1), " - ", round(max_T2, 1))) +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  p1_T2 <- add_boundaries(p1_T2)
  
  # ============================================================================
  # SECTION 12: PAGE 2 LEFT - NEW OCCURRENCES (LINEAR TIME SCALE, 5 BINS)
  # ============================================================================
  
  msg("Creating Page 2 Left: New occurrences map...")
  
  if (!is.null(new_records_sf) && nrow(new_records_sf) > 0) {
    
    # Sort by year so newer records are plotted on top
    new_records_sf <- new_records_sf[order(new_records_sf$year), ]
    
    # Create 5 bins with equal time intervals (linear scale)
    year_range <- range(new_records_sf$year, na.rm = TRUE)
    year_breaks <- seq(year_range[1], year_range[2], length.out = 6)
    
    new_records_sf$year_bin <- cut(
      new_records_sf$year,
      breaks = year_breaks,
      include.lowest = TRUE,
      labels = paste0(
        round(year_breaks[-length(year_breaks)]), 
        "-", 
        round(year_breaks[-1])
      )
    )
    
    n_bins <- 5  # Always 5 bins
    
    # 5-color palette for years (light to dark = old to new)
    year_colors <- c("#FDE725", "#5DC863", "#21908C", "#3B528B", "#440154")
    
    p2_left <- ggplot2::ggplot() +
      ggplot2::coord_equal() +
      ggplot2::theme_classic() +
      ggplot2::geom_sf(data = site_plot, fill = "#F5F5F5", color = "black", linewidth = 1) +
      ggplot2::geom_sf(
        data = new_records_sf,
        ggplot2::aes(size = uncertainty, color = year_bin),
        alpha = 0.6
      ) +
      ggplot2::scale_size_continuous(
        name = "Spatial\nUncertainty",
        range = c(0.5, 4),
        breaks = c(10, 100, 1000, 10000, 50000, 100000),
        labels = c("10 m", "100 m", "1 km", "10 km", "50 km", "100 km"),
        trans = "log10",
        guide = ggplot2::guide_legend(
          keyheight = grid::unit(1.2, "lines"),
          keywidth = grid::unit(1.2, "lines"),
          nrow = 6,
          order = 2,
          override.aes = list(alpha = 1)
        )
      ) +
      ggplot2::scale_color_manual(
        name = "Year",
        values = year_colors,
        drop = FALSE,
        guide = ggplot2::guide_legend(
          keyheight = grid::unit(1.2, "lines"),
          keywidth = grid::unit(1.2, "lines"),
          nrow = 5,
          order = 1,
          override.aes = list(alpha = 1)
        )
      ) +
      ggplot2::ggtitle("New Occurrences (T1 -> T2)") +
      ggplot2::labs(subtitle = paste0(nrow(new_records_sf), " records")) +
      ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") +
      ggplot2::theme(legend.position = "right")
    
    if (!is.null(excl_plot)) {
      p2_left <- p2_left +
        ggplot2::geom_sf(data = excl_plot, fill = "gray", color = "black",
                         linewidth = 0.5, linetype = "dotted", alpha = 0.3)
    }
  } else {
    p2_left <- ggplot2::ggplot() +
      ggplot2::coord_equal() +
      ggplot2::theme_classic() +
      ggplot2::geom_sf(data = site_plot, fill = "#F5F5F5", color = "black", linewidth = 1) +
      ggplot2::ggtitle("New Occurrences (T1 -> T2)") +
      ggplot2::labs(subtitle = "0 records") +
      ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  }
  
  # ============================================================================
  # SECTION 13: PAGE 2 RIGHT - 4-CLASS EFFICIENCY MAP
  # ============================================================================
  
  msg("Creating Page 2 Right: 4-class targeting efficiency map...")
  
  efficiency_df <- as.data.frame(efficiency_raster, xy = TRUE)
  names(efficiency_df) <- c("x", "y", "class")
  
  efficiency_df$class <- as.integer(as.character(efficiency_df$class))
  
  efficiency_df$class <- factor(
    efficiency_df$class,
    levels = c(1, 4, 3, 2),  # Reordered: Optimal, Additional, Stable, Missed
    labels = c("Optimal", "Additional", "Stable", "Missed")
  )
  
  # 4-class color scheme
  efficiency_colors <- c(
    "Optimal" = "#1A5490",      # Dark blue
    "Additional" = "#87CEEB",   # Light blue
    "Stable" = "#FFFFFF",       # White
    "Missed" = "#FF8C00"        # Dark orange
  )
  
  p2_right <- ggplot2::ggplot(efficiency_df) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = class),
                       alpha = 0.9, colour = "black", linewidth = 0.1) +
    ggplot2::scale_fill_manual(
      values = efficiency_colors,
      na.value = "transparent",
      name = "Sampling Efficiency",
      labels = c(
        paste0("Optimal: ", eff_counts["1"], " (", round(eff_pct["1"], 1), "%)"),
        paste0("Additional: ", eff_counts["4"], " (", round(eff_pct["4"], 1), "%)"),
        paste0("Stable: ", eff_counts["3"], " (", round(eff_pct["3"], 1), "%)"),
        paste0("Missed: ", eff_counts["2"], " (", round(eff_pct["2"], 1), "%)")
      ),
      drop = FALSE,
      guide = ggplot2::guide_legend(
        keyheight = grid::unit(1.2, "lines"),
        keywidth = grid::unit(1.2, "lines")
      )
    ) +
    ggplot2::ggtitle("Targeting Efficiency") +
    ggplot2::labs(subtitle = paste0("Q75 threshold: ", round(q75_T1, 1))) +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") +
    ggplot2::theme(legend.position = "right")
  p2_right <- add_boundaries(p2_right)
  
  # ============================================================================
  # SECTION 14: PAGE 3 - SIDE-BY-SIDE HISTOGRAMS
  # ============================================================================
  
  msg("Creating Page 3: Distribution analysis...")
  
  # Create clean data frames
  irfi_T1_clean <- irfi_T1_values[!is.na(irfi_T1_values)]
  irfi_T2_clean <- irfi_T2_values[!is.na(irfi_T2_values)]
  
  if (length(irfi_T1_clean) == 0) {
    stop("No valid IRFI values for T1")
  }
  if (length(irfi_T2_clean) == 0) {
    stop("No valid IRFI values for T2")
  }
  
  hist_data_T1 <- data.frame(
    IRFI = irfi_T1_clean,
    stringsAsFactors = FALSE
  )
  
  hist_data_T2 <- data.frame(
    IRFI = irfi_T2_clean,
    stringsAsFactors = FALSE
  )
  
  msg(paste0("  T1 histogram data: ", nrow(hist_data_T1), " values"))
  msg(paste0("  T2 histogram data: ", nrow(hist_data_T2), " values"))
  
  p3_hist_T1 <- ggplot2::ggplot(
    data = hist_data_T1,
    mapping = ggplot2::aes(x = IRFI)
  ) +
    ggplot2::geom_histogram(
      fill = "#66C2A5",
      alpha = 0.7,
      bins = 30,
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("T1 Distribution (", year_initial, ")"),
      subtitle = paste0("Mean: ", round(mean(irfi_T1_clean), 1), 
                        " | Median: ", round(median(irfi_T1_clean), 1)),
      x = "IRFI",
      y = "Frequency"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"))
  
  p3_hist_T2 <- ggplot2::ggplot(
    data = hist_data_T2,
    mapping = ggplot2::aes(x = IRFI)
  ) +
    ggplot2::geom_histogram(
      fill = "#FC8D62",
      alpha = 0.7,
      bins = 30,
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("T2 Distribution (", year_updated, ")"),
      subtitle = paste0("Mean: ", round(mean(irfi_T2_clean), 1), 
                        " | Median: ", round(median(irfi_T2_clean), 1)),
      x = "IRFI",
      y = "Frequency"
    ) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"))
  
  # Distribution metrics table
  p3_metrics <- gridExtra::tableGrob(
    distribution_metrics,
    rows = NULL,
    theme = gridExtra::ttheme_default(
      base_size = 10,
      core = list(fg_params = list(hjust = 0, x = 0.05)),
      colhead = list(fg_params = list(hjust = 0, x = 0.05))
    )
  )
  
  # ============================================================================
  # SECTION 15: PAGE 4 - SUMMARY STATISTICS
  # ============================================================================
  
  msg("Creating Page 4: Summary statistics...")
  
  end_time <- Sys.time()
  
  summary_df <- data.frame(
    Statistic = c(
      "Analysis completed",
      "Elapsed time (secs)",
      "Study area (cells)",
      "Optimal (Q4 + sampled)",
      "Additional (Q1-Q3 + sampled)",
      "Stable (Q1-Q3 + not sampled)",
      "Missed (Q4 + not sampled)",
      "Number of Q4 cells",
      "Total taxa",
      "Total occurrences",
      "Mean richness per cell",
      "Mean IRFI",
      "Median IRFI",
      "Max IRFI"
    ),
    Value = c(
      as.character(end_time),
      round(as.numeric(difftime(end_time, start_time, units = "secs"))),
      sum(!is.na(terra::values(mrfi_T1))),
      paste0(eff_counts["1"], " (", round(eff_pct["1"], 1), "%)"),
      paste0(eff_counts["4"], " (", round(eff_pct["4"], 1), "%)"),
      paste0(eff_counts["3"], " (", round(eff_pct["3"], 1), "%)"),
      paste0(eff_counts["2"], " (", round(eff_pct["2"], 1), "%) PRIORITY"),
      paste0(n_q4_T1, " -> ", n_q4_T2,
             " (", distribution_metrics$Change[4], ")"),
      paste0(total_taxa_T1, " -> ", total_taxa_T2, " (+", new_taxa, ")"),
      paste0(total_occ_T1, " -> ", total_occ_T2, " (+", new_occ, ")"),
      paste0(round(mean_rich_T1, 1), " -> ", round(mean_rich_T2, 1),
             " (+", round(mean_rich_change_pct, 1), "%)"),
      paste0(round(mean_irfi_T1, 1), " -> ", round(mean_irfi_T2, 1),
             " (", ifelse(mean_irfi_change_pct > 0, "+", ""), round(mean_irfi_change_pct, 1), "%)"),
      paste0(round(median(irfi_T1_clean), 1), " -> ",
             round(median(irfi_T2_clean), 1)),
      paste0(round(max(irfi_T1_clean), 1), " -> ",
             round(max(irfi_T2_clean), 1))
    )
  )
  
  # ============================================================================
  # SECTION 16: SAVE PDF
  # ============================================================================
  
  msg("Saving PDF output...")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_report.pdf"))
  
  terra::writeRaster(efficiency_raster,
                     file.path(output_dir, paste0(output_prefix, "_efficiency.tif")),
                     overwrite = TRUE)
  terra::writeRaster(new_occ_raster,
                     file.path(output_dir, paste0(output_prefix, "_new_occurrences.tif")),
                     overwrite = TRUE)
  
  site_bbox <- sf::st_bbox(site_plot)
  aspect_ratio <- (site_bbox$xmax - site_bbox$xmin) / (site_bbox$ymax - site_bbox$ymin)
  layout_ncol <- if (aspect_ratio > 1.5) 1 else 2
  
  grDevices::pdf(pdf_path, width = 11.69, height = 8.27, onefile = TRUE)
  
  # Page 1
  gridExtra::grid.arrange(
    p1_T1, p1_T2,
    ncol = layout_ncol,
    top = grid::textGrob(
      "Page 1: Comparative Context (Separate Scales)",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  
  # Page 2
  gridExtra::grid.arrange(
    p2_left, p2_right,
    ncol = 2,
    top = grid::textGrob(
      "Page 2: Targeting Efficiency Assessment",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  
  # Page 3
  gridExtra::grid.arrange(
    grobs = list(p3_hist_T1, p3_hist_T2, p3_metrics),
    layout_matrix = rbind(c(1, 2), c(3, 3)),
    heights = c(2, 1),
    top = grid::textGrob(
      "Page 3: Distribution Analysis",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  
  # Page 4
  grid::grid.draw(
    gridExtra::grid.arrange(
      top = grid::textGrob(
        "Page 4: Summary Statistics",
        gp = grid::gpar(fontsize = 14, fontface = "bold")
      ),
      gridExtra::tableGrob(summary_df)
    )
  )
  
  grDevices::dev.off()
  
  msg(paste0("PDF saved: ", pdf_path))
  msg(paste0("Dark blue (Optimal): ", eff_counts["1"], " cells with sampling"))
  msg(paste0("Orange (Missed): ", eff_counts["2"], " cells needing attention"))
  
  # ============================================================================
  # SECTION 17: RETURN RESULTS
  # ============================================================================
  
  results <- list(
    mrfi_initial = mrfi_initial,
    mrfi_updated = mrfi_updated,
    efficiency_map = efficiency_raster,
    new_occurrences_map = new_occ_raster,
    efficiency_summary = data.frame(
      Class = c("Optimal", "Missed", "Stable", "Additional"),
      Count = as.vector(eff_counts[c("1", "2", "3", "4")]),
      Percentage = as.vector(eff_pct[c("1", "2", "3", "4")]),
      Interpretation = c(
        "High-ignorance areas with sampling",
        "High-ignorance areas needing attention",
        "Well-known areas appropriately maintained",
        "Well-known areas with supplemental coverage"
      )
    ),
    distribution_metrics = distribution_metrics,
    summary = summary_df
  )
  
  class(results) <- c("mrfi_comparison", "list")
  
  return(results)
}