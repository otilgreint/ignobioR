#' @title Virtual Floristic List (Modernized)
#'
#' @description
#' Computes a Virtual Floristic List (VFL): a list of taxa potentially
#' occurring within a study site, each assigned a probability of occurrence
#' based on spatial uncertainty (buffer overlap) and temporal decay.
#'
#' @param data_flor Data frame with 5 columns:
#'   'Taxon', 'Long', 'Lat', 'uncertainty' (radius m), and 'year'.
#' @param site sf polygon (or SpatialPolygonsDataFrame) of the study area.
#' @param year_study Numeric year of the analysis (default = current year).
#' @param excl_areas Optional sf polygon(s) of unsuitable zones.
#' @param CRS.new Numeric EPSG code for projected CRS (default = 3035).
#' @param tau Percent taxa loss per 100 years (0 ≤ τ < 100).
#' @param upperlimit Max records per taxon used in inclusion–exclusion.
#' @param verbose Logical; print progress messages.
#' @param output_dir Directory for saving outputs.
#' @param output_prefix File-name prefix (default = "VFL").
#'
#' @return A list with:
#' \itemize{
#'   \item{\code{VFL}}{Data frame of taxa and probabilities.}
#'   \item{\code{Statistics}}{Run-metadata table.}
#'   \item{\code{Plots}}{List of ggplot objects.}
#' }
#' @import sf terra future.apply ggplot2 gridExtra
#' @export
#'
virtual_list_mod <- function(data_flor, site, year_study = NULL,
                             excl_areas = NULL, CRS.new = 3035,
                             tau, upperlimit = 20, verbose = TRUE,
                             output_dir = getwd(),
                             output_prefix = "VFL") {
  msg <- function(...) if (verbose) message(...)
  start_time <- Sys.time()
  
  ## ---- 1. Checks ----
  if (is.null(year_study))
    year_study <- as.numeric(format(Sys.Date(), "%Y"))
  stopifnot(tau >= 0, tau < 100)
  req <- c("Taxon", "Long", "Lat", "uncertainty", "year")
  if (!all(req %in% names(data_flor)))
    stop("Missing columns: ", paste(setdiff(req, names(data_flor)), collapse = ", "))
  
  ## ---- 2. CRS handling ----
  crs_sf <- sf::st_crs(CRS.new)
  if (inherits(site, "Spatial")) site <- sf::st_as_sf(site)
  if (is.na(sf::st_crs(site))) sf::st_crs(site) <- 4326
  site_proj <- sf::st_transform(sf::st_make_valid(site), crs_sf)
  
  if (!is.null(excl_areas)) {
    if (inherits(excl_areas, "Spatial")) excl_areas <- sf::st_as_sf(excl_areas)
    if (is.na(sf::st_crs(excl_areas))) sf::st_crs(excl_areas) <- 4326
    excl_proj <- sf::st_union(sf::st_transform(sf::st_make_valid(excl_areas), crs_sf))
  } else {
    excl_proj <- NULL
  }
  
  pts_sf <- sf::st_as_sf(data_flor, coords = c("Long", "Lat"), crs = 4326, remove = FALSE)
  pts_proj <- sf::st_transform(pts_sf, crs_sf)
  
  ## ---- 3. Buffering and intersection ----
  msg("Creating buffers...")
  buffers_sf <- sf::st_buffer(pts_proj, dist = pts_proj$uncertainty)
  if (!is.null(excl_proj))
    buffers_sf <- sf::st_difference(buffers_sf, excl_proj)
  
  msg("Computing intersections with study area...")
  inter_sf <- sf::st_intersection(buffers_sf, site_proj)
  inter_sf$area_buffer <- as.numeric(sf::st_area(buffers_sf)[match(inter_sf$Taxon, buffers_sf$Taxon)])
  inter_sf$area_intersection <- as.numeric(sf::st_area(inter_sf))
  inter_sf$p_occurrence_spatial <- inter_sf$area_intersection / inter_sf$area_buffer
  inter_sf$p_occurrence_temporal <- (1 - (tau / 100))^((year_study - inter_sf$year) / 100)
  inter_sf$p_occurrence_spatiotemporal <- inter_sf$p_occurrence_spatial * inter_sf$p_occurrence_temporal
  
  ## ---- 4. Per-taxon aggregation (inclusion–exclusion) ----
  msg("Aggregating probabilities by taxon (parallel)...")
  compute_taxon_prob <- function(v, upperlimit) {
    l <- length(v)
    if (l == 0) return(c(0, 0, 0, 0))
    if (l == 1) return(c(v[1], 1, v[1], v[1]))
    v <- sort(v, decreasing = TRUE)[1:min(l, upperlimit)]
    sum_terms <- sum(v)
    inter_terms <- sapply(2:length(v), function(i) sum(combn(v, i, FUN = prod)))
    signs <- (-1)^(1:length(inter_terms))
    total <- sum(sum_terms, sum(signs * inter_terms))
    c(total, length(v), max(v), min(v))
  }
  
  taxa <- unique(inter_sf$Taxon)
  future::plan(future::multisession)
  results <- future.apply::future_lapply(
    taxa,
    function(tn) compute_taxon_prob(
      inter_sf$p_occurrence_spatiotemporal[inter_sf$Taxon == tn],
      upperlimit
    )
  )
  future::plan(future::sequential)
  
  VFL <- data.frame(
    Taxon = taxa,
    Estimated_Spatiotemporal_probability = sapply(results, `[`, 1),
    Number_of_records = sapply(results, `[`, 2),
    Max_probability = sapply(results, `[`, 3),
    Min_probability = sapply(results, `[`, 4)
  ) |>
    dplyr::filter(Estimated_Spatiotemporal_probability > 0) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x * 100, 1)))
  
  ## ---- 5. Statistics ----
  end_time <- Sys.time()
  stats <- data.frame(
    Statistic = c("Started", "Finished", "Elapsed_s", "CRS", "Tau", "Upperlimit",
                  "Records_used", "Taxa"),
    Value = c(
      as.character(start_time), as.character(end_time),
      round(as.numeric(end_time - start_time, units = "secs")),
      as.character(CRS.new), tau, upperlimit,
      nrow(data_flor), length(taxa)
    )
  )
  
  ## ---- 6. Diagnostic plots ----
  p1 <- ggplot2::ggplot(VFL, ggplot2::aes(x = Estimated_Spatiotemporal_probability)) +
    ggplot2::geom_histogram(fill = "#FF6666", alpha = 0.7, bins = 30) +
    ggplot2::theme_classic() + ggplot2::xlab("Spatio-temporal probability (%)")
  
  p2 <- ggplot2::ggplot(VFL, ggplot2::aes(x = Number_of_records,
                                          y = Estimated_Spatiotemporal_probability)) +
    ggplot2::geom_point(alpha = 0.6, color = "darkblue") +
    ggplot2::theme_classic() + ggplot2::xlab("Records per taxon") +
    ggplot2::ylab("Probability (%)")
  
  p3 <- ggplot2::ggplot(VFL, ggplot2::aes(x = Max_probability, y = Min_probability)) +
    ggplot2::geom_point(alpha = 0.6, color = "darkred") +
    ggplot2::theme_classic() + ggplot2::xlab("Max probability (%)") +
    ggplot2::ylab("Min probability (%)")
  
  ## ---- 7. Save outputs ----
  csv_path <- file.path(output_dir, paste0(output_prefix, "_VFL.csv"))
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_VFL.pdf"))
  utils::write.csv(VFL, csv_path, row.names = FALSE)
  grDevices::pdf(pdf_path, onefile = TRUE)
  print(p1); print(p2); print(p3)
  grid::grid.draw(gridExtra::grid.arrange(top = "Summary statistics",
                                          gridExtra::tableGrob(stats)))
  grDevices::dev.off()
  msg(paste0("Files saved to: ", output_dir))
  
  ## ---- 8. Return ----
  list(VFL = VFL, Statistics = stats, Plots = list(p1, p2, p3))
}
