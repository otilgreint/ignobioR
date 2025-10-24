#' Standardize point inputs
#' @noRd
#' @keywords internal


source("R/utils_input.R")

pts_sf <- standardize_points(data_flor, verbose = TRUE)
site_sf <- standardize_polygon(site, name = "site", verbose = TRUE)
excl_sf <- standardize_polygon(excl_areas, name = "excl_areas", verbose = TRUE)

# utils_input.R

# Standardize point data (data_flor) to sf
standardize_points <- function(data_flor, coords = c("Long", "Lat"), assumed_crs = 4326, verbose = TRUE) {
  msg <- function(...) if (verbose) message(...)
  msg("Standardizing point input (data_flor)...")
  
  if (inherits(data_flor, "sf")) {
    msg(" - input is already sf")
    return(data_flor)
  }
  
  if (inherits(data_flor, "Spatial")) {
    msg(" - converting from sp to sf")
    sf_obj <- sf::st_as_sf(data_flor)
    if (is.na(sf::st_crs(sf_obj))) sf::st_crs(sf_obj) <- assumed_crs
    return(sf_obj)
  }
  
  if (is.data.frame(data_flor)) {
    if (!all(coords %in% names(data_flor))) {
      stop("data_flor is a data.frame but missing coordinate columns: ", paste(coords, collapse = ", "))
    }
    msg(" - converting from data.frame to sf using coords: ", paste(coords, collapse = ", "))
    sf_obj <- sf::st_as_sf(data_flor, coords = coords, crs = assumed_crs, remove = FALSE)
    return(sf_obj)
  }
  
  stop("Unsupported data_flor type. Provide an sf, Spatial* or data.frame.")
}

# Standardize polygon inputs (site, excl_areas)
standardize_polygon <- function(x, name = "polygon", target_crs = 4326, verbose = TRUE) {
  msg <- function(...) if (verbose) message(...)
  if (is.null(x)) {
    msg(name, ": input is NULL -> returning NULL")
    return(NULL)
  }
  msg("Standardizing polygon input: ", name)
  
  if (inherits(x, "sf")) {
    poly_sf <- x
  } else if (inherits(x, "Spatial")) {
    msg(" - converting sp to sf for ", name)
    poly_sf <- sf::st_as_sf(x)
  } else {
    stop(name, " must be an sf or Spatial object (or NULL).")
  }
  
  if (is.na(sf::st_crs(poly_sf))) {
    msg(" - no CRS found for ", name, "; assuming EPSG:", target_crs)
    sf::st_crs(poly_sf) <- target_crs
  }
  
  # Validate geometries
  poly_valid <- sf::st_make_valid(poly_sf)
  return(poly_valid)
}
