#' Load DEM Example Raster
#'
#' @description Loads the example Digital Elevation Model (DEM) raster for the study area.
#' 
#' @details
#' Source: European Space Agency (2024). Copernicus Global Digital Elevation Model.
#' Distributed by OpenTopography.
#' Available at: https://portal.opentopography.org/
#' 
#' @return A terra SpatRaster object with elevation values in meters.
#' @export
#' @examples
#' \dontrun{
#' dem <- load_dem_example()
#' terra::plot(dem)
#' 
#' # Calculate slope
#' slope <- terra::terrain(dem, "slope", unit = "degrees")
#' terra::plot(slope)
#' }
load_dem_example <- function() {
  file_path <- system.file("extdata/dem_example.tif", package = "ignobioR")
  
  if (file_path == "" || !file.exists(file_path)) {
    stop("DEM example file not found. Please reinstall the package.")
  }
  
  terra::rast(file_path)
}