#' Load NDVI Example Raster
#'
#' @description Loads the example NDVI raster for the study area.
#' 
#' @details
#' Source: Copernicus Land Monitoring Service (2020-present). 
#' Normalised Difference Vegetation Index. 300m resolution raster product.
#' Available at: https://land.copernicus.eu/
#' 
#' @return A terra SpatRaster object with NDVI values.
#' @export
#' @examples
#' \dontrun{
#' ndvi <- load_ndvi_example()
#' terra::plot(ndvi)
#' }
load_ndvi_example <- function() {
  file_path <- system.file("extdata/ndvi_example.tif", package = "ignobioR")
  
  if (file_path == "" || !file.exists(file_path)) {
    stop("NDVI example file not found. Please reinstall the package.")
  }
  
  terra::rast(file_path)
}