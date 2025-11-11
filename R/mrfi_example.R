#' Load MRFI Example Raster
#'
#' @description Loads the example Map of Relative Floristic Ignorance (MRFI) raster.
#' 
#' @details
#' This MRFI was generated using the ignorance_map() function from the ignobioR package,
#' based on floristic occurrence data with spatial and temporal uncertainty.
#' Higher values indicate areas with greater floristic knowledge gaps.
#' 
#' @return A terra SpatRaster object with MRFI values.
#' @export
#' @examples
#' \dontrun{
#' mrfi <- load_mrfi_example()
#' terra::plot(mrfi)
#' }
load_mrfi_example <- function() {
  file_path <- system.file("extdata/mrfi_example.tif", package = "ignobioR")
  
  if (file_path == "" || !file.exists(file_path)) {
    stop("MRFI example file not found. Please reinstall the package.")
  }
  
  terra::rast(file_path)
}