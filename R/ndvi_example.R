#' Example NDVI Raster
#'
#' A sample normalized difference vegetation index (NDVI) raster derived from
#' Landsat 8 OLI imagery (Path 192, Row 030; acquired June 11, 2017).
#'
#' This dataset is provided for testing the \code{sampleboost_mod()} function.
#'
#' @format A \code{SpatRaster} with one layer:
#' \describe{
#'   \item{values}{NDVI values ranging from -1 (non-vegetated) to +1 (dense vegetation).}
#' }
#' @source Processed from Landsat 8 imagery (U.S. Geological Survey)
#' @examples
#' data(ndvi_example)
#' terra::plot(ndvi_example)
"ndvi_example"
