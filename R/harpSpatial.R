##' @useDynLib harpSpatial
##' @importFrom Rcpp sourceCpp
.onUnload <- function (libpath) {
  library.dynam.unload("harpSpatial", libpath)
}
