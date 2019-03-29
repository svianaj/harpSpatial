# Speedup SpatialVx using this plug-in replacement for hoods2dsmooth()
# 'setups and '...' are only there for compatibility. 
# USE: hoods2d(...,smooth.fun=fastSmooth)
# Storing the result of cumsum2d wouldn't work, so we calculate it again every time...
fastSmooth <- function(x, lambda, setup=FALSE, ...) {
  if (setup) return(NA)
  if (!is.matrix(x)) stop("Input must be a matrix")
  if (floor(lambda) != lambda) {
        warning("fastSmooth: attempting to give an illegal value for the neighborhood length.  Flooring lambda.")
    lambda <- floor(lambda)
  }
  if (lambda%%2 == 0) {
    warning("fastSmooth: attempting to give an even neighborhood length, subtracting one from it.")
    lambda <- lambda - 1
  }
  if (lambda < 1) {
    warning("fastSmooth: attempting to give an illegal value for the neighborhood length.  Setting to one, and returning x untouched.")
    lambda <- 1
  }
  if (lambda == 1) return(x)
  windowMean(x, (lambda-1)/2) 

}

