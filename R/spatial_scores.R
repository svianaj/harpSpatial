# various spatial verification scores
# for the moment, we often call SpatialVx
# but sometimes we add some speed (e.g. new smoother for fuzzy scores)

#' Run "fuzzy" spatial verification for 1 case
#'
#' @param obfield Observation grid.
#' @param fcfield Forecast field
#' @param thresholds A vector of thresholds
#' @param window_sizes A vector of (odd!) window sizes
#' @return A tibble with columns for threshold, window_size and various scores.
#' @export
verify_fuzzy <- function(obfield, fcfield, thresholds, window_sizes) {
  # fuzzy is needed for ETS
  # fss is needed for ... FSS :-)
  # multi.event is needed for HK (?)
  # you  might as well calculate a fixed set: they're 'cheap'
  vx_methods <- c("fuzzy", "multi.event", "fss")

  oo <- SpatialVx::make.SpatialVx(obfield, fcfield, thresholds=thresholds)
  # TODO: field.type="Preciptation" ???
  # we use our optimised Rcpp fastSmooth code
  vv <- SpatialVx::hoods2d(oo, which.methods = vx_methods,
                           smooth.fun = fastSmooth,
                           levels = window_sizes)

  nthresh <- length(thresholds)
  nwin <- length(window_sizes)

  result <- tibble::tibble(
                   threshold = rep(thresholds, each=nwin),
                   scale     = rep(window_sizes, nthresh),
                   ets       = as.vector(vv$fuzzy$ets),
                   fss       = as.vector(vv$fss$values),
                   hk        = as.vector(vv$multi.event$hk)
  )
  result
}

#' Run spatial verification for 1 case
#'
#' @param obfield Observation grid.
#' @param fcfield Forecast field
#' return A 1-row tibble of scores
#' @export
verify_basic <- function(obfield, fcfield) {
  ## all spatial scores that do not require a threshold, scale etc.
  ## so for a given case (date, time, leadtime), every score is a single number.
  ## we store MSE, not RMSE, because eventually we may want to sum over a period, too.
  VXstats = c("ets", "hk", "f", "bias", "mse")
  s1 <- SpatialVx::vxstats(obfield, fcfield, which.stats = VXstats)
#  mse <- sum((fcfield - obfield)^2)
#  bias <- sum(fcfield - obfield)
# TODO ets, hk, f not yet written out
# TODO: baserate, anomaly correlation (needs climatology...)
  try(sal <- SAL(fcfield, obfield, min.rain=0.1))

  # put all together in a tibble
  result <- tibble::tibble(bias  = s1$bias,
                           mse   = s1$mse,
                           S     = sal$S,
                           A     = sal$A,
                           L     = sal$L)

  result
}

