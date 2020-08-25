# various spatial verification scores

#' Run "fuzzy" spatial verification for 1 case
#'
#' @param obfield Observation grid.
#' @param fcfield Forecast field
#' @param thresholds A vector of thresholds
#' @param window_sizes A vector of (odd!) window sizes
#' @return A tibble with columns for threshold, window_size and various scores.
#' @export
verify_fuzzy <- function(obfield, fcfield, thresholds, window_sizes, scores=list("fss")) {
  # you  might as well calculate a fixed set: they're 'cheap'
  if (is.character(scores)) scores <- list(scores)
  if (any(window_sizes %% 2 != 1)) stop("Window sizes must be odd.")

  # basic preparation
  # TODO: field.type="Preciptation" ???
  # we use our optimised Rcpp fastSmooth code

  nthresh <- length(thresholds)
  nwin <- length(window_sizes)

#  result <- tibble::tibble(
#    threshold = rep(thresholds, each=nwin),
#    scale     = rep(window_sizes, nthresh),
#    ets       = as.vector(vv$fuzzy$ets),
#    fss       = as.vector(vv$fss$values),
#    hk        = as.vector(vv$multi.event$hk)
#  )
  # other scores are probably best added in the same call, so we only calculate fractions once
  # but it may not matter, fact
  score_fss(obfield, fcfield, thresholds, window_sizes)
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
#  VXstats = c("ets", "hk", "f", "bias", "mse")
#  s1 <- SpatialVx::vxstats(obfield, fcfield, which.stats = VXstats)
  dimxy <- prod(dim(obfield))
  mse <- sum((fcfield - obfield)^2)/dimxy
  bias <- sum(fcfield - obfield)/dimxy
# TODO ets, hk, f not yet written out
# TODO: baserate, anomaly correlation (needs climatology...)
  try(sal <- SAL(fcfield, obfield, min.rain = 0.1))

  # put all together in a tibble
  result <- tibble::tibble(
    bias  = bias,
    mse   = mse,
    S     = sal$S,
    A     = sal$A,
    L     = sal$L
  )

  result
}

