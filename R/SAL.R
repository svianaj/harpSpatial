###################################################
### SAL library for precipitation verification  ###
###################################################
# implementation by 
# Daan Degrauwe (RMIB)
# modified for harpSpatial by
# Alex Deckmyn (RMIB)
# March 2019

#' SAL score
#' @param fcfield The forecast field
#' @param obfield The observation field
#' @param threshscale Used for threshold
#' @param threshold A 2-vector containing (precipitation) thresholds for fc and obs.
#' @param sameThreshold If TRUE, threshold is a single number for both fields
#' @param maxobj Maximum number of objects to identify
#' @param min.rain The minimum value to be considered
#' @param add.objs If TRUE, the object matrices are returned as attributes of the result.
#' @return a vector with S, A, L values.
#' @export
"SAL" <- function(fcfield, obfield, threshScale=15.,
                  threshold = pmax(.1, c(max(fcfield, na.rm=TRUE), 
                                         max(obfield, na.rm=TRUE)))/threshScale,
                  sameThreshold=FALSE, maxobj=1000, min.rain=1.0, add.objs=FALSE) {
    # set common threshold, if required
    if (sameThreshold) {
      threshold <- rep(max(threshold),2)
    }
#    message('threshold (max/threshScale): ',threshold)

    # check if sizes are equal
    if ( ! all(dim(fcfield)==dim(obfield)) ) {
      stop('Model and Observation precipitation field should have the same dimensions')
    }
	
    # set NA's in both fields the same
    fcfield[is.na(obfield)] <- NA
    obfield[is.na(fcfield)] <- NA

    if ( max(fcfield, na.rm=TRUE) <= min.rain && max(obfield, na.rm=TRUE) <= min.rain ) {
       ##stop('Stop. No precipitation in model and observation above given threshold')
       message('No precipitation in model and observation above given threshold')
       return(data.frame(S=0, A=0, L=0))
    }
	
    # identify objects
    fc_objects <- sal_identify_objects(fcfield, threshold=threshold[1], maxobj=maxobj)
    ob_objects <- sal_identify_objects(obfield, threshold=threshold[2], maxobj=maxobj)

    # combine ob and fc to get S,A,L :
    S <- 2*(fc_objects$stats$s - ob_objects$stats$s) / (fc_objects$stats$s + ob_objects$stats$s)
    if (is.na(S)) S <- 0

    A <- 2*(fc_objects$stats$a - ob_objects$stats$a) / (fc_objects$stats$a + ob_objects$stats$a)
    if (is.na(A)) A <- 0

    nx <- dim(fcfield)[1]
    ny <- dim(fcfield)[2]
    d <-  sqrt((nx-1)^2 + (ny-1)^2)
    L1 <- sqrt((fc_objects$stats$lx - ob_objects$stats$lx)^2 + (fc_objects$stats$ly - ob_objects$stats$ly)^2)
    L2 <- 2 * abs(fc_objects$stats$lr - ob_objects$stats$lr)
    L <- (L1 + L2) / d

    scores <- data.frame("S"=S, "A"=A, "L"=L)
    if (add.objs) {
      attr(scores, "fc_objects") <- fc_objects
      attr(scores, "ob_objects") <- ob_objects
    }
    return(scores)
}

