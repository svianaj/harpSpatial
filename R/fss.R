### Fractions Skill Score
### Fast algorithm based on
### N. Faggian et al. (2015)

### fss calculates for a list of scales and thresholds
### 1. get the 'summed_area_table' for a given threshold
### 2. use this for fast calculation of FSS at different scales

### TO DO: 1. check dimensioning, C-code compatibility ...
##         2. sat.fc is a vector (from C code), not a matrix!
##            so either reformat or change indexing -> what is faster?
##         3. return the 3 components in stead of fss? Then you can accumulate over time.



fss <- function(forecast, obs, 
                windows=rbind(c(3,3),c(5,5)),
                thresholds=c(1,5,15) ) {
  if (any(dim(obs) != dim(forecast))) stop("FC and OBS must have same dimensions.")
  nx <- dim(forecast)[1]
  ny <- dim(forecast)[2]
  nscales <- nrows(windows)
  result <- NULL

  for (thr in thresholds) {
#    sat.fc <- .C("summed_area_table",
#                indat=as.integer(forecast>=thr),
#                nx=as.integer(nx),
#                ny=as.integer(ny))$indat
#    sat.fc <- matrix(sat.fc, nrow=nx)
    sat.fc <- apply(apply(forecast>=thr, 1, cumsum), 1, cumsum)
#    sat.obs <- .C("summed_area_table",
#                indat=as.integer(obs>=thr),
#                nx=as.integer(nx),
#                ny=as.integer(ny))$indat
#    sat.obs <- matrix(sat.obs, nrow=nx)
    sat.obs <- apply(apply(obs>=thr, 1, cumsum), 1, cumsum)

    res <- data.frame(threshold=rep(thr, nscales),
                      nx=rep(nx, nscales), ny=rep(ny, nscales),
                      fss=rep(NA, nscales),
                      fss1=rep(NA, nscales),fss2=rep(NA, nscales))

    for (i in 1:nrows(windows)) {
      # works best with odd numbers (symmetric around center)
      sx2 <- floor(windows[i,1]/2)
      sy2 <- floor(windows[i,2]/2)

      x0 <- 1:nx - sx2
      x0[x0<=0] <- NA
      X0 <- rep(x0, ny)
#      X0 <- rep(pmax(1:nx - sx2, 1), ny)
      X1 <- rep(pmin(1:nx + sx2, nx), ny)
      Y0 <- rep(pmax(1:ny - sy2, 1), each=nx)
      Y1 <- rep(pmin(1:ny + sy2, ny), each=nx)

      frac.fc <- (sat.fc[rbind(X0, Y0)] + sat.fc[rbind(X1, Y1)]
                 -sat.fc[rbind(X0, Y1)] - sat.fc[rbind(X0, Y1)])/(sx*sy)
      frac.obs <- (sat.obs[rbind(X0, Y0)] + sat.obs[rbind(X1, Y1)]
                  -sat.obs[rbind(X0, Y1)] - sat.obs[rbind(X0, Y1)])/(sx*sy)
      res$fss1[i] <- sum( (frac.fc - frac.obs)^2 ) # /(nx*ny)
      res$fss2[i] <- sum(frac.fc^2  + frac.obs^2) #/(nx*ny) 
  
    }
    res$fss <- 1 - res$fss1/res$fss2
    if (is.null(result)) result <- res
    else result <- rbind(result, res)
  }
  result
}




