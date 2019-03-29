// for fast FSS & other "fuzzy" score calculation
// reference: N. Faggian et al. "Fast calculation of the franctions skill score"
//            MAUSAM, 66 (2015) 457-466
// 1. get the 'summed_area_table' for a given threshold (cumsum2d)
// 2. use this for fast calculation for Fractions tables

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cumsum2d(NumericMatrix indat) {
  int i,j, ni=indat.nrow(), nj=indat.ncol() ;
  NumericMatrix result(ni,nj);

  for (i=0; i< ni ; i++) {
    result(i, 0) = indat(i, 0);
    for (j=1; j< nj; j++) result(i, j) = indat(i, j) + result(i, j-1);
  }
  for (j=0; j< nj ; j++) for (i=1; i< ni; i++) result(i, j) += result(i-1, j);
  return result;
}

// [[Rcpp::export]]
NumericMatrix windowMeanFromCumsum(NumericMatrix indat, NumericVector radius) {
  int i, j, N, ni=indat.nrow(), nj=indat.ncol(), rad=(int) radius[0] ;
  int imax, jmax;
  NumericMatrix result(ni, nj);

  N = (2*rad+1)*(2*rad+1);
// Version that does "zero padding":
  for (i=0 ; i < ni; i++) {
    imax = std::min(i+rad, ni-1) ;
    for(j=0 ; j < nj; j++) {
      jmax = std::min(j+rad, nj-1) ;
      result(i, j) = indat(imax, jmax) ;
      if (i > rad) {
        result(i,j) -= indat(i-rad-1, jmax);
        if (j > rad) 
          result(i,j) += indat(i-rad-1, j-rad-1) - indat(imax, j-rad-1);
      }
      else if (j > rad) result(i,j) -= indat(imax, j-rad-1);
      result(i, j) /= N;
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix windowMean(NumericMatrix indat, NumericVector radius) {
  int rad=(int) radius[0] ;
  if (rad==0) return indat;
  return windowMeanFromCumsum(cumsum2d(indat), radius);
}

