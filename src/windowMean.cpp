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
NumericMatrix window_mean_from_cumsum(NumericMatrix indat, int wsize) {
  int i, j, N, ni=indat.nrow(), nj=indat.ncol(), rad=(int) (wsize-1)/2 ;
  int imax, jmax;
  NumericMatrix result(ni, nj);

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
      result(i, j) /= wsize*wsize;
    }
  }
  return result;
}

// a fast windowing call that can be used with spatialVx
// No longer used
// [[Rcpp::export]]
NumericMatrix windowMean(NumericMatrix indat, NumericVector radius) {
  int rad=(int) radius[0] ;
  if (rad==0) return indat;
  return window_mean_from_cumsum(cumsum2d(indat), 2*rad+1);
}


// [[Rcpp::export]]
DataFrame score_fss(NumericMatrix fc, NumericMatrix ob,
                    NumericVector thresholds, NumericVector window_sizes) {
  int i, j, k, x, y;
  int n_thresholds=thresholds.length(), n_sizes=window_sizes.length();
  int ni=fc.ncol(), nj=fc.nrow();
  double fss1, fss2;
  NumericVector res_fss(n_thresholds * n_sizes);
  NumericVector res_thresh(n_thresholds * n_sizes);
  NumericVector res_size(n_thresholds * n_sizes);
  NumericMatrix fc2(ni,nj), ob2(ni,nj), frac_fc(ni,nj), frac_ob(ni,nj);
  NumericMatrix cum_fc(ni,nj), cum_ob(ni,nj);

  // TODO:
  // if (ob.nrow() != ni || ob.ncol != nj) ERROR
  //
  for (i=0 ; i < n_thresholds ; i++) {
    // turn into a field value 0/1
    for (x=0; x<ni; x++) {
      for (y=0; y<nj; y++) {
        fc2(x,y) = (fc(x,y) >= thresholds(i)) ? 1 : 0 ;
        ob2(x,y) = (ob(x,y) >= thresholds(i)) ? 1 : 0 ;
      }
    }
//    fc2(_,_) = int(fc(_,_) > thresholds(i)) ;
//    ob2(_,_) = int(obs(_,_) > thresholds(i)) ;
    // calculate cumsum matrices
    cum_fc = cumsum2d(fc2);
    cum_ob = cumsum2d(ob2);
    for (j=0 ; j < n_sizes ; j++) {
      k = i*n_sizes + j;
      res_thresh(k) = thresholds(i);
      res_size(k) = window_sizes(j);
      // fraction matrices
      frac_fc = window_mean_from_cumsum(cum_fc, (int) window_sizes[j]);
      frac_ob = window_mean_from_cumsum(cum_ob, (int) window_sizes[j]);
      // FSS
      // FIXME: boundary points? currently: zero padding
      //   to skip incomplete windows : adapt sum below to ~[rad, ni-rad]
      fss1 = 0. ;
      fss2 = 0. ;
      for (x=0; x<ni; x++) {
        for (y=0; y<nj; y++) {
          fss1 += (frac_fc(x,y)-frac_ob(x,y))*(frac_fc(x,y)-frac_ob(x,y)) ;
          fss2 += frac_fc(x,y)*frac_fc(x,y) + frac_ob(x,y)*frac_ob(x,y) ;
        }
//        Rcout << x << y << fss1 << fss2 ;
      }
      res_fss(k) = 1. - fss1/fss2 ;
        //mean( (frac_fc(_,_)-frac_ob(_,_))^2) /
        //mean(frac_fc^2 + frac_ob^2);
      // other "fuzzy" scores: ETS, ...
    }
  }

  return Rcpp::DataFrame::create(Named("threshold")=res_thresh,
                                 Named("scale")=res_size,
                                 Named("fss")=res_fss);
}
  

