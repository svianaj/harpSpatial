#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

typedef std::vector<int> int_vector;
//' @export
// [[Rcpp::export]]
List sal_identify_objects(NumericMatrix indat, double threshold,
    NumericVector maxobj ) {
  // original fortran code by Daan Degrauwe
  // based on original idea of Wernli
  int i, j, ni=indat.nrow(), nj=indat.ncol() ;
  double pmax;
  int nobj = 0, kb;
  int_vector ib(ni*nj), jb(ni*nj),  ibn(ni*nj), jbn(ni*nj);
  int nb, nbn, ii, jj;
  IntegerMatrix objects(ni, nj);

//  Rcout << "SAL" << std::endl ;
//  Rcout <<  "dim " << ni << " " << nj << std::endl ;
  // 1. a matrix with object number
  while (1) {
    // find starting point as maximum value of points that do not belong to an object yet
    nobj++;
    pmax = 0.;
    // TODO: check data ordering. for R-objects, I guess matrix storage is R-style, not C
    // so we take columns (2nd index) as outer loop
    for (j=0 ; j < nj ; j++) {
      for (i=0 ; i < ni ; i++) {
        if (objects(i, j) == 0 && indat(i, j) > pmax) {
          pmax = indat(i, j);
          ib[0] = i;
          jb[0] = j;
        }
      }
    }
    if (pmax < threshold) break;

    objects(ib[0], jb[0]) = nobj;
    nb= 1;
    // determine the object around this peak value
    while (1) {
      nbn = 0;
      // loop over current object boundary { IB(1:NB), JB(1:NB) }
      // storing neighboring points that are not assigned yet and have value > threshold
      // in new object boundary { IBN(1:NBN), JBN(1:NBN) }
      for (kb = 0; kb < nb ; kb++) {
        i = ib[kb];
        j = jb[kb];
        if (i > 0 && objects(i-1, j) == 0 && indat(i-1, j) > threshold ) {
          ibn[nbn] = i-1;
          jbn[nbn] = j;
          objects(i-1, j) = nobj;
          nbn++;
        }
        if (i < ni - 1 && objects(i+1, j) == 0 && indat(i+1, j) > threshold ) {
          ibn[nbn] = i+1;
          jbn[nbn] = j;
          objects(i+1, j) = nobj;
          nbn++;
        }
        if (j > 0 && objects(i, j-1) == 0 && indat(i, j-1) > threshold) {
          ibn[nbn] = i;
          jbn[nbn] = j-1;
          objects(i, j-1) = nobj;
          nbn++;
        }
        if (j < nj - 1 && objects(i, j+1) == 0 && indat(i, j+1) > threshold) {
          ibn[nbn] = i;
          jbn[nbn] = j+1;
          objects(i, j+1) = nobj;
          nbn++;
        }
      }
      // no new points added: object is finished
      if (nbn == 0) break;
      // the new boundary is formed by the added points
      nb = nbn;
      for (kb=0; kb < nb ; kb++) {
        ib[kb] = ibn[kb];
        jb[kb] = jbn[kb];
      }
    }
  }
//  Rcout << "found " << nobj << " objects" << std::endl;
  // 2. find max value, total value, (x, y) centre of gravity for every object
  NumericMatrix sal_stats(nobj, 4);
  //

  for (j=0 ; j < nj ; j++) {
    for (i=0 ; i < ni ; i++) {
      if (objects(i, j) == 0) {
        if (indat(i, j) >= 0.) objects(i, j) = nobj; // last object is "below threshold"
        else objects(i, j) = NA_INTEGER ; // check R code: neg values are from NA
      }
      if (objects(i, j) > 0) {
        kb = objects(i, j) - 1;
        // col 1: total
        sal_stats(kb, 0) += indat(i, j);
        // col 2: max
        if (indat(i, j) > sal_stats(kb, 1))
          sal_stats(kb, 1) = indat(i, j);
        // col 3: c.o.g. X
        sal_stats(kb, 2) += (i+1) * indat(i, j);
        // col 4: c.o.g. Y
        sal_stats(kb, 3) += (j+1) * indat(i, j);
      }
    }
  }

  // 3. the actual (partial) scores (the final combination obs/model is done in R)
  double s0=0., a0=0., lx=0., ly=0., lr=0., rtot;
  // S = 2*(smod-sobs)/(smod+sobs)
  // A <- 2*(dMod-dObs)/(dMod+dObs)
  //    for S, A: treat remainder points as a final object
  //    BUT: for L: NOT, so only nobj-1 objects.
//  v = 0;
//  v1 = 0;
  for (i = 0; i < nobj - 1 ; i++) {
    s0 += sal_stats(i, 0) * sal_stats(i, 0) / sal_stats(i, 1);
    a0 += sal_stats(i, 0) ;
    lx += sal_stats(i, 2) ;
    ly += sal_stats(i, 3) ;
  }
  rtot = a0;
  // for S & A: also include "non-object points" as final object i = nobj-1
  s0 += sal_stats(i, 0) * sal_stats(i, 0) / sal_stats(i, 1);
  a0 += sal_stats(i, 0) ;
  s0 /= a0;

  // L:
  if (rtot == 0) {
    lx = ni/2. ;
    ly = nj/2. ;
    rtot = 1.E-6 ;
  } else {
    lx /= rtot;
    ly /= rtot;
  }
  //
  for (i = 0 ; i < nobj -1 ; i++) {
    lr += sqrt( pow(sal_stats(i, 2) - sal_stats(i, 0) * lx, 2) +
               pow(sal_stats(i, 3) - sal_stats(i, 0) * ly, 2) );
  }
  lr /= rtot;
  // return objects & obj_stat
  return List::create(Rcpp::Named("objects") = objects,
      Rcpp::Named("obj_list") = sal_stats,
      Rcpp::Named("stats") = List::create(
        Rcpp::Named("s") = s0,
        Rcpp::Named("a") = a0,
        Rcpp::Named("lx") = lx,
        Rcpp::Named("ly") = ly,
        Rcpp::Named("lr") = lr)
    );
}
