#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "qfc.h"

using namespace Rcpp;

#define PI 3.14159265358979

/**
 * Cauchy combination test
 *
 * @param P Vector of p-values
 * @return p-value
 */
// [[Rcpp::export(.cauchyP)]]
double cauchyP(NumericVector P) {
  int n = P.size() - 1;
  double pval;
  NumericVector T (n);

  for (int i=0; i < n; ++i) {
    if (P[i] < 1e-16) {
      T[i] = 1 / (P[i] * PI);
    } else {
      T[i] = tanpi(0.5 - P[i]);
    }
  }

  double z = sum(T) / n; // the last spot is 0!
  if(z > 1e+15) {
    pval = (1 / z) / PI;
  } else {
    pval = R::pcauchy(z, 0, 1, false, false);
  }

  return(pval);
}
