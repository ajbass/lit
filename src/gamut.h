#ifndef GAMUT_CPP_H
#define GAMUT_CPP_H
#include <RcppArmadillo.h>
using namespace Rcpp;

NumericVector gamut_internal(arma::vec Xs, arma::mat Ys, arma::mat Hs);

#endif
