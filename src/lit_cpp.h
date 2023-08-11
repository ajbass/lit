#ifndef LIT_CPP_H
#define LIT_CPP_H
#include <RcppArmadillo.h>
using namespace Rcpp;

NumericVector lit_cpp(SEXP Xs, SEXP Ys, SEXP Hs);
NumericVector lit_internal(arma::vec Xs, arma::mat Ys, arma::mat Hs);

#endif
