#ifndef MARGINAL_CPP_H
#define MARGINAL_CPP_H
#include <RcppArmadillo.h>
using namespace Rcpp;

double marginal_test_F(const arma::mat y, const arma::mat X, const arma::mat H);
NumericVector marginal_internal(arma::vec Xs, arma::mat Ys, arma::mat Hs);

#endif
