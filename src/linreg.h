#ifndef LINREG_H
#define LINREG_H
#include <RcppArmadillo.h>

arma::mat quick_lm(arma::mat X, arma::mat Y);
arma::mat quick_lm_marg(arma::mat X, arma::colvec Y);
arma::mat quick_lm_cpp(SEXP Xs, SEXP Ys);

#endif
