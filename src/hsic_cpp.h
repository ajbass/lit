#ifndef HSIC_CPP_H
#define HSIC_CPP_H
#include <RcppArmadillo.h>
using namespace Rcpp;

double hsic_cpp(arma::mat Kx_v, arma::vec Kx_e,
                arma::mat Ky_v, arma::vec Ky_e);

#endif
