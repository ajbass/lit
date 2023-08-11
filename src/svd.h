#ifndef SVD_H
#define SVD_H
#include <RcppArmadillo.h>
using namespace Rcpp;

List quick_svd(arma::mat X);

#endif
