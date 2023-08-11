#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

/**
 * Linear regression (R env)
 *
 * @param Xs matrix of covariates
 * @param Ys matrix of responses
 * @return matrix of residuals
 */
// [[Rcpp::export(.quick_lm_cpp)]]
arma::mat quick_lm_cpp(SEXP Xs, SEXP Ys) {
  Rcpp::NumericMatrix Xr(Xs);
  Rcpp::NumericMatrix Yr(Ys);
  int n = Xr.nrow(), k = Xr.ncol(), p = Yr.ncol();

  arma::mat X(Xr.begin(), n, k, false);
  arma::mat Y(Yr.begin(), n, p, false);

  // fit model y ~ X, extract residuals
  arma::mat coef = arma::solve(X.t() * X, X.t() * Y);
  arma::mat res  = Y - X * coef;
  return(res);
}

/**
 * Linear regression (internal)
 *
 * @param X matrix of covariates
 * @param Y matrix of responses
 * @return matrix of residuals
 */
arma::mat quick_lm(arma::mat X, arma::mat Y) {
  // fit model y ~ X, extract residuals
  arma::mat coef = arma::solve(X.t() * X, X.t() * Y);
  arma::mat res  = Y - X * coef;
  return(res);
}

/**
 * Linear regression (internal)
 *
 * @param X matrix of covariates
 * @param Y vector of responses
 * @return matrix of residuals
 */
arma::mat quick_lm_marg(arma::mat X, arma::colvec Y) {
  // fit model y ~ X, extract residuals
  arma::colvec coef = arma::solve(X.t() * X, X.t() * Y);
  arma::colvec res  = Y - X * coef;
  return(res);
}
