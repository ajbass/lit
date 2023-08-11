#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "linreg.h"
#include "crossprod.h"

using namespace Rcpp;

/**
 * Performs an F-test
 *
 * @param X Squared residual or cross product
 * @param y trait
 * @param H population structure
 * @return p-value
 */
double marginal_test_F(const arma::mat & y, const arma::mat & X, const arma::mat & H) {
  int n = X.n_rows, k1 = X.n_cols, k2 = H.n_cols;
  // regress out additive effects
  arma::colvec null_res = quick_lm_marg(H, y);
  arma::colvec alt_res = quick_lm_marg(arma::join_rows(H, X), y);

  // F-test
  double sig02 = arma::as_scalar(arma::trans(null_res) * null_res);
  double sig12 = arma::as_scalar(arma::trans(alt_res) * alt_res);
  double F = ((sig02 - sig12) / (k1)) / (sig12 / (n - k2 - k1));

  double p = R::pf(F, k1, n - k1 - k2, false, false);
  return(p);
}

/**
 * Marginal testing procedure
 *
 * @param Xs Genotype Matrix
 * @param Ys Trait Matrix
 * @param Hs Matrix of PCs
 * @return A vector of p-values for squared residuals and cross products
 */
// [[Rcpp::export(.marginal_internal)]]
NumericVector marginal_internal(arma::vec Xs, arma::mat Ys, arma::mat Hs) {
  // adjust mean effect of structure on genotypes
  arma::mat res_x = quick_lm(Hs, Xs);
  arma::mat design(Ys.n_rows, 2, arma::fill::ones);
  design.col(1) = res_x;

  // Adjust traits for additive effect
  arma::mat res_y = quick_lm(design, Ys);

  // standardize residuals
  int pp = res_y.n_cols;
  arma::rowvec meany(pp);
  arma::rowvec sigmay(pp);
  //
  meany = mean(res_y, 0);
  sigmay = stddev(res_y, 0);
  //
  for(int j=0;j<pp;j++)
  {
    res_y.col(j) = res_y.col(j) - meany(j);
    res_y.col(j) = res_y.col(j) / sigmay(j);
  }

  // SQ/CP
  arma::mat cp = pairwise_prod(res_y);
  res_y = arma::pow(res_y, 2);

  // marginal test for cross products
  int n = cp.n_cols;
  const int n2 = res_y.n_cols;
  NumericVector p (cp.n_cols + res_y.n_cols);
  for (int i=0; i < n; ++i) {
    p[i] = marginal_test_F(cp.col(i), res_x, Hs);
  }

  // marginal test for squared residuals
  for (int i=0; i < n2; ++i) {
    p[i + n] = marginal_test_F(res_y.col(i), res_x, Hs);
  }

  return(p);
}
