#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "linreg.h"
#include "svd.h"
#include "hsic_cpp.h"
#include "crossprod.h"
#include "cauchy.h"

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/**
 * Latent Interaction Testing (R env)
 *
 * @param Xs Genotype Matrix
 * @param Ys Trait Matrix
 * @param Hs Matrix of PCs
 * @return A vector of 3 p-values: wLIT, uLIT, and aLIT.
 */
// [[Rcpp::export(.lit_cpp)]]
NumericVector lit_cpp(SEXP Xs, SEXP Ys, SEXP Hs) {
  // convert R objects
  Rcpp::NumericMatrix Xr(Xs);
  Rcpp::NumericMatrix Yr(Ys);
  Rcpp::NumericMatrix Hr(Hs);
  int n = Xr.nrow(), k = Xr.ncol(), m = Yr.ncol(), r = Hr.ncol();
  arma::mat Xrr(Xr.begin(), n, k, false);
  arma::mat Y(Yr.begin(), n, m, false);
  arma::mat H(Hr.begin(), n, r, false);
  NumericVector p (3);

  // adjust genotypes for structure
  arma::mat res_x = quick_lm(H, Xrr);
  arma::mat design(n, 2, arma::fill::ones);
  design.col(1) = res_x;

  // adjust traits for additive effect
  arma::mat res_y = quick_lm(design, Y);

  // svd genotypes
  List svd_x = quick_svd(res_x);
  arma::mat Ux = svd_x["U"];
  arma::mat dx = svd_x["d"];

  // SQ/CP and remove structure
  arma::mat cp = pairwise_prod(res_y);
  arma::mat combined  = arma::join_rows(cp, arma::pow(res_y, 2));
  res_y = quick_lm(H, combined);

  // unequal eigenvalues
  List svd_y = quick_svd(res_y);
  arma::mat Uy = svd_y["U"];
  arma::mat dy = svd_y["d"];
  p[0] = hsic_cpp(Ux, dx, Uy, dy);

  // equal eigenvalues
  dy.ones();
  p[1] = hsic_cpp(Ux, dx, Uy,dy);

  // CCT
  p[2] = cauchyP(p);

  return(p);
}

/**
 * Latent Interaction Testing (Internal)
 *
 * @param Xs Genotype Matrix
 * @param Ys Trait Matrix
 * @param Hs Matrix of PCs
 * @return A vector of 3 p-values: wLIT, uLIT, and aLIT.
 */
NumericVector lit_internal(arma::vec Xs, arma::mat Ys, arma::mat Hs) {
  NumericVector p (3);

  // adjust genotypes for structure
  arma::mat res_x = quick_lm(Hs, Xs);
  arma::mat design(Ys.n_rows, 2, arma::fill::ones);
  design.col(1) = res_x;//X;

  // adjust traits
  arma::mat res_y = quick_lm(design, Ys);

  // svd genotypes
  List svd_x = quick_svd(res_x);
  arma::mat Ux = svd_x["U"];
  arma::mat dx = svd_x["d"];

  // SQ/CP and remove structure
  arma::mat cp = pairwise_prod(res_y);
  arma::mat combined  = arma::join_rows(cp, arma::pow(res_y, 2));
  res_y = quick_lm(Hs, combined);

  // unequal eigenvalues
  List svd_y = quick_svd(res_y);
  arma::mat Uy = svd_y["U"];
  arma::mat dy = svd_y["d"];
  p[0] = hsic_cpp(Ux, dx, Uy, dy);

  // equal eigenvalues
  dy.ones();
  p[1] = hsic_cpp(Ux, dx, Uy,dy);

  // CCT
  p[2] = cauchyP(p);
  return(p);
}
