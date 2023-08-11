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
 * GAMuT
 *
 * @param Xs Genotype Matrix
 * @param Ys Trait Matrix
 * @param Hs Matrix of PCs
 * @return p-value
 */
NumericVector gamut_internal(arma::vec Xs, arma::mat Ys, arma::mat Hs) {
  NumericVector p (3);

  // mean center first, scale genotypes
  arma::mat res_x = quick_lm(Hs, Xs);
  arma::mat design(Ys.n_rows, 1, arma::fill::ones);

  // Adjust traits
  arma::mat res_y = quick_lm(design, Ys);

  // svd genotypes
  List svd_x = quick_svd(res_x);
  arma::mat Ux = svd_x["U"];
  arma::mat dx = svd_x["d"];

  // svd traits
  List svd_y = quick_svd(res_y);
  arma::mat Uy = svd_y["U"];
  arma::mat dy = svd_y["d"];

  p[0] = hsic_cpp(Ux, dx, Uy, dy);

  // equal eigenvalues
  dy.ones();
  p[1] = hsic_cpp(Ux, dx, Uy, dy);

  // CCT
  p[2] = cauchyP(p);
  return(p);
}
