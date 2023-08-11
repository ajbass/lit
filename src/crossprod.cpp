#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

/**
 * Calculates the pairwise products
 *
 * Takes a n by k matrix and calculates all possible pairwise products between
 * the columns of the matrix.
 *
 * @param x Matrix
 * @return A matrix of dimensions n by (all possible pairwise products)
 */
// [[Rcpp::export(.pairwise_prod)]]
arma::mat pairwise_prod(arma::mat x) {
  const int nr = x.n_rows;
  const int nc = x.n_cols;
  int k = 0;

  arma::mat y(nr, nc * (nc - 1) / 2);

  for (int col = 0; col < nc; ++col) {
    for (int col2 = col + 1; col2 < nc; ++col2) {
      for (int i = 0; i < nr; ++i) {
        y(i, k) = x(i, col) * x(i, col2);
      }
      k++;
    }
  }
  return y;
}
