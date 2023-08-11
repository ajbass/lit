#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

/**
 * SVD
 *
 * @param X Matrix
 * @return eigenvectors and eigenvalues
 */
List quick_svd(arma::mat X) {
  int p = X.n_cols;

  List ret;

  arma::rowvec meanx(p);
  arma::rowvec sigmax(p);

  meanx=mean(X, 0);
  sigmax=stddev(X, 0);
  for(int j=0; j<p; j++)
  {
    X.col(j) = X.col(j) - meanx(j);
    X.col(j) = X.col(j) / sigmax(j);
  }
  arma::mat U;
  arma::vec s;
  arma::mat V;
  svd_econ(U, s, V, X, "left") ;
  s = pow(s, 2);
  arma::colvec index (s.n_rows, arma::fill::ones);

  for (int i = s.size() - 1; i>=0; --i) {
    if (s(i) > 1e-10) {
      break;
    } else {
      index(i) = 0;
    }
  }
  arma::uvec index2 = find(index);
  arma::vec s_sub = s.rows(index2);
  arma::mat U_sub = U.cols(index2);
  ret["d"] = s_sub;
  ret["U"] = U_sub ;
  return(ret) ;
}
