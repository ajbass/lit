#include <RcppArmadillo.h>
#include "Rcpp.h"
using namespace Rcpp;

// Define the method signature
extern "C" {
#include "qfc.h"
}

/**
 * Calculates a p-value for the HSIC test statistic using Liu approach
 *
 * @param q The HSIC test statistic
 * @param lambda A vector of eigenvalue products
 * @return p-value
 */
double liu_cpp(double q, NumericVector lambda) {
  // This is simplified for LIT and will not generalize.
  double p, sigmaX, a, muX, delta, l, c1, c2, c3, c4, s1, s2, sigmaQ;
  c1 = sum(lambda); //+ sum(lambda * delta)
  c2 = sum(pow(lambda, 2));// + 2 * sum(lambda ^ 2 * delta)
  c3 = sum(pow(lambda, 3)); //+ 3 * sum(lambda ^ 3 * delta)
  c4 = sum(pow(lambda, 4));// + 4 * sum(lambda ^ 4 * delta)
  s1 = c3 / (pow(c2, 1.5));
  s2 = c4 / pow(c2, 2);
  sigmaQ = sqrt(2 * c2);
  q = (q - c1) / sigmaQ;

  if ((s1 * s1) > s2) {
    a = 1 / (s1 - sqrt((s1 * s1) - s2));
    delta = s1 * (a * a * a) - (a * a);
    l = (a * a) - (2 * delta);
  } else {
    a = 1 / s1;
    delta = 0;
    l = (c2 * c2 * c2) / (c3 * c3);
  }
  muX = l + delta;
  sigmaX = sqrt(2) * a;
  p = R::pnchisq(q * sigmaX + muX, l, delta, false, false);
  return(p);
}

/**
 * The HSIC testing procedure
 *
 * The inputs are the SVD of two kernel matrices and the HSIC test is
 * applied to determine whether they are statistically independent.
 *
 * @param Kx_v Kernel eigenvectors for matrix x
 * @param Kx_e Kernel eigenvalues for matrix x
 * @param Ky_v Kernel eigenvectors for matrix y
 * @param Ky_e Kernel eigenvalues for matrix y
 * @return A p-value
 */
double hsic_cpp(arma::mat Kx_v, arma::vec Kx_e,
                arma::mat Ky_v, arma::vec Ky_e) {
  //initialization
  int n = Ky_v.n_rows;
  int lim = 10000, id = 0;
  double sigma = 0, p, eps = 0.000001;
  int k = Ky_e.n_rows;

  //calculate statistic
  arma::mat gamma = Ky_v.t() * Kx_v;
  gamma = gamma % gamma;
  arma::vec gamma_sq = gamma * Kx_e;
  double stat = sum(gamma_sq.t() * Ky_e) * n;
  arma::mat Z = Ky_e * Kx_e.t(); // eigenvalue products
  NumericMatrix y = wrap(Z);
  y.attr("dim") = R_NilValue;
  y.sort(true);

  // initialize for C function
  double* yp = y.begin();
  k = y.length();
  NumericVector delta (k);
  IntegerVector df (k, 1);
  double trace [7];
  int* df0 = &df[0];
  double* delta2 = &delta[0];

  // This function is from the CompQuadForm package
  qfc(yp, delta2, df0, &k, &sigma, &stat, &lim, &eps, trace, &id, &p);

  // If non-sensical p-value, use Liu's method to approximate
  p = 1 - p;
  if (p > 1 || p <= 0) {
    p = liu_cpp(stat, y);
  }
  return(p);
}
