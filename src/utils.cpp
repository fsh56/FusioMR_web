#include "utils.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector my_rinvgamma(int n, double shape, double rate) {
  Function ff("rinvgamma", Environment::namespace_env("invgamma"));
  NumericVector res = ff(n, Named("shape") = shape, Named("rate") = rate);
  return res;
}

// [[Rcpp::export]]
arma::mat my_rinvwishart(double nu, arma::mat S) {
  return arma::iwishrnd(S, nu);
}

// Helper function for multivariate normal sampling with robust Cholesky
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  sigma = 0.5 * (sigma + sigma.t());
  arma::mat L;
  bool success = arma::chol(L, sigma, "lower");
  if (!success) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, sigma);
    eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
    L = eigvec * arma::diagmat(arma::sqrt(eigval));
  }
  return arma::repmat(mu, 1, n).t() + Y * L.t();
}

// [[Rcpp::export]]
NumericVector my_rdirichlet(int n, NumericVector alpha) {
  Function ff("rdirichlet", Environment::namespace_env("gtools"));
  NumericVector res = ff(Named("n") = n, Named("alpha") = alpha);
  return res;
}

// Implementation of my_rinvwishart
// [[Rcpp::export]]
arma::mat my_rinvwishart(double nu, arma::mat S) {
  return arma::iwishrnd(S, nu);
}

// Helper function for multivariate normal sampling with robust Cholesky
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  sigma = 0.5 * (sigma + sigma.t());
  arma::mat L;
  bool success = arma::chol(L, sigma, "lower");
  if (!success) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, sigma);
    eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
    L = eigvec * arma::diagmat(arma::sqrt(eigval));
  }
  return arma::repmat(mu, 1, n).t() + Y * L.t();
}

// [[Rcpp::export]]
NumericVector my_rdirichlet(int n, NumericVector alpha) {
  Function ff("rdirichlet");
  NumericVector res = ff(Named("n") = n, _["alpha"] = alpha);
  return res;
}
