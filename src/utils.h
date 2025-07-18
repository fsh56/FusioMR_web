#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

// Function declaration
Rcpp::NumericVector my_rinvgamma(int n, double shape, double rate);
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
arma::mat my_rinvwishart(double nu, arma::mat S);
Rcpp::NumericVector my_rdirichlet(int n, Rcpp::NumericVector alpha);

#endif
