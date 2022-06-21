#ifndef _FILTERS_H
#define _FILTERS_H
#include <RcppArmadillo.h>
#include <RcppDist.h>

Rcpp::List bsts_posterior_predict(Rcpp::NumericVector , arma::rowvec , arma::mat , arma::mat , arma::mat , arma::mat , arma::mat , arma::mat , arma::mat , arma::mat , arma::vec , arma::mat);
    
#endif
