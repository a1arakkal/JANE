#include <RcppArmadillo.h>
#include "helper_funs.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Inverse logit function
// [[Rcpp::export]]
 double logit_inv(double eta) {
  if(eta >= 0){
    return 1.0/(1.0 + std::exp(-1.0*eta));
  } else {
    double exp_eta = std::exp(eta);
    return exp_eta/(1.0 + exp_eta);
  }
}

// Truncated Poisson density 
// [[Rcpp::export]]
double trunc_poisson_density(double w, double mean, double log){
  
  double temp = -1.0*arma::datum::inf;

  if (w > 0.0){

    temp = w*std::log(mean) - mean - std::lgamma(w + 1.0) - std::log(-1.0*std::expm1(-1.0*mean));

    if(log <= 0.0){

      temp = std::exp(temp);

    } 
    
  }

  return temp;

}

// Log-normal density 
// [[Rcpp::export]]
double lognormal_density(double w, double precision, double mean, double log){

  double temp = -1.0*arma::datum::inf;

  if (w > 0.0){

    temp = 0.5*std::log(precision) - std::log(w) - 0.5*std::log(2.0*arma::datum::pi) - 0.5*precision*std::pow(std::log(w) - mean, 2);

    if(log <= 0.0){

      temp = std::exp(temp);

    } 
    
  }

  return temp;

}

// solve_only_sympd
// [[Rcpp::export]]
arma::colvec solve_only_sympd(arma::mat A, arma::colvec b){

    if (!A.is_sympd(100.0 * arma::datum::eps)) {
      Rcpp::stop("Matrix is not symmetric positive definite.");
    }

    arma::colvec x = arma::solve(A, b, arma::solve_opts::no_approx + arma::solve_opts::likely_sympd);

    return x;

  }

