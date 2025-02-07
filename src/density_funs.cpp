
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double trunc_poisson_density(double w, double mean){

  return (std::pow(mean, w)/std::tgamma(w+1.0))*(1.0/(std::exp(mean) - 1));
}

// [[Rcpp::export]]
double lognormal_density(double w, double precision, double mean){

  return (std::pow(precision, 0.5)/(w*std::pow(2*arma::datum::pi, 0.5)))*std::exp(-0.5*precision*std::pow(std::log(w) - mean, 2));

}
