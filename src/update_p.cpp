
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void update_p(arma::mat prob_matrix, arma::colvec& p, arma::colvec nu){

  int K = prob_matrix.n_cols;
  int N = prob_matrix.n_rows;
  
  arma::rowvec sums_prob_mat = arma::ones<arma::rowvec>(N) * prob_matrix;
  arma::rowvec p_vec = sums_prob_mat + nu.t() - 1.0;
  
  p = (p_vec.t())/(N - K + sum(nu));

  
}
