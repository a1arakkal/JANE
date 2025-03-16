
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void update_q_prob(arma::colvec& q_prob, arma::mat prob_matrix_W, Rcpp::String model, double N, double h, double l){

  int M = prob_matrix_W.n_rows;
  double p1 = arma::sum(prob_matrix_W.col(4));
  double c = 1.0;

  if (model != "RSR"){
    c = 0.5;    
  } 

  double p1_1 = p1 + h - 1.0;
  double p2 = p1 + (c*N*(N-1.0)) - M + h + l - 2.0;
  q_prob(0) = p1_1/p2;
  
}
