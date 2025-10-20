
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void update_mus_omegas(arma::mat prob_matrix, arma::mat U, double b, arma::rowvec a, double c, arma::mat G, arma::mat& mus, arma::cube& omegas){
  int K = mus.n_rows;
  int D = mus.n_cols;
  int N = U.n_rows;
  
  arma::rowvec sums_prob_mat = arma::ones<arma::rowvec>(N) * prob_matrix;
  
  for (int k = 0; k < K; k++){
  
    arma::rowvec p1_mu = arma::zeros<arma::rowvec>(D);
    arma::mat p2_1_omega(D, D, arma::fill::zeros);
    
    for(int i = 0; i < N; i++){
      arma::rowvec temp = prob_matrix(i,k)*U.row(i);
      p1_mu = p1_mu + temp;
      p2_1_omega = p2_1_omega + (U.row(i).t()*temp);
    }
    
    mus.row(k) = ( (p1_mu + b*a) / (sums_prob_mat(k) + b) );
    
    arma::mat temp = G + p2_1_omega + (b*a.t()*a) - ( (sums_prob_mat(k) + b)*mus.row(k).t()*mus.row(k) );

    if (!temp.is_sympd()) {
      Rcpp::stop("Matrix is not symmetric positive definite.");
    }

    arma::mat p2_omega = arma::inv_sympd(temp);

    omegas.slice(k) =  (sums_prob_mat(k) + c - D) * p2_omega;
    
  }
  
}