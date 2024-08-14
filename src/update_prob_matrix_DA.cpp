
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void update_prob_matrix_DA(arma::mat& prob_matrix, arma::mat mus, arma::cube omegas, arma::colvec p, arma::mat U, double temp_beta){
  int K = mus.n_rows;
  int N = U.n_rows;
  
  for(int i = 0; i < N; i++){
  
   arma::rowvec p_vec = arma::zeros<arma::rowvec>(K);
   
   for(int k = 0; k < K; k++){
   
     arma::rowvec temp = U.row(i) - mus.row(k);
     arma::rowvec temp2 = -0.5 * temp * omegas.slice(k) * temp.t();
     p_vec(k) = std::pow( p(k) * std::pow(arma::det(omegas.slice(k)), 0.5) * std::exp(temp2(0)), temp_beta );
     
   }
   
   prob_matrix.row(i) = p_vec/sum(p_vec);

   if(prob_matrix.row(i).has_nan()){
     prob_matrix.row(i) = (1.0/K) * arma::ones<arma::rowvec>(K);
   } 
  
  }
  
}
