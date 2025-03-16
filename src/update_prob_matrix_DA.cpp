
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void update_prob_matrix_DA(arma::mat& prob_matrix, arma::mat mus, arma::cube omegas, arma::colvec p, arma::mat U, double temp_beta){
  int K = mus.n_rows;
  int N = U.n_rows;
  int D = U.n_cols;
  
  for(int i = 0; i < N; i++){
  
   arma::rowvec p_vec = arma::zeros<arma::rowvec>(K);
   
   for(int k = 0; k < K; k++){
   
     arma::rowvec temp = U.row(i) - mus.row(k);
     arma::rowvec temp2 = -0.5 * temp * omegas.slice(k) * temp.t();
     double log_mvn = (-0.5*D)*std::log(2.0*arma::datum::pi) + 0.5*std::log(arma::det(omegas.slice(k))) + temp2(0);
     p_vec(k) = temp_beta*std::log(p(k)) + temp_beta*log_mvn;
     
   }
   
   double max_val = arma::max(p_vec);
   double p2 = std::log( arma::sum( arma::exp(p_vec - arma::ones<arma::rowvec>(K)*max_val) ) ); 

   prob_matrix.row(i) = arma::exp( p_vec - arma::ones<arma::rowvec>(K)*(max_val + p2) );

   if(prob_matrix.row(i).has_nan()){
     prob_matrix.row(i) = (1.0/(K*1.0)) * arma::ones<arma::rowvec>(K);
   } 
  
  }
  
}
