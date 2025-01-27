
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void update_prob_matrix_W_DA(arma::mat& prob_matrix_W, Rcpp::String model, arma::colvec beta, arma::mat U, arma::mat X, double q, double temp_beta){

  int M = prob_matrix_W.n_rows;

  for(int m = 0; m < M; m++){
      
      int i = prob_matrix_W(m, 0) - 1.0;
      int j = prob_matrix_W(m, 1) - 1.0;
      arma::rowvec temp = U.row(i) - U.row(j);
      arma::rowvec cross_prod = temp * temp.t();
      double pij = 0;

      if (model == "NDH"){
         
         double eta_exp = std::exp(beta(0) - cross_prod(0));
         pij = eta_exp/(1.0 + eta_exp);
      
      } else if (model == "RS"){

         arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
         x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
         arma::rowvec x_ij_beta = x_ij*beta; 
         double eta_exp = std::exp(x_ij_beta(0) - cross_prod(0));
         pij = eta_exp/(1.0 + eta_exp);

      } else {
       
         arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
         x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
         arma::rowvec x_ij_beta = x_ij*beta; 
         double eta_exp = std::exp(x_ij_beta(0) - cross_prod(0));
         pij = eta_exp/(1.0 + eta_exp);

      }
     
     double z_hat = std::pow(pij, temp_beta)/(std::pow(pij, temp_beta) + std::pow(q*(1.0-pij), temp_beta));
   
     if(std::isnan(z_hat)){

       z_hat = 0.5;
     
     }

     prob_matrix_W(m, 3) = z_hat;
     prob_matrix_W(m, 4) = 1.0-z_hat;
  
  }

}
