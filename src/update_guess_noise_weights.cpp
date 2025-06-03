
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
void update_guess_noise_weights(arma::colvec& guess_noise_weights, arma::mat prob_matrix_W, double precision_noise_weights, double prior_mean_mean_noise, double prior_var_mean_noise, Rcpp::String family){

  int M = prob_matrix_W.n_rows;
 
  
  for(int m = 0; m < M; m++){

    if (family == "poisson"){
    
   
    } else if (family == "exp_lognormal"){

      double o_4 = prior_mean_noise

 
    } else {
      
      double p1 = (precision_noise_weights * arma::sum(prob_matrix_W.col(4) % arma::log(prob_matrix_W.col(2)))) + ((1.0/prior_var_mean_noise) * prior_mean_mean_noise);
      double p2 = (precision_noise_weights * arma::sum(prob_matrix_W.col(4))) + (1.0/prior_var_mean_noise));
      guess_noise_weights(0) = p1/p2;
   
    }

  }

      
}

