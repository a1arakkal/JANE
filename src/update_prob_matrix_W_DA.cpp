
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double trunc_poisson_density(double w, double mean, double log){
  
  double temp = -1.0*arma::datum::inf;

  if (w > 0.0){

    temp = w*std::log(mean) - mean - std::lgamma(w + 1.0) - std::log(1.0 - std::exp(-1.0*mean));

    if(log <= 0.0){

      temp = std::exp(temp);

    } 
    
  }

  return temp;

}

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

// [[Rcpp::export]]
void update_prob_matrix_W_DA(arma::mat& prob_matrix_W, Rcpp::String model, Rcpp::String family, arma::colvec beta, arma::colvec beta2, double precision_weights, double precision_noise_weights, double guess_noise_weights, arma::mat U, arma::mat X, arma::mat X2, double q, double temp_beta){

  int M = prob_matrix_W.n_rows;

  for(int m = 0; m < M; m++){
      
      int i = prob_matrix_W(m, 0) - 1.0;
      int j = prob_matrix_W(m, 1) - 1.0;
      arma::rowvec temp = U.row(i) - U.row(j);
      arma::rowvec cross_prod = temp * temp.t();
      double pij = 0.0;
      double eta_w = 0.0;
      double z_hat = 0.0;
      
      if (model == "NDH"){
         
         double eta_exp = std::exp(beta(0) - cross_prod(0));
         pij = 1.0/(1.0 + (1.0/eta_exp));
        
         if (family != "bernoulli"){
           
           eta_w =  beta2(0);

         }

      } else if (model == "RS"){

         arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
         x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
         arma::rowvec x_ij_beta = x_ij*beta; 
         double eta_exp = std::exp(x_ij_beta(0) - cross_prod(0));
         pij = 1.0/(1.0 + (1.0/eta_exp));

         if (family != "bernoulli"){
           
           arma::rowvec x2_ij = arma::ones<arma::rowvec>(1+X2.n_cols);
           x2_ij(arma::span(1, X2.n_cols)) = X2.row(i) + X2.row(j);
           arma::rowvec x2_ij_beta = x2_ij*beta2; 
           eta_w = x2_ij_beta(0);

         }


      } else {
       
         arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
         x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
         arma::rowvec x_ij_beta = x_ij*beta; 
         double eta_exp = std::exp(x_ij_beta(0) - cross_prod(0));
         pij = 1.0/(1.0 + (1.0/eta_exp));

         if (family != "bernoulli"){
           
           arma::rowvec x2_ij = arma::ones<arma::rowvec>(1+X2.n_cols);
           x2_ij(arma::span(1, X2.n_cols)) = arma::join_rows(X2.row(i).subvec(0, (X2.n_cols*0.5) - 1), X2.row(j).subvec(X2.n_cols*0.5, X2.n_cols - 1));
           arma::rowvec x2_ij_beta = x2_ij*beta2; 
           eta_w = x2_ij_beta(0);

         }


      }
     
     if (family == "bernoulli"){

        double log_density = 0.0;
        double log_density_noise = 0.0;
        double temp1 = temp_beta*(std::log(pij) + log_density);
        double temp2 = temp_beta*(std::log(q*(1.0-pij)) + log_density_noise);
        double max_val = std::max(temp1, temp2);

        z_hat = std::exp( temp1 - max_val - std::log( std::exp(temp1 - max_val) + std::exp(temp2 - max_val) ) );
     
     } else if (family == "poisson"){
        
        double w = prob_matrix_W(m, 2);
        double log_density = trunc_poisson_density(w, std::exp(eta_w), 1.0);
        double log_density_noise = trunc_poisson_density(w, guess_noise_weights, 1.0);
        double temp1 = temp_beta*(std::log(pij) + log_density);
        double temp2 = temp_beta*(std::log(q*(1.0-pij)) + log_density_noise);
        double max_val = std::max(temp1, temp2);

        z_hat = std::exp( temp1 - max_val - std::log( std::exp(temp1 - max_val) + std::exp(temp2 - max_val) ) );
     
     } else {
        
        double w = prob_matrix_W(m, 2);
        double log_density = lognormal_density(w, precision_weights, eta_w, 1.0);
        double log_density_noise = lognormal_density(w, precision_noise_weights, guess_noise_weights, 1.0);
        double temp1 = temp_beta*(std::log(pij) + log_density);
        double temp2 = temp_beta*(std::log(q*(1.0-pij)) + log_density_noise);
        double max_val = std::max(temp1, temp2);

        z_hat = std::exp( temp1 - max_val - std::log( std::exp(temp1 - max_val) + std::exp(temp2 - max_val) ) );
     
     }

     if(std::isnan(z_hat)){

       z_hat = 0.5;
     
     }

     prob_matrix_W(m, 3) = z_hat;
     prob_matrix_W(m, 4) = 1.0-z_hat;
  
  }

}
