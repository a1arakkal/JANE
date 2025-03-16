
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
void update_beta2(arma::colvec& beta2, arma::mat prob_matrix_W, arma::mat f_2, arma::colvec e_2, arma::mat X2, Rcpp::String model, Rcpp::String family){

  int M = prob_matrix_W.n_rows;
 
  arma::colvec p_1 = arma::zeros<arma::colvec>(1+X2.n_cols);
  arma::mat p_2(1+X2.n_cols, 1+X2.n_cols, arma::fill::zeros);

  double p_1_NDH = 0.0;
  double p_2_NDH = 0.0;
  
  for(int m = 0; m < M; m++){

    if (family == "poisson"){
      
      if (model == "NDH"){
          
        double lambda = std::exp(beta2(0));
        double g_1 =  lambda / (1.0 - std::exp(-1.0*lambda));
        double g_2 = ( lambda * (1.0 - (std::exp(-1.0*lambda) * (1.0 + lambda))) )/std::pow(1.0 - std::exp(-1.0*lambda), 2);

        p_1_NDH = p_1_NDH + ( prob_matrix_W(m, 3)*(prob_matrix_W(m, 2) + g_2*beta2(0) - g_1) );
        p_2_NDH = p_2_NDH + ( prob_matrix_W(m, 3)*g_2 );  
          
      } else {
      
        int i = prob_matrix_W(m, 0) - 1.0;
        int j = prob_matrix_W(m, 1) - 1.0;
        
        arma::rowvec x2_ij = arma::ones<arma::rowvec>(1+X2.n_cols);
        arma::rowvec x2_ij_beta = arma::ones<arma::rowvec>(1);

        if (model == "RS"){

          x2_ij(arma::span(1, X2.n_cols)) = X2.row(i) + X2.row(j);
          x2_ij_beta = x2_ij*beta2; 

        } else {

          x2_ij(arma::span(1, X2.n_cols)) = arma::join_rows(X2.row(i).subvec(0, (X2.n_cols*0.5) - 1), X2.row(j).subvec(X2.n_cols*0.5, X2.n_cols - 1));
          x2_ij_beta = x2_ij*beta2; 

        }
        
        double lambda = std::exp(x2_ij_beta(0));
        arma::colvec g_1 = (lambda / (1.0 - std::exp(-1.0*lambda))) * x2_ij.t();
        arma::mat g_2 = (( lambda * (1.0 - (std::exp(-1.0*lambda) * (1.0 + lambda))) )/std::pow(1.0 - std::exp(-1.0*lambda), 2)) * x2_ij.t()*x2_ij;
        
        p_1 = p_1 + ( prob_matrix_W(m, 3)*(prob_matrix_W(m, 2)*x2_ij.t() + g_2*beta2 - g_1) );
        p_2 = p_2 + ( prob_matrix_W(m, 3)*g_2 );  
      
     }
   
   } else {
      
      if (model == "NDH"){

        p_1_NDH = p_1_NDH + ( prob_matrix_W(m, 3) * std::log(prob_matrix_W(m, 2)) );
        p_2_NDH = p_2_NDH + ( prob_matrix_W(m, 3) );  
          
      } else {
      
        int i = prob_matrix_W(m, 0) - 1.0;
        int j = prob_matrix_W(m, 1) - 1.0;
        
        arma::rowvec x2_ij = arma::ones<arma::rowvec>(1+X2.n_cols);

        if (model == "RS"){

          x2_ij(arma::span(1, X2.n_cols)) = X2.row(i) + X2.row(j);

        } else {

          x2_ij(arma::span(1, X2.n_cols)) = arma::join_rows(X2.row(i).subvec(0, (X2.n_cols*0.5) - 1), X2.row(j).subvec(X2.n_cols*0.5, X2.n_cols - 1));

        }
                
        p_1 = p_1 + ( prob_matrix_W(m, 3)*std::log(prob_matrix_W(m, 2))*x2_ij.t() );
        p_2 = p_2 + ( prob_matrix_W(m, 3)*x2_ij.t()*x2_ij );  
      
     
     }
   
   }

  }

  if (model != "NDH") {

    beta2 = arma::inv_sympd(f_2 + p_2) *  (f_2*e_2 + p_1);

  } else {

    beta2(0) = (f_2(0)*e_2(0) + p_1_NDH)/(f_2(0) + p_2_NDH);

  }
  
      
}

