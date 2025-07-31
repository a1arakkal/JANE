
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double trunc_poisson_density_BIC(double w, double mean, double log){
  
  double temp = -1.0*arma::datum::inf;

  if (w > 0.0){

    temp = w*std::log(mean) - mean - std::lgamma(w + 1.0) - std::log(-1.0*std::expm1(-1.0*mean));

    if(log <= 0.0){

      temp = std::exp(temp);

    } 
    
  }

  return temp;

}

// [[Rcpp::export]]
double lognormal_density_BIC(double w, double precision, double mean, double log){

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
double BIC_logit_NDH(arma::sp_mat A, Rcpp::List object){
   
   arma::mat U = object["U"];
   arma::colvec beta = object["beta"];
   int N = U.n_rows;
   double n_edges = 0.5*arma::accu(A);
   double p_1 = 0;
   
   for(int i = 0; i < N; i++){
  
    for(int j = 0; j < N; j++){
     
       if (j > i){
       
          arma::rowvec temp = U.row(i) - U.row(j);
          arma::rowvec cross_prod = temp * temp.t();
          double eta = beta(0) - cross_prod(0);
          double max_val = std::max(0.0, eta);
          p_1 = p_1 + ( A(i,j)*eta - max_val - std::log( std::exp(-1.0*max_val) + std::exp(eta - max_val) ) );
          
       }
       
    }
   
   }
   
   return (-2.0*p_1) + ( beta.n_elem * std::log(n_edges));
   
}

// [[Rcpp::export]]
double BIC_logit_RS(arma::sp_mat A, Rcpp::List object){
   
   arma::mat U = object["U"];
   arma::mat X = object["X"];
   arma::colvec beta = object["beta"];
   int N = U.n_rows;
   double n_edges = 0.5*arma::accu(A);
   double p_1 = 0;
   
   for(int i = 0; i < N; i++){
  
    for(int j = 0; j < N; j++){
     
       if (j > i){
       
         arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
         x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
         arma::rowvec x_ij_beta = x_ij*beta; 
         arma::rowvec temp = U.row(i) - U.row(j);
         arma::rowvec cross_prod = temp * temp.t();
         double eta = x_ij_beta(0) - cross_prod(0);
         double max_val = std::max(0.0, eta);
         p_1 = p_1 + ( A(i,j)*eta - max_val - std::log( std::exp(-1.0*max_val) + std::exp(eta - max_val) ) );
       
       }
       
    }
   
   }
   
   return (-2.0*p_1) + ( beta.n_elem * std::log(n_edges));
   
}

// [[Rcpp::export]]
double BIC_logit_RSR(arma::sp_mat A, Rcpp::List object){
   
   arma::mat U = object["U"];
   arma::mat X = object["X"];
   arma::colvec beta = object["beta"];
   int N = U.n_rows;
   double n_edges = arma::accu(A);
   double p_1 = 0;
   
   for(int i = 0; i < N; i++){
  
    for(int j = 0; j < N; j++){
     
       if (j != i){
       
         arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
         x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
         arma::rowvec x_ij_beta = x_ij*beta; 
         arma::rowvec temp = U.row(i) - U.row(j);
         arma::rowvec cross_prod = temp * temp.t();
         double eta = x_ij_beta(0) - cross_prod(0);
         double max_val = std::max(0.0, eta);
         p_1 = p_1 + ( A(i,j)*eta - max_val - std::log( std::exp(-1.0*max_val) + std::exp(eta - max_val) ) );
       
       }
       
    }
   
   }
   
   return (-2.0*p_1) + ( beta.n_elem * std::log(n_edges));
   
}

// [[Rcpp::export]]
Rcpp::List BIC_ICL_MBC(Rcpp::List object){
   
   arma::mat U = object["U"];
   arma::mat mus = object["mus"];
   arma::cube omegas = object["omegas"];
   arma::colvec p = object["p"];
   arma::mat prob_mat = object["prob_matrix"];
   
   int N = U.n_rows;
   int D = U.n_cols;
   int K = p.n_elem;
   
   double n_params = (K-1.0) + (K*D) + ( K * ( 0.5 * (D*(D+1)) ) );
   double p_1 = 0;
   double expected_entropy = 0;
      
   for(int i = 0; i < N; i++){
    
    arma::rowvec p_vec = arma::zeros<arma::rowvec>(K);
    arma::rowvec z_star = arma::zeros<arma::rowvec>(K);
    arma::rowvec z = prob_mat.row(i);
    z_star(arma::index_max(z)) = 1.0;
    
    for(int k = 0; k < K; k++){
     
      arma::rowvec temp = U.row(i) - mus.row(k);
      arma::rowvec temp2 = -0.5 * temp * omegas.slice(k) * temp.t();
      double log_mvn = (-0.5*D)*std::log(2.0*arma::datum::pi) + 0.5*std::log(arma::det(omegas.slice(k))) + temp2(0);
      p_vec(k) = std::log(p(k)) + log_mvn;
        
      double temp_log_z_k = 0.0;
      if (z(k) > 0.0){
        temp_log_z_k = std::log(z(k));
      } 
      
      expected_entropy = expected_entropy + (z_star(k)*temp_log_z_k);
      
    }
    
     double max_val = arma::max(p_vec);
     double p_2 = std::log( arma::sum( arma::exp(p_vec - arma::ones<arma::rowvec>(K)*max_val) ) ); 
     p_1 = p_1 + max_val + p_2;
   
   }
   
   double BIC = (-2.0*p_1) + ( n_params * std::log(N));
   
   expected_entropy = -1.0*expected_entropy;
   double ICL = BIC + (2.0*(expected_entropy));
   
   return Rcpp::List::create(Rcpp::Named("BIC_mbc")=BIC,
                             Rcpp::Named("ICL_mbc")=ICL);
                             
}

// [[Rcpp::export]]
double BIC_hurdle(arma::sp_mat W, Rcpp::List object){
   
   arma::mat U = object["U"];
   arma::mat X = object["X"];
   arma::colvec beta = object["beta"];
   Rcpp::String model = object["model"];
   Rcpp::String family = object["family"];
   arma::mat prob_matrix_W = object["prob_matrix_W"];
   double q_prob = object["q_prob"];

   int N = U.n_rows;
   double n_edges = prob_matrix_W.n_rows;

   double p_1 = 0;
   
   for(int i = 0; i < N; i++){
  
    for(int j = 0; j < N; j++){

      if(model == "NDH"){
        
        if (j > i){
       
           arma::rowvec temp = U.row(i) - U.row(j);
           arma::rowvec cross_prod = temp * temp.t();
           double eta = beta(0) - cross_prod(0);
           double w = W(i,j);
           double max_val = std::max(0.0, eta);
           p_1 = p_1 + ( -max_val - std::log(std::exp(-1.0*max_val) + std::exp(eta - max_val)) );
          
           if (w > 0.0){
              
              double log_density = 0.0;
              double log_density_noise = 0.0;

              if (family != "bernoulli"){
           
                arma::colvec beta2 = object["beta2"];
                arma::mat X2 = object["X2"];
                double guess_noise_weights = object["guess_noise_weights"];
                double eta_w =  beta2(0);

                if (family == "lognormal"){

                  double precision_weights = object["precision_weights"];
                  double precision_noise_weights = object["precision_noise_weights"];
                  log_density = lognormal_density_BIC(w, precision_weights, eta_w, 1.0);
                  log_density_noise = lognormal_density_BIC(w, precision_noise_weights, guess_noise_weights, 1.0);

                } else {

                  log_density = trunc_poisson_density_BIC(w, std::exp(eta_w), 1.0);
                  log_density_noise = trunc_poisson_density_BIC(w, guess_noise_weights, 1.0);

                }

              }

              double temp1 = eta + log_density;
              double temp2 = std::log(q_prob) + log_density_noise;
              double max_val = std::max(temp1, temp2);

              p_1 = p_1 + ( max_val + std::log( std::exp(temp1 - max_val) + std::exp(temp2 - max_val) ) );

           } else {

             p_1 = p_1 + std::log(1.0 - q_prob);

           }

          
        }

      } else if (model == "RS"){

        if (j > i){
        
          arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
          x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
          arma::rowvec x_ij_beta = x_ij*beta; 
          arma::rowvec temp = U.row(i) - U.row(j);
          arma::rowvec cross_prod = temp * temp.t();
          double eta = x_ij_beta(0) - cross_prod(0);
          double w = W(i,j);
          double max_val = std::max(0.0, eta);
          p_1 = p_1 + ( -max_val - std::log(std::exp(-1.0*max_val) + std::exp(eta - max_val)) );
          
           if (w > 0.0){
              
              double log_density = 0.0;
              double log_density_noise = 0.0;

              if (family != "bernoulli"){

                arma::colvec beta2 = object["beta2"];
                arma::mat X2 = object["X2"];
                double guess_noise_weights = object["guess_noise_weights"];

                arma::rowvec x2_ij = arma::ones<arma::rowvec>(1+X2.n_cols);
                x2_ij(arma::span(1, X2.n_cols)) = X2.row(i) + X2.row(j);
                arma::rowvec x2_ij_beta = x2_ij*beta2; 
                double eta_w = x2_ij_beta(0);

                if (family == "lognormal"){

                  double precision_weights = object["precision_weights"];
                  double precision_noise_weights = object["precision_noise_weights"];
                  log_density = lognormal_density_BIC(w, precision_weights, eta_w, 1.0);
                  log_density_noise = lognormal_density_BIC(w, precision_noise_weights, guess_noise_weights, 1.0);

                } else {

                  log_density = trunc_poisson_density_BIC(w, std::exp(eta_w), 1.0);
                  log_density_noise = trunc_poisson_density_BIC(w, guess_noise_weights, 1.0);

                }

              }

              double temp1 = eta + log_density;
              double temp2 = std::log(q_prob) + log_density_noise;
              double max_val = std::max(temp1, temp2);

              p_1 = p_1 + ( max_val + std::log( std::exp(temp1 - max_val) + std::exp(temp2 - max_val) ) );

           } else {

             p_1 = p_1 + std::log(1.0 - q_prob);

           }
       
        }

      } else {

         if (j != i){
       
          arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
          x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
          arma::rowvec x_ij_beta = x_ij*beta; 
          arma::rowvec temp = U.row(i) - U.row(j);
          arma::rowvec cross_prod = temp * temp.t();
          double eta = x_ij_beta(0) - cross_prod(0);
          double w = W(i,j);
          double max_val = std::max(0.0, eta);
          p_1 = p_1 + ( -max_val - std::log(std::exp(-1.0*max_val) + std::exp(eta - max_val)) );
          
           if (w > 0.0){
              
              double log_density = 0.0;
              double log_density_noise = 0.0;

              if (family != "bernoulli"){

                arma::colvec beta2 = object["beta2"];
                arma::mat X2 = object["X2"];
                double guess_noise_weights = object["guess_noise_weights"];

                arma::rowvec x2_ij = arma::ones<arma::rowvec>(1+X2.n_cols);
                x2_ij(arma::span(1, X2.n_cols)) = arma::join_rows(X2.row(i).subvec(0, (X2.n_cols*0.5) - 1), X2.row(j).subvec(X2.n_cols*0.5, X2.n_cols - 1));
                arma::rowvec x2_ij_beta = x2_ij*beta2; 
                double eta_w = x2_ij_beta(0);

                if (family == "lognormal"){

                  double precision_weights = object["precision_weights"];
                  double precision_noise_weights = object["precision_noise_weights"];
                  log_density = lognormal_density_BIC(w, precision_weights, eta_w, 1.0);
                  log_density_noise = lognormal_density_BIC(w, precision_noise_weights, guess_noise_weights, 1.0);

                } else {

                  log_density = trunc_poisson_density_BIC(w, std::exp(eta_w), 1.0);
                  log_density_noise = trunc_poisson_density_BIC(w, guess_noise_weights, 1.0);

                }

              }

              double temp1 = eta + log_density;
              double temp2 = std::log(q_prob) + log_density_noise;
              double max_val = std::max(temp1, temp2);

              p_1 = p_1 + ( max_val + std::log( std::exp(temp1 - max_val) + std::exp(temp2 - max_val) ) );

           } else {

             p_1 = p_1 + std::log(1.0 - q_prob);

           }
             
        }

      }
            
    }
   
   }

   double n_params = 0.0;

   if (family == "bernoulli"){

      n_params = beta.n_elem + 1.0;

   } else if (family == "poisson"){

      arma::colvec beta2 = object["beta2"];
      n_params = beta.n_elem + beta2.n_elem + 1.0;

   } else {

      arma::colvec beta2 = object["beta2"];
      n_params = beta.n_elem + beta2.n_elem + 2.0 + 1.0;

   }
   
   return (-2.0*p_1) + ( (n_params) * std::log(n_edges));
   
}

