
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

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
          p_1 = p_1 + ( (A(i,j)*eta - std::log(1.0+std::exp(eta)) ));
          
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
         p_1 = p_1 + ( (A(i,j)*eta - std::log(1.0+std::exp(eta)) ));
       
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
         p_1 = p_1 + ( (A(i,j)*eta - std::log(1.0+std::exp(eta)) ));
       
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
    
    double p_2 = 0;
    arma::rowvec z_star = arma::zeros<arma::rowvec>(K);
    arma::rowvec z = prob_mat.row(i);
    z_star(arma::index_max(z)) = 1.0;
    
    for(int k = 0; k < K; k++){
     
      arma::rowvec temp = U.row(i) - mus.row(k);
      arma::rowvec temp2 = -0.5 * temp * omegas.slice(k) * temp.t();
      double den = std::pow( 2.0*arma::datum::pi , -1.0*0.5*D) * std::pow(arma::det(omegas.slice(k)), 0.5) * std::exp(temp2(0));
      p_2 = p_2 + (p(k)*den);
        
        
      double temp_log_z_k = 0.0;
      if (z(k) > 0){
        temp_log_z_k = std::log(z(k));
      } 
      
      expected_entropy = expected_entropy + (z_star(k)*temp_log_z_k);
      
    }
    
    p_1 = p_1 + std::log(p_2);
   
   }
   
   double BIC = (-2.0*p_1) + ( n_params * std::log(N));
   
   expected_entropy = -1.0*expected_entropy;
   double ICL = BIC + (2.0*(expected_entropy));
   
   return Rcpp::List::create(Rcpp::Named("BIC_mbc")=BIC,
                             Rcpp::Named("ICL_mbc")=ICL);
                             
}

// [[Rcpp::export]]
double BIC_hurdle(arma::sp_mat A, Rcpp::List object){
   
   arma::mat U = object["U"];
   arma::mat X = object["X"];
   arma::colvec beta = object["beta"];
   Rcpp::String model = object["model"];
   double q_prob = object["q_prob"];

   int N = U.n_rows;
   double n_edges = 0.0;

   if(model == "RSR"){
     n_edges = arma::accu(A);
   } else {
     n_edges = 0.5*arma::accu(A);
   }

   double p_1 = 0;
   
   for(int i = 0; i < N; i++){
  
    for(int j = 0; j < N; j++){

      if(model == "NDH"){
        
        if (j > i){
       
           arma::rowvec temp = U.row(i) - U.row(j);
           arma::rowvec cross_prod = temp * temp.t();
           double eta = beta(0) - cross_prod(0);
           p_1 = p_1 + ( ( (1.0 - A(i,j))*std::log(1.0 - q_prob) ) - std::log(1.0 + std::exp(eta)) + ( A(i,j)*std::log(std::exp(eta) + q_prob ) ) );
          
        }

      } else if (model == "RS"){

        if (j > i){
        
          arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
          x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
          arma::rowvec x_ij_beta = x_ij*beta; 
          arma::rowvec temp = U.row(i) - U.row(j);
          arma::rowvec cross_prod = temp * temp.t();
          double eta = x_ij_beta(0) - cross_prod(0);
          p_1 = p_1 + ( ( (1.0 - A(i,j))*std::log(1.0 - q_prob) ) - std::log(1.0 + std::exp(eta)) + ( A(i,j)*std::log(std::exp(eta) + q_prob ) ) );
       
        }

      } else {

         if (j != i){
       
          arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
          x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
          arma::rowvec x_ij_beta = x_ij*beta; 
          arma::rowvec temp = U.row(i) - U.row(j);
          arma::rowvec cross_prod = temp * temp.t();
          double eta = x_ij_beta(0) - cross_prod(0);
          p_1 = p_1 + ( ( (1.0 - A(i,j))*std::log(1.0 - q_prob) ) - std::log(1.0 + std::exp(eta)) + ( A(i,j)*std::log(std::exp(eta) + q_prob ) ) );
             
        }

      }
            
    }
   
   }
   
   return (-2.0*p_1) + ( (beta.n_elem + 1.0) * std::log(n_edges));
   
}

