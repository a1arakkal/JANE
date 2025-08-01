
#include <RcppArmadillo.h>
#include "helper_funs.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void update_U(arma::mat& U, arma::sp_mat A, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, void * X, void * n_control, void * model){

  int N = U.n_rows;
  int D = U.n_cols;
  int K = prob_matrix.n_cols;

  for(int i = 0; i < N; i++){
    
    arma::mat p1_1(D, D, arma::fill::zeros);
    arma::colvec p2_1 = arma::zeros<arma::colvec>(D);
    
    for(int k = 0; k < K; k++){
      arma::mat temp = prob_matrix(i,k)*omegas.slice(k);
      p1_1 = p1_1 + temp;
      p2_1 = p2_1 + (temp * mus.row(k).t());
    }
    
    
   arma::mat p1_2(D, D, arma::fill::zeros);  
   arma::colvec p2_2 = arma::zeros<arma::colvec>(D);
   arma::colvec p2_3 = arma::zeros<arma::colvec>(D);
   
   for(int j = 0; j < N; j++){
   
     if (j != i){
     
       p1_2 = p1_2 + (2.0*A(i,j))*arma::eye<arma::mat>(D, D);
       
       p2_2 = p2_2 + (2.0*A(i,j))*U.row(j).t();
       
       arma::rowvec diff = U.row(i) - U.row(j);
       arma::rowvec cross_prod = diff * diff.t();
       double eta = beta(0) - cross_prod(0);
       double p = logit_inv(eta);
       
       p2_3 = p2_3 + 2.0*p*diff.t();
    
     }
     
   }

   arma::colvec new_ui = arma::solve((p1_1 + p1_2), (p2_1 + p2_2 + p2_3));

   U.row(i) = new_ui.t();
     
  }
  
}

// [[Rcpp::export]]
void update_U_CC(arma::mat& U, double n_control, arma::sp_mat A, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, void * X, void * model){

  int N = U.n_rows;
  int D = U.n_cols;
  int K = prob_matrix.n_cols;

  for(int i = 0; i < N; i++){
    
    arma::mat p1_1(D, D, arma::fill::zeros);
    arma::colvec p2_1 = arma::zeros<arma::colvec>(D);
    
    for(int k = 0; k < K; k++){
      arma::mat temp = prob_matrix(i,k)*omegas.slice(k);
      p1_1 = p1_1 + temp;
      p2_1 = p2_1 + (temp * mus.row(k).t());
    }
    
   
   // Identify indices of 1.0s 
    
    arma::sp_mat A_i = A.row(i);   
    arma::sp_mat::const_iterator start = A_i.begin();
    arma::sp_mat::const_iterator end = A_i.end();
    
    int n_its = std::distance(start, end);
    
    // Build a location storage matrix
    arma::uvec A_1_red(n_its);
    
    // Start collecting locations
    arma::sp_mat::const_iterator it = start; 
    
    for(int ii = 0; ii < n_its; ++ii){
      A_1_red(ii) = it.col();
      ++it; // increment
    }
   
   int N_A_1_i = A_1_red.n_elem;
   arma::mat p1_2 =  (2.0*arma::accu(A_i))*arma::eye<arma::mat>(D, D); 
   arma::colvec p2_2 = 2.0 * arma::sum( ( A_i.cols(A_1_red).t()*arma::ones<arma::rowvec>(2.0) ) % U.rows(A_1_red), 0).t();
   arma::colvec p2_3 = arma::zeros<arma::colvec>(D);
   arma::colvec p2_4 = arma::zeros<arma::colvec>(D);   
   
   for(int l = 0; l < N_A_1_i; l++){
       
       double j = A_1_red(l);
       
       arma::rowvec diff = U.row(i) - U.row(j);
       arma::rowvec cross_prod = diff * diff.t();
       double eta = beta(0) - cross_prod(0);
       double p = logit_inv(eta);
       p2_3 = p2_3 + 2.0*p*diff.t();
  
   }
   
   // Randomly sample 0.0s
    
   arma::colvec temp = Rcpp::runif(round(n_control * 1.1), 0.0, N - 1.0);
   arma::colvec temp2 = arma::round(temp);
   arma::colvec A_0_red = temp2(find(temp2 != i));
    
   for (size_t j = 0; j < A_1_red.n_elem; j++) {
      arma::uvec q1 = arma::find(A_0_red == A_1_red(j));
      if (!q1.empty()) {
          A_0_red.shed_rows(q1);
      }
   }
    
   if (A_0_red.n_elem > 0.0){
     A_0_red = A_0_red.subvec(0, std::min(n_control - 1.0, (A_0_red.n_elem*1.0) - 1.0));
   } 

   int n_A_0_i = A_0_red.n_elem;
   double N_A_0_i = ((N-1.0)*1.0) - (N_A_1_i*1.0);
   
   for(int l = 0; l < n_A_0_i; l++){
       
       double j = A_0_red(l);
       
       arma::rowvec diff = U.row(i) - U.row(j);
       arma::rowvec cross_prod = diff * diff.t();
       double eta = beta(0) - cross_prod(0);
       double p = logit_inv(eta);       
       p2_4 = p2_4 + ((N_A_0_i)/(n_A_0_i*1.0))*2.0*p*diff.t();
  
     
   }

   arma::colvec new_ui = arma::solve((p1_1 + p1_2), (p2_1 + p2_2 + p2_3 + p2_4));
   
   U.row(i) = new_ui.t();
     
  }
  
}

// [[Rcpp::export]]
void update_U_RE(arma::mat& U, arma::sp_mat A, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, arma::mat X, Rcpp::String model, void * n_control){

  int N = U.n_rows;
  int D = U.n_cols;
  int K = prob_matrix.n_cols;

  for(int i = 0; i < N; i++){
    
    arma::mat p1_1(D, D, arma::fill::zeros);
    arma::colvec p2_1 = arma::zeros<arma::colvec>(D);
    
    for(int k = 0; k < K; k++){
      arma::mat temp = prob_matrix(i,k)*omegas.slice(k);
      p1_1 = p1_1 + temp;
      p2_1 = p2_1 + (temp * mus.row(k).t());
    }
    
    
   arma::mat p1_2(D, D, arma::fill::zeros);  
   arma::colvec p2_2 = arma::zeros<arma::colvec>(D);
   arma::colvec p2_3 = arma::zeros<arma::colvec>(D);

   for(int j = 0; j < N; j++){
     
     if (j != i){
       
         p1_2 = p1_2 + (2.0*A(i,j))*arma::eye<arma::mat>(D, D);
       
         p2_2 = p2_2 + (2.0*A(i,j))*U.row(j).t();
       
         arma::rowvec diff = U.row(i) - U.row(j);
         arma::rowvec cross_prod = diff * diff.t();
       
         arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
         
         if (model == "RS"){
           x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
         } else {
           x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
         }
  
         arma::rowvec x_ij_beta = x_ij*beta; 
       
         double eta = x_ij_beta(0)- cross_prod(0);
         double p = logit_inv(eta);
       
         p2_3 = p2_3 + 2.0*p*diff.t();
         
         if (model != "RS"){
           
           p1_2 = p1_2 + (2.0*A(j,i))*arma::eye<arma::mat>(D, D);

           p2_2 = p2_2 + (2.0*A(j,i))*U.row(j).t();
         
           diff = U.row(i) - U.row(j);
           cross_prod = diff * diff.t();
       
           x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
           x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(j).subvec(0, (X.n_cols*0.5) - 1), X.row(i).subvec(X.n_cols*0.5, X.n_cols - 1));
           x_ij_beta = x_ij*beta; 
           eta = x_ij_beta(0)- cross_prod(0);
           p = logit_inv(eta);
           p2_3 = p2_3 + 2.0*p*diff.t();

       }
       
     }
     
   }

   arma::colvec new_ui = arma::solve((p1_1 + p1_2), ( p2_1 + p2_2  + p2_3 ));
   
   U.row(i) = new_ui.t();
     
  }
  
}

// [[Rcpp::export]]
void update_U_RE_CC(arma::mat& U, double n_control, arma::sp_mat A, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, arma::colvec beta, arma::mat X, Rcpp::String model){

  int N = U.n_rows;
  int D = U.n_cols;
  int K = prob_matrix.n_cols;

  for(int i = 0; i < N; i++){
    
    arma::mat p1_1(D, D, arma::fill::zeros);
    arma::colvec p2_1 = arma::zeros<arma::colvec>(D);
    
    for(int k = 0; k < K; k++){
      arma::mat temp = prob_matrix(i,k)*omegas.slice(k);
      p1_1 = p1_1 + temp;
      p2_1 = p2_1 + (temp * mus.row(k).t());
    }
    
   // Identify row indices of 1.0s 
    
    arma::sp_mat A_i = A.row(i);   
    arma::sp_mat::const_iterator start = A_i.begin();
    arma::sp_mat::const_iterator end = A_i.end();
    
    int n_its = std::distance(start, end);
    
    // Build a location storage matrix
    arma::uvec A_1_red(n_its);
    
    // Start collecting locations
    arma::sp_mat::const_iterator it = start; 
    
    for(int ii = 0; ii < n_its; ++ii){
      A_1_red(ii) = it.col();
      ++it; // increment
    }
   
    int N_A_1_i = A_1_red.n_elem;
    arma::mat p1_2 =  (2.0*arma::accu(A_i))*arma::eye<arma::mat>(D, D); 
    arma::colvec p2_2 = 2.0 * arma::sum( ( A_i.cols(A_1_red).t()*arma::ones<arma::rowvec>(2.0) ) % U.rows(A_1_red), 0).t();
    arma::colvec p2_3 = arma::zeros<arma::colvec>(D);
    arma::colvec p2_4 = arma::zeros<arma::colvec>(D);   
   
    for(int l = 0; l < N_A_1_i; l++){
       
       double j = A_1_red(l);
       
       arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
       
       if (model == "RS"){
           x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
       } else {
           x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
       }
         
       arma::rowvec x_ij_beta = x_ij*beta; 
       
       arma::rowvec diff = U.row(i) - U.row(j);
       arma::rowvec cross_prod = diff * diff.t();
       
       double eta = x_ij_beta(0)- cross_prod(0);
       double p = logit_inv(eta);
       p2_3 = p2_3 + 2.0*p*diff.t();
  
    }
   
    // Randomly sample 0.0s in rows
  
    arma::colvec random_samp = Rcpp::runif(round(n_control * 1.1), 0.0, N - 1.0);
    random_samp = arma::round(random_samp);
    arma::colvec A_0_red = random_samp(find(random_samp != i));
    
    for (size_t j = 0; j < A_1_red.n_elem; j++) {
      arma::uvec q1 = arma::find(A_0_red == A_1_red(j));
      if (!q1.empty()) {
          A_0_red.shed_rows(q1);
      }
    }
    
    if (A_0_red.n_elem > 0.0){
     A_0_red = A_0_red.subvec(0, std::min(n_control - 1.0, (A_0_red.n_elem*1.0) - 1.0));
    } 
   
    int n_A_0_i = A_0_red.n_elem;
    double N_A_0_i = ((N-1.0)*1.0) - (N_A_1_i*1.0);
   
    for(int l = 0; l < n_A_0_i; l++){
       
       double j = A_0_red(l);
       
       arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
       
       if (model == "RS"){
           x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
       } else {
           x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
       }
       
       arma::rowvec x_ij_beta = x_ij*beta; 
       
       arma::rowvec diff = U.row(i) - U.row(j);
       arma::rowvec cross_prod = diff * diff.t();
       
       double eta = x_ij_beta(0)- cross_prod(0);
       double p = logit_inv(eta);
       p2_4 = p2_4 + ((N_A_0_i)/(n_A_0_i*1.0))*2.0*p*diff.t();
  
    }
   
    if(model != "RS"){
   
      // Identify col indices of 1.0s 
    
      arma::sp_mat A_i_col = A.col(i).t();   
      start = A_i_col.begin();
      end = A_i_col.end();
    
      n_its = std::distance(start, end);
    
      // Build a location storage matrix
      arma::uvec A_1_red_col(n_its);
     
      // Start collecting locations
      it = start; 
    
      for(int ii = 0; ii < n_its; ++ii){
        A_1_red_col(ii) = it.col();
        ++it; // increment
      }
   
      int N_A_1_col_i = A_1_red_col.n_elem;
     
      p1_2 = p1_2 + ((2.0*arma::accu(A_i_col))*arma::eye<arma::mat>(D, D)); 
      p2_2 = p2_2 + (2.0 * arma::sum( ( A_i_col.cols(A_1_red_col).t()*arma::ones<arma::rowvec>(2.0) ) % U.rows(A_1_red_col), 0).t());
    
      for(int l = 0; l < N_A_1_col_i; l++){
       
        double j = A_1_red_col(l);
       
        arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
        x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(j).subvec(0, (X.n_cols*0.5) - 1), X.row(i).subvec(X.n_cols*0.5, X.n_cols - 1));
     
        arma::rowvec x_ij_beta = x_ij*beta; 
       
        arma::rowvec diff = U.row(i) - U.row(j);
        arma::rowvec cross_prod = diff * diff.t();
       
        double eta = x_ij_beta(0)- cross_prod(0);
        double p = logit_inv(eta);
        p2_3 = p2_3 + 2.0*p*diff.t();
  
      }
      
      // Randomly sample 0.0s in cols 
      
      arma::colvec A_0_red_col = random_samp(find(random_samp != i));
     
      for (size_t j = 0; j < A_1_red_col.n_elem; j++) {
        arma::uvec q1 = arma::find(A_0_red_col ==  A_1_red_col(j));
        if (!q1.empty()) {
          A_0_red_col.shed_rows(q1);
        }
      }
    
      if (A_0_red_col.n_elem > 0.0){
        A_0_red_col = A_0_red_col.subvec(0, std::min(n_control - 1.0, (A_0_red_col.n_elem*1.0) - 1.0));
      } 

      n_A_0_i = A_0_red_col.n_elem;
      N_A_0_i = ((N-1.0)*1.0) - (N_A_1_col_i*1.0);
   
      for(int l = 0; l < n_A_0_i; l++){
       
        double j = A_0_red_col(l);
       
        arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
        x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(j).subvec(0, (X.n_cols*0.5) - 1), X.row(i).subvec(X.n_cols*0.5, X.n_cols - 1));
       
        arma::rowvec x_ij_beta = x_ij*beta; 
       
        arma::rowvec diff = U.row(i) - U.row(j);
        arma::rowvec cross_prod = diff * diff.t();
       
        double eta = x_ij_beta(0)- cross_prod(0);
        double p = logit_inv(eta);
       
        p2_4 = p2_4 + ((N_A_0_i)/(n_A_0_i*1.0))*2.0*p*diff.t();
  
     }
    
   }

   arma::colvec new_ui = arma::solve((p1_1 + p1_2), (p2_1 + p2_2 + p2_3 + p2_4));
   
   U.row(i) = new_ui.t();
     
  }
  
}