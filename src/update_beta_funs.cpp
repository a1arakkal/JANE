
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void update_beta(arma::colvec& beta, arma::sp_mat A , arma::mat U, double f, double e, void * X, void * n_control, void * model){

  int N = U.n_rows;
  
  double sum_A = arma::accu(A);
  double p_1 = 0;
  double p_2 = 0;
  
  for(int i = 0; i < N; i++){
  
   for(int j = 0; j < N; j++){
     
       if (j > i){
       
          arma::rowvec temp = U.row(i) - U.row(j);
          arma::rowvec cross_prod = temp * temp.t();
          double eta_exp = std::exp(beta(0) - cross_prod(0));
          double p = eta_exp/(1.0 + eta_exp);
          double p_1_p = p*(1.0 - p);
          p_1 = p_1 + (p_1_p*beta(0) - p);
          p_2 = p_2 + p_1_p;
       }
       
   }
   
  }
  
  beta(0) = (f*e + (0.5*sum_A) + (1.0*p_1))/(f + (1.0*p_2));
  
}

// [[Rcpp::export]]
void update_beta_CC(arma::colvec& beta, arma::sp_mat A, double n_control, arma::mat U, double f, double e, void * X, void * model){

   int N = U.n_rows;
  
   double sum_A = arma::accu(A);
   double p_1 = 0;
   double p_2 = 0;
   double p_3 = 0;
   double p_4 = 0;
  
   for(int i = 0; i < N; i++){
  
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

    for(int l = 0; l < N_A_1_i; l++){
     
      double j = A_1_red(l);
      arma::rowvec temp = U.row(i) - U.row(j);
      arma::rowvec cross_prod = temp * temp.t();
      double eta_exp = std::exp(beta(0) - cross_prod(0));
      double p = eta_exp/(1.0 + eta_exp);
      double p_1_p = p*(1.0 - p);
      p_1 = p_1 + (p_1_p*beta(0) - p);
      p_2 = p_2 + p_1_p;

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
      arma::rowvec temp = U.row(i) - U.row(j);
      arma::rowvec cross_prod = temp * temp.t();
      double eta_exp = std::exp(beta(0) - cross_prod(0));
      double p = eta_exp/(1.0 + eta_exp);
      double p_1_p = p*(1.0 - p);
      p_3 = p_3 + ((N_A_0_i)/(n_A_0_i*1.0))*(p_1_p*beta(0) - p);
      p_4 = p_4 + ((N_A_0_i)/(n_A_0_i*1.0))*p_1_p;

   }
   
  }
  
  beta(0) = (f*e + (0.5*sum_A) + (0.5*p_1) + (0.5*p_3))/(f + (0.5*p_2) + (0.5*p_4));
  
}


// [[Rcpp::export]]
void update_beta_RE(arma::colvec& beta, arma::sp_mat A, arma::mat U, arma::mat f, arma::colvec e, arma::mat X, Rcpp::String model, void * n_control){

  int N = U.n_rows;
 
  arma::colvec p_1 = arma::zeros<arma::colvec>(1+X.n_cols);
  arma::mat p_2(1+X.n_cols, 1+X.n_cols, arma::fill::zeros);
  arma::colvec p_3 = arma::zeros<arma::colvec>(1+X.n_cols);
  
  for(int i = 0; i < N; i++){
   
      if (model == "RS"){
      
        for(int j = 0; j < N; j++){
           
          if(j > i){
          
            arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
            x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
            arma::rowvec x_ij_beta = x_ij*beta; 
            arma::rowvec temp = U.row(i) - U.row(j);
            arma::rowvec cross_prod = temp * temp.t();
            double eta_exp = std::exp(x_ij_beta(0) - cross_prod(0));
            double p = eta_exp/(1.0 + eta_exp);
            arma::mat p_1_p = p*(1.0 - p)*(x_ij.t() * x_ij);
            p_1 = p_1 + (p_1_p*beta - (p*x_ij.t()));
            p_2 = p_2 + p_1_p;
            p_3 = p_3 + A(i,j)*(x_ij.t());
            
          }
         
        }
      
      } else {
      
        for(int j = 0; j < N; j++){
     
          if (j != i){
        
           arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
           x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
           arma::rowvec x_ij_beta = x_ij*beta; 
           arma::rowvec temp = U.row(i) - U.row(j);
           arma::rowvec cross_prod = temp * temp.t();
           double eta_exp = std::exp(x_ij_beta(0) - cross_prod(0));
           double p = eta_exp/(1.0 + eta_exp);
           arma::mat p_1_p = p*(1.0 - p)*(x_ij.t() * x_ij);
           p_1 = p_1 + (p_1_p*beta - (p*x_ij.t()));
           p_2 = p_2 + p_1_p;
           p_3 = p_3 + A(i,j)*(x_ij.t());
          
         }
        
      }

     }
   
  }
       
  beta = arma::inv_sympd(f + (p_2)) *  (f*e + (p_3) + (p_1));
      
}


// [[Rcpp::export]]
void update_beta_RE_CC(arma::colvec& beta, arma::sp_mat A, double n_control, arma::mat U, arma::mat f, arma::colvec e, arma::mat X, Rcpp::String model){

   int N = U.n_rows;

   arma::colvec p_1 = arma::zeros<arma::colvec>(1+X.n_cols);
   arma::mat p_2(1+X.n_cols, 1+X.n_cols, arma::fill::zeros);
   arma::colvec p_3 = arma::zeros<arma::colvec>(1+X.n_cols);
   arma::colvec p_4 = arma::zeros<arma::colvec>(1+X.n_cols);
   arma::mat p_5(1+X.n_cols, 1+X.n_cols, arma::fill::zeros);
  
   for(int i = 0; i < N; i++){
  
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

    for(int l = 0; l < N_A_1_i; l++){
     
      double j = A_1_red(l);
      
      arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
      
      if (model == "RS"){
        x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
      } else {
        x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
      }
      
      arma::rowvec x_ij_beta = x_ij*beta; 
          
      arma::rowvec temp = U.row(i) - U.row(j);
      arma::rowvec cross_prod = temp * temp.t();
      double eta_exp = std::exp(x_ij_beta(0) - cross_prod(0));
      double p = eta_exp/(1.0 + eta_exp);
      arma::mat p_1_p = p*(1.0 - p)*(x_ij.t() * x_ij);
      
      p_1 = p_1 + (p_1_p*beta - (p*x_ij.t()));
      p_2 = p_2 + p_1_p;
      p_3 = p_3 + (A(i,j)*(x_ij.t()));

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
      
      arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
      
      if (model == "RS"){
        x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
      } else {
        x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
      }
      
      arma::rowvec x_ij_beta = x_ij*beta; 
          
      arma::rowvec temp = U.row(i) - U.row(j);
      arma::rowvec cross_prod = temp * temp.t();
      double eta_exp = std::exp(x_ij_beta(0) - cross_prod(0));
      double p = eta_exp/(1.0 + eta_exp);
      arma::mat p_1_p = p*(1.0 - p)*(x_ij.t() * x_ij);
 
      p_4 = p_4 + ((N_A_0_i)/(n_A_0_i*1.0))*(p_1_p*beta - (p*x_ij.t()));
      p_5 = p_5 + ((N_A_0_i)/(n_A_0_i*1.0))*p_1_p;

   }
   
  }
  
   if (model == "RS"){
   
     beta = arma::inv_sympd(f + (0.5*p_2) + (0.5*p_5)) *  (f*e + (0.5*p_3) + (0.5*p_1) + (0.5*p_4));
      
   } else {

     beta = arma::inv_sympd(f + (p_2) + (p_5)) *  (f*e + (p_3) + (p_1) + (p_4));
      
   }
   
}
