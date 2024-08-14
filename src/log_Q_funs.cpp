
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double log_Q(arma::sp_mat A, arma::mat U, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, 
             arma::colvec beta, arma::colvec p, arma::rowvec a, double b, double c, arma::mat G, arma::colvec nu, double e, double f,
             void * X, void * n_control, void * model){

  int K = mus.n_rows;
  int D = mus.n_cols;
  int N = U.n_rows;
  
  double sum_A = arma::accu(A);
  double p1_1 = 0.5*sum_A*beta(0);
  double p1_2 = 0;
  double p2 = 0;
  
  for(int i = 0; i < N; i++){
  
   for(int j = 0; j < N; j++){
     
       if (j != i){
       
          arma::rowvec temp = U.row(i) - U.row(j);
          arma::rowvec cross_prod = temp * temp.t();
          double eta = beta(0) - cross_prod(0);
          p1_2 = p1_2 + ( (0.5)*(-A(i,j)*cross_prod(0) - std::log(1.0 + std::exp(eta)) ) );
          
       }

   }
   
   for (int k = 0; k < K; k++){
     arma::rowvec temp = U.row(i) - mus.row(k);
     arma::rowvec cross_prod = temp * omegas.slice(k) * temp.t(); 
     p2  = p2 + ( prob_matrix(i,k)*(std::log(p(k)) - ( (0.5*D)*std::log(2*arma::datum::pi) ) + 
                                   ( (0.5)*std::log(arma::det(omegas.slice(k))) ) -
                                   ( (0.5)*cross_prod(0) ) ) );
   }
   
  }
  
  double p1 = p1_1 + p1_2;
  
  double p3 = 0;
  for (int k = 0; k < K; k++){
  
    arma::rowvec temp = mus.row(k) - a;
    arma::rowvec cross_prod = temp * (b*omegas.slice(k)) * temp.t(); 
    double log_det = std::log(arma::det(omegas.slice(k))); 
    
    p3 = p3 + ( ( (nu(k) - 1.0)*std::log(p(k)) ) +  ( (0.5)*log_det) -
                ( (0.5)*cross_prod(0) ) + ( ((c-D-1.0)*0.5)*log_det ) -
                ( (0.5)*arma::trace(omegas.slice(k)*G) )
              );
              
  }
  
  double log_Q = p1 + p2 + p3 + ( (-f*0.5)*std::pow((beta(0)-e), 2.0) );
  
  return log_Q;
  
}

// [[Rcpp::export]]
double log_Q_CC(arma::sp_mat A, arma::mat U, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, 
                arma::colvec beta, arma::colvec p, arma::rowvec a, double b, double c,
                arma::mat G, arma::colvec nu, double e, double f,
                double n_control, void * X, void * model){

  int K = mus.n_rows;
  int D = mus.n_cols;
  int N = U.n_rows;

  double p1 = 0;
  double p2 = 0;
  
  for(int i = 0; i < N; i++){
  
    double p1_1 = 0;
    double p1_2 = 0;
    
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
      double eta = beta(0) - cross_prod(0);
      p1_1 = p1_1 + (eta - std::log(1.0 + std::exp(eta)));
    
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
      p1_2 = p1_2 + ((N_A_0_i)/(n_A_0_i*1.0))*(-1.0*std::log(1.0 + eta_exp));
    
    }
  
   p1 = p1 + (0.5*p1_1) + (0.5*p1_2);

   for (int k = 0; k < K; k++){
     arma::rowvec temp = U.row(i) - mus.row(k);
     arma::rowvec cross_prod = temp * omegas.slice(k) * temp.t(); 
     p2  = p2 + ( prob_matrix(i,k)*(std::log(p(k)) - ( (0.5*D)*std::log(2*arma::datum::pi) ) + 
                                   ( (0.5)*std::log(arma::det(omegas.slice(k))) ) -
                                   ( (0.5)*cross_prod(0) ) ) );
   }
   
  }
  

  double p3 = 0;
  for (int k = 0; k < K; k++){
  
    arma::rowvec temp = mus.row(k) - a;
    arma::rowvec cross_prod = temp * (b*omegas.slice(k)) * temp.t(); 
    double log_det = std::log(arma::det(omegas.slice(k))); 
    
    p3 = p3 + ( ( (nu(k) - 1.0)*std::log(p(k)) ) +  ( (0.5)*log_det) -
                ( (0.5)*cross_prod(0) ) + ( ((c-D-1.0)*0.5)*log_det ) -
                ( (0.5)*arma::trace(omegas.slice(k)*G) )
              );
              
  }
  
  double log_Q = p1 + p2 + p3 + ( (-f*0.5)*std::pow((beta(0)-e), 2.0) );
  
  return log_Q;
  
}

// [[Rcpp::export]]
double log_Q_RE(arma::sp_mat A, arma::mat U, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, 
               arma::colvec beta, arma::colvec p, arma::rowvec a, double b, double c, arma::mat G, 
               arma::colvec nu, arma::colvec e, arma::mat f,
               arma::mat X, Rcpp::String model, void * n_control){

  int K = mus.n_rows;
  int D = mus.n_cols;
  int N = U.n_rows;

  double p1 = 0;
  double p2 = 0;
  
  for(int i = 0; i < N; i++){
  
   for(int j = 0; j < N; j++){
     
       if (j != i){
       
         arma::rowvec x_ij = arma::ones<arma::rowvec>(1+X.n_cols);
         if (model == "RS"){
           x_ij(arma::span(1, X.n_cols)) = X.row(i) + X.row(j);
         } else {
           x_ij(arma::span(1, X.n_cols)) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1), X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1));
         }
         
         arma::rowvec x_ij_beta = x_ij*beta; 
      
         arma::rowvec temp = U.row(i) - U.row(j);
         arma::rowvec cross_prod = temp * temp.t();
         
         double eta = x_ij_beta(0) - cross_prod(0);
         
         if (model == "RS"){
            p1 = p1 + ( (0.5)*( (A(i,j)*eta) - std::log(1.0 + std::exp(eta)) ) );
         } else {
            p1 = p1 + ( (1.0)*( (A(i,j)*eta) - std::log(1.0 + std::exp(eta)) ) );
         }
         
       }

   }
   
   for (int k = 0; k < K; k++){
     arma::rowvec temp = U.row(i) - mus.row(k);
     arma::rowvec cross_prod = temp * omegas.slice(k) * temp.t(); 
     p2  = p2 + ( prob_matrix(i,k)*(std::log(p(k)) - ( (0.5*D)*std::log(2*arma::datum::pi) ) + 
                                   ( (0.5)*std::log(arma::det(omegas.slice(k))) ) -
                                   ( (0.5)*cross_prod(0) ) ) );
   }
   
  }
  
  
  double p3 = 0;
  for (int k = 0; k < K; k++){
  
    arma::rowvec temp = mus.row(k) - a;
    arma::rowvec cross_prod = temp * (b*omegas.slice(k)) * temp.t(); 
    double log_det = std::log(arma::det(omegas.slice(k))); 
    
    p3 = p3 + ( ( (nu(k) - 1.0)*std::log(p(k)) ) +  ( (0.5)*log_det) -
                ( (0.5)*cross_prod(0) ) + ( ((c-D-1.0)*0.5)*log_det ) -
                ( (0.5)*arma::trace(omegas.slice(k)*G) )
              );
              
  }
  
  arma::colvec temp = beta - e;
  arma::colvec cross_prod = temp.t() * f * temp; 
  double p4 = -0.5*(cross_prod(0));
  
  double log_Q = p1 + p2 + p3 + p4;
  
  return log_Q;
  
}

// [[Rcpp::export]]
double log_Q_RE_CC(arma::sp_mat A, arma::mat U, arma::mat mus, arma::cube omegas, arma::mat prob_matrix, 
               arma::colvec beta, arma::colvec p, arma::rowvec a, double b, double c, arma::mat G,
               arma::colvec nu, arma::colvec e, arma::mat f,
               arma::mat X, Rcpp::String model, double n_control){

  int K = mus.n_rows;
  int D = mus.n_cols;
  int N = U.n_rows;
  
  double p1 = 0;
  double p2 = 0;
  
  for(int i = 0; i < N; i++){
  
    double p1_1 = 0;
    double p1_2 = 0;
    
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
      double eta = x_ij_beta(0) - cross_prod(0);
      p1_1 = p1_1 + (eta - std::log(1.0 + std::exp(eta)));
    
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
      p1_2 = p1_2 + ((N_A_0_i)/(n_A_0_i*1.0))*(-1.0*std::log(1.0 + eta_exp));
    
    }

    if (model == "RS"){
      p1 = p1 + (0.5*p1_1) + (0.5*p1_2);
    } else {
      p1 = p1 + (p1_1) + (p1_2);
    }

   for (int k = 0; k < K; k++){
     arma::rowvec temp = U.row(i) - mus.row(k);
     arma::rowvec cross_prod = temp * omegas.slice(k) * temp.t(); 
     p2  = p2 + ( prob_matrix(i,k)*(std::log(p(k)) - ( (0.5*D)*std::log(2*arma::datum::pi) ) + 
                                   ( (0.5)*std::log(arma::det(omegas.slice(k))) ) -
                                   ( (0.5)*cross_prod(0) ) ) );
    }
   
  }
  
  
  double p3 = 0;
  for (int k = 0; k < K; k++){
  
    arma::rowvec temp = mus.row(k) - a;
    arma::rowvec cross_prod = temp * (b*omegas.slice(k)) * temp.t(); 
    double log_det = std::log(arma::det(omegas.slice(k))); 
    
    p3 = p3 + ( ( (nu(k) - 1.0)*std::log(p(k)) ) +  ( (0.5)*log_det) -
                ( (0.5)*cross_prod(0) ) + ( ((c-D-1.0)*0.5)*log_det ) -
                ( (0.5)*arma::trace(omegas.slice(k)*G) )
              );
              
  }
  
  arma::colvec temp = beta - e;
  arma::colvec cross_prod = temp.t() * f * temp; 
  double p4 = -0.5*(cross_prod(0));
  
  double log_Q = p1 + p2 + p3 + p4;
  
  return log_Q;
  
}