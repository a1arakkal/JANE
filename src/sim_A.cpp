
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::sp_mat draw_A_NDH_c(arma::mat U, double beta0){
  int N = U.n_rows;
  arma::sp_mat A(N, N);
  
  for(int i = 0; i < N; i++){
   for(int j = 0; j < N; j++){
     if (i < j){
       arma::rowvec temp = U.row(i) - U.row(j);
       arma::rowvec temp2 = temp * temp.t();
       double eta = beta0 - temp2(0);
       double p = 0;
       if(eta >= 0){
         p = 1.0/(1.0 + std::exp(-1.0*eta));
       } else {
         p = std::exp(eta)/(1.0 + std::exp(eta));
       }
       if (arma::randu<double>() < p){
         A(i,j) = 1.0;
         A(j,i) = A(i,j);
       } 
     }
   }
  }
  
  return(A);
}

// [[Rcpp::export]]
arma::sp_mat draw_A_RS_c(arma::mat U, double beta0, arma::colvec s){
  int N = U.n_rows;
  arma::sp_mat A(N, N);
  
  for(int i = 0; i < N; i++){
   for(int j = 0; j < N; j++){
     if (i < j){
       arma::rowvec temp = U.row(i) - U.row(j);
       arma::rowvec temp2 = temp * temp.t();
       double eta = beta0 - temp2(0) + s(i) + s(j);
       double p = 0;
       if(eta >= 0){
         p = 1.0/(1.0 + std::exp(-1.0*eta));
       } else {
         p = std::exp(eta)/(1.0 + std::exp(eta));
       }
       if (arma::randu<double>() < p){
         A(i,j) = 1.0;
         A(j,i) = A(i,j);
       } 
     }
   }
  }
  
  return(A);
}

// [[Rcpp::export]]
arma::sp_mat draw_A_RSR_c(arma::mat U, double beta0, arma::colvec s, arma::colvec r){
  int N = U.n_rows;
  arma::sp_mat A(N, N);
  
  for(int i = 0; i < N; i++){
   for(int j = 0; j < N; j++){
     if (i != j){
       arma::rowvec temp = U.row(i) - U.row(j);
       arma::rowvec temp2 = temp * temp.t();
       double eta = beta0 - temp2(0) + s(i) + r(j);
       double p = 0;
       if(eta >= 0){
         p = 1.0/(1.0 + std::exp(-1.0*eta));
       } else {
         p = std::exp(eta)/(1.0 + std::exp(eta));
       }
       if (arma::randu<double>() < p){
         A(i,j) = 1.0;
       } 
     }
   }
  }
  
  return(A);
}


// [[Rcpp::export]]
void compute_mean_edge_weight(arma::mat& edge_indices, double beta0, arma::mat RE, Rcpp::String model){

  int N = edge_indices.n_rows;
  
  for(int n = 0; n < N; n++){

    int i = edge_indices(n, 0) - 1.0;
    int j = edge_indices(n, 1) - 1.0;
    
  if(model == "RS") {

    edge_indices(n, 2) = beta0 + RE(i, 0) + RE(j, 0);

  } else {

    edge_indices(n, 2) = beta0 + RE(i, 0) + RE(j, 1);

  }
    
 }

}
