
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double log_like_C(arma::colvec par, arma::mat X, arma::colvec y){
  
  arma::colvec eta = arma::zeros<arma::colvec>(X.n_rows);
  eta = X*par;
  double log_like = sum(y%(eta) - arma::log(1.0 + arma::exp(eta)));

  return(-log_like);
  
}

// [[Rcpp::export]]
arma::colvec gradient_C(arma::colvec par, arma::mat X, arma::colvec y){

  arma::colvec eta = arma::zeros<arma::colvec>(X.n_rows);
  arma::colvec out = arma::zeros<arma::colvec>(X.n_cols);
  eta = X*par;
  arma::colvec p = 1.0 / ( 1.0 + arma::exp(-1.0*eta));
  out = X.t()*(y - p); 
  
  return(-out);
  
}

// [[Rcpp::export]]
void compute_dist(arma::mat U, arma::mat& distances, std::string model, arma::mat X, arma::mat indices, bool downsampling){
  
  int N = U.n_rows;
  int index = 0;
  
  if (model == "NDH") {
  
    if (downsampling){
      int N_index = indices.n_rows;
     
      for(int l = 0; l < N_index; l++){
        int i = indices(l, 0);
        int j = indices(l, 1);
        arma::rowvec temp = U.row(i) - U.row(j);
        arma::rowvec cross_prod = temp * temp.t();
        distances.row(index) = cross_prod;
        index = index + 1;
      }
   
    } else {
      for(int i = 0; i < (N-1); i++){
       for(int j = (i+1); j < N; j++){
         arma::rowvec temp = U.row(i) - U.row(j);
         arma::rowvec cross_prod = temp * temp.t();
         distances.row(index) = cross_prod;
         index = index + 1;
       }
      }
    }
    
  } else if (model == "RS"){
  
     if (downsampling){
      int N_index = indices.n_rows;
     
      for(int l = 0; l < N_index; l++){
        int i = indices(l, 0);
        int j = indices(l, 1);
         arma::rowvec temp = U.row(i) - U.row(j);
         arma::rowvec cross_prod = temp * temp.t();
         distances.row(index) = arma::join_rows(X.row(i) + X.row(j),
                                              cross_prod);
         index = index + 1;
      }
     
     } else {
      for(int i = 0; i < (N-1); i++){
       for(int j = (i+1); j < N; j++){
         arma::rowvec temp = U.row(i) - U.row(j);
         arma::rowvec cross_prod = temp * temp.t();
         distances.row(index) = arma::join_rows(X.row(i) + X.row(j),
                                              cross_prod);
         index = index + 1;
       }
     }
    }
  
  } else {
        
    if (downsampling){
     int N_index = indices.n_rows;
     
      for(int l = 0; l < N_index; l++){
        int i = indices(l, 0);
        int j = indices(l, 1);
        arma::rowvec temp = U.row(i) - U.row(j);
        arma::rowvec cross_prod = temp * temp.t();
        distances.row(index) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1),
                                                X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1),
                                                cross_prod);
        index = index + 1;
      }
    } else {
      for(int i = 0; i < N; i++){
       for(int j = 0; j < N; j++){
     
         if (i != j){
        
          arma::rowvec temp = U.row(i) - U.row(j);
          arma::rowvec cross_prod = temp * temp.t();
          distances.row(index) = arma::join_rows(X.row(i).subvec(0, (X.n_cols*0.5) - 1),
                                                X.row(j).subvec(X.n_cols*0.5, X.n_cols - 1),
                                                cross_prod);
          index = index + 1;
       
         }
     
       }
     }
   
   } 
   
  }

  
}
