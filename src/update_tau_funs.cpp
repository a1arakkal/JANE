
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void update_precision_weights(arma::colvec& precision_weights, arma::mat prob_matrix_W, Rcpp::String model, arma::colvec beta2, arma::mat X2, double m_1, double o_1, arma::mat f_2, arma::colvec e_2){

  int M = prob_matrix_W.n_rows;
  double p1 = arma::sum(prob_matrix_W.col(3));
  double p2 = arma::sum(prob_matrix_W.col(3) % arma::pow(arma::log(prob_matrix_W.col(2)), 2));
  arma::mat p3(1+X2.n_cols, 1+X2.n_cols, arma::fill::zeros);

  if (model == "NDH"){

    precision_weights(0) = (p1 + m_1 - 1.0) / (p2 - std::pow(beta2(0), 2)*(p1 + f_2(0)) + (f_2(0) * std::pow(e_2(0), 2)) + o_1);
  
  } else {

    for(int m = 0; m < M; m++){
       
      int i = prob_matrix_W(m, 0) - 1.0;
      int j = prob_matrix_W(m, 1) - 1.0;
        
      arma::rowvec x2_ij = arma::ones<arma::rowvec>(1+X2.n_cols);

      if (model == "RS"){

        x2_ij(arma::span(1, X2.n_cols)) = X2.row(i) + X2.row(j);

      } else {

        x2_ij(arma::span(1, X2.n_cols)) = arma::join_rows(X2.row(i).subvec(0, (X2.n_cols*0.5) - 1), X2.row(j).subvec(X2.n_cols*0.5, X2.n_cols - 1));

      }
      
      p3 = p3 +  prob_matrix_W(m, 3)*x2_ij.t()*x2_ij;


    }

    arma::rowvec temp1 = beta2.t()*(p3 + f_2)*beta2;
    arma::rowvec temp2 = e_2.t()*f_2*e_2;

    precision_weights(0) = (p1 + m_1 + beta2.n_elem - 2.0)/(p2 - temp1(0) + temp2(0) + o_1);

  }

}

// [[Rcpp::export]]
void update_precision_noise_weights(arma::colvec& precision_noise_weights, arma::mat prob_matrix_W, double guess_noise_weights, double m_2, double o_2){

  double p1 = arma::sum(prob_matrix_W.col(4));
  double p2 = arma::sum(prob_matrix_W.col(4) % arma::pow(arma::log(prob_matrix_W.col(2)) - guess_noise_weights*arma::ones<arma::colvec>(prob_matrix_W.n_rows), 2));
  
  precision_noise_weights(0) = (p1 + m_2 - 2.0)/(p2 + o_2);

}
