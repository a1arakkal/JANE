
BICL <- function(A, object){
  
  if(!object$noise_weights){
    
    if(object$model == "NDH"){
      BIC_model <- BIC_logit_NDH(A, object)
    } else if(object$model == "RS" ){
      BIC_model <- BIC_logit_RS(A, object)
    } else {
      BIC_model <- BIC_logit_RSR(A, object)
    }
    
  } else {
    
    A[object$prob_matrix_W[, 1:2, drop = FALSE]] <- object$prob_matrix_W[, 3]
    
    if(object$model != "RSR"){
      A[object$prob_matrix_W[, 2:1, drop = FALSE]] <- object$prob_matrix_W[, 3]
    }
    
    BIC_model <- BIC_hurdle(A, object)
    
  }
  
  BIC_ICL_MBC_res <- BIC_ICL_MBC(object)
  
  out <- list(BIC_model = BIC_model,
              BIC_mbc = BIC_ICL_MBC_res$BIC_mbc,
              ICL_mbc = BIC_ICL_MBC_res$ICL_mbc,
              Total_BIC = BIC_model + BIC_ICL_MBC_res$BIC_mbc,
              Total_ICL = BIC_model + BIC_ICL_MBC_res$ICL_mbc)
  
  return(out)
}

#' @useDynLib JANE  
BIC_logit_NDH <- function(A, object) {
.Call('_JANE_BIC_logit_NDH', PACKAGE = 'JANE', A, object)
}

#' @useDynLib JANE  
BIC_logit_RS <- function(A, object) {
  .Call('_JANE_BIC_logit_RS', PACKAGE = 'JANE', A, object)
}

#' @useDynLib JANE  
BIC_logit_RSR <- function(A, object) {
  .Call('_JANE_BIC_logit_RSR', PACKAGE = 'JANE', A, object)
}

#' @useDynLib JANE  
BIC_ICL_MBC <- function(object) {
  .Call('_JANE_BIC_ICL_MBC', PACKAGE = 'JANE', object)
}

#' @useDynLib JANE  
BIC_hurdle <- function(W, object) {
  .Call('_JANE_BIC_hurdle', PACKAGE = 'JANE', W, object)
}

#' @useDynLib JANE  
trunc_poisson_density_BIC <- function(w, mean, log) {
  .Call('_JANE_trunc_poisson_density_BIC', PACKAGE = 'JANE', w, mean, log)
}

#' @useDynLib JANE  
lognormal_density_BIC <- function(w, precision, mean, log) {
  .Call('_JANE_lognormal_density_BIC', PACKAGE = 'JANE', w, precision, mean, log)
}
