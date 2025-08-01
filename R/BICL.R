
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
