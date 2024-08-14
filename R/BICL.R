
#' @export
BICL <- function(A, object){
  
  if(object$model == "NDH"){
    BIC_logit <- BIC_logit_NDH(A, object)
  } else if(object$model == "RS" ){
    BIC_logit <- BIC_logit_RS(A, object)
  } else {
    BIC_logit <- BIC_logit_RSR(A, object)
  }
  
  BIC_ICL_MBC_res <- BIC_ICL_MBC(object)
  
  out <- list(BIC_logit = BIC_logit,
              BIC_mbc = BIC_ICL_MBC_res$BIC_mbc,
              ICL_mbc = BIC_ICL_MBC_res$ICL_mbc,
              Total_BIC = BIC_logit + BIC_ICL_MBC_res$BIC_mbc,
              Total_ICL = BIC_logit + BIC_ICL_MBC_res$ICL_mbc)
  
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
