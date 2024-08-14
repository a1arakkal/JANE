#' sim_A - Simulate unweighted network
#' @param N Number of nodes/actors
#' @param mus Cluster centers (KxD matrix)
#' @param omegas Precision matrices for cluster (DxDxK array)
#' @param p_k Mixture probabilities (Kx1 vector)
#' @param beta0 Intercept parameter
#' @param model Model to simulate data from ("NDH", "RS", or "RSR)
#' @param precision_R_effects Precision for random effects (a scalar for "RS" or 2x2 matrix for "RSR")
#' @param remove_isolates A logical indicating whether or not to exclude isolates from network 
#' @return Add details
#' @details
#'    More summary needed here. 
#' @export

sim_A <- function(N, mus, omegas, p_k, beta0, 
                  model = "NDH",
                  precision_R_effects,
                  remove_isolates = F){
  
  if(!model %in% c("NDH", "RS", "RSR")){
    stop("Model needs to be one of the following: NDH, RS, or RSR")
  }
  
  # draw Z
  Z <-  t(stats::rmultinom(n = N, size = 1, prob = p_k))
  
  # draw U
  U <- matrix(stats::rnorm(N*ncol(mus)), nrow = N, ncol = ncol(mus))
  n_k <- colSums(Z)
  for (k in 1:length(n_k)){
    U[Z[, k] == 1,] <- (U[Z[, k] == 1, ] %*% t(solve(chol(omegas[,,k])))) + ( rep(1, n_k[k]) %*% t(mus[k,]) )
  }
  
  if(model == "RS"){
    
    if(!missing(precision_R_effects)){
      if(!(length(precision_R_effects) == 1 && precision_R_effects > 0)){
        stop("For RS Model, please supply a positive scalar precision value for omega_s^2")
      }
    } else {
      precision_R_effects <- stats::rgamma(n = 1, shape = 1, scale = 1)
    }
    
    # draw s
    s <- stats::rnorm(n = N, mean = 0, sd = sqrt(1/precision_R_effects))
    # draw A
    A <- draw_A_RS_c(U = U, beta0 = beta0, s = s)
    
    RE <- matrix(s, nrow = N, ncol = 1)
    colnames(RE) <- "s"
    
  } else if(model == "RSR"){
    
    if(!missing(precision_R_effects)){
      if(!(length(precision_R_effects) == 4 && all(dim(precision_R_effects) == c(2,2)) && all(eigen(precision_R_effects)$values > 0))) {
        stop("For RSR Model, please supply a 2x2 p.d. precision matrix")
      }
    } else {
      precision_R_effects <- stats::rWishart(n = 1, df  = 2 + 1, Sigma = diag(2))[,,1]
    }

    RE <- matrix(stats::rnorm(n = N*2), ncol = 2) %*% t(solve(chol(precision_R_effects)))
    
    # draw A
    A <- draw_A_RSR_c(U = U, beta0 = beta0, s = RE[,1], r = RE[,2])
    colnames(RE) <- c("s", "r")
    
  } else {
    
    # draw A
    A <- draw_A_NDH_c(U = U, beta0 = beta0)
    RE <- NULL
    precision_R_effects <- NULL
    
  }
  
  if (remove_isolates){
    
    isolates <- which(rowSums(A)==0 & colSums(A)==0)
    
    if(length(isolates) == 0){
      message("No isolates to remove")
      return(list(A = A,
                  Z = Z,
                  U = U,
                  RE = RE,
                  precision_R_effects = precision_R_effects,
                  model = model))
    } else {
      message(paste0(length(isolates), " isolate(s) removed"))
      return(list(A = A[-isolates, -isolates],
                  Z = Z[-isolates, ],
                  U = U[-isolates, ],
                  RE = RE[-isolates, , drop = F],
                  precision_R_effects = precision_R_effects,
                  model = model))
    }
    
    
  } else {
    return(list(A = A,
                Z = Z,
                U = U,
                RE = RE,
                precision_R_effects = precision_R_effects,
                model = model))
  }
  
}

#' @useDynLib JANE   
draw_A_NDH_c <- function(U, beta0) {
  .Call('_JANE_draw_A_NDH_c', PACKAGE = 'JANE', U, beta0)
}

#' @useDynLib JANE   
draw_A_RS_c <- function(U, beta0, s) {
  .Call('_JANE_draw_A_RS_c', PACKAGE = 'JANE', U, beta0, s)
}

#' @useDynLib JANE   
draw_A_RSR_c <- function(U, beta0, s, r) {
  .Call('_JANE_draw_A_RSR_c', PACKAGE = 'JANE', U, beta0, s, r)
}
