#' Simulate unweighted or weighted networks, with or without noise edges, from latent space cluster models
#' @description Simulate an unweighted or weighted network, with or without noise edges, from a \eqn{D}-dimensional latent space cluster model with \eqn{K} clusters and \eqn{N} actors. The \emph{squared} euclidean distance is used (i.e., \eqn{dist(U_i,U_j)^2}), where \eqn{U_i} and \eqn{U_j} are the respective actor's positions in an unobserved social space.
#' @param N An integer specifying the number of actors in the network.
#' @param mus A numeric \eqn{K \times D} matrix specifying the mean vectors of the multivariate normal distribution for the latent positions of the \eqn{K} clusters.
#' @param omegas A numeric \eqn{D \times D \times K} array specifying the precision matrices of the multivariate normal distribution for the latent positions of the \eqn{K} clusters.
#' @param p A numeric vector of length \eqn{K} specifying the mixture weights of the finite multivariate normal mixture distribution for the latent positions.
#' @param model A character string to specify the model to simulate the network from:
#'  \itemize{
#'   \item{'NDH': generates an \strong{undirected} network with no degree heterogeneity}
#'   \item{'RS': generates an \strong{undirected} network with degree heterogeneity, specifically by including actor specific random sociality effects}
#'   \item{'RSR': generates a \strong{directed} network with degree heterogeneity, specifically by including actor specific random sender and receiver effects}
#'   }
#' @param family A character vector specifying the distribution of the edge weights.
#'  \itemize{
#'   \item{'bernoulli': Expects an unweighted network; utilizes a logit link (default)}
#'   \item{'lognormal': Expects a weighted network with non-negative continuous weights; utilizes an identity link}
#'   \item{'poisson': Expects a weighted network with weights representing counts; utilizes an log link}
#'   }
#' @param params_LR A list containing the parameters of the logistic regression model, including:
#'  \itemize{
#'   \item{'beta0': A numeric value specifying the intercept parameter for the logistic regression model}
#'   \item{'precision_R_effects': Precision parameters for random degree heterogeneity effects, specific to the logistic regression model:
#'     \itemize{
#'       \item{'NDH': does not apply, can leave as missing}
#'       \item{'RS': a numeric value specifying the precision parameter of the normal distribution of the random sociality effect, if missing will generate from a gamma(shape = 1, rate = 1)}
#'       \item{'RSR': a numeric matrix specifying the precision matrix of the multivariate normal distribution of the random sender and receiver effects, if missing will generate from a Wishart(df = 3, Sigma = \eqn{I_2})}
#'     } 
#'    }
#'   }
#' @param params_weights Only relevant if \code{family \%in\% c('lognormal', 'poisson')}. A list containing the parameters of the models for the edge weights, including:
#'  \itemize{
#'   \item{'beta0': A numeric value specifying the intercept parameter for the Poisson or Log-normal regression model}
#'   \item{'precision_R_effects': Precision parameters for random degree heterogeneity effects, specific to the Poisson or Log-normal regression model:
#'     \itemize{
#'       \item{'NDH': does not apply, can leave as missing}
#'       \item{'RS': a numeric value specifying the precision parameter of the normal distribution of the random sociality effect, if missing will generate from a gamma(shape = 1, rate = 1)}
#'       \item{'RSR': a numeric matrix specifying the precision matrix of the multivariate normal distribution of the random sender and receiver effects, if missing will generate from a Wishart(df = 3, Sigma = \eqn{I_2})}
#'     } 
#'    }
#'   \item{'precision_weights': A positive, non-zero, numeric representing the precision (on the log scale) of the log-normal weight distribution. Only relevant for \code{family = 'lognormal'}}
#'   }
#' @param q_prob A numeric in \[0,1\] representing the probability of a non-edge to be converted to noise edge (default is 0.0).
#' @param mean_noise_weights A numeric representing the mean of the noise weight distribution. Only relevant for \code{family \%in\% c('lognormal', 'poisson')} and \code{q_prob>0.0}. For family = 'poisson' value has to be > 0.0, for family = "lognormal" the mean is on the log scale.
#' @param precision_noise_weights A positive, non-zero, numeric representing the precision of the log-normal noise weight distribution. Only relevant for \code{family = 'lognormal'} and \code{q_prob>0.0}.
#' @param remove_isolates A logical; if \code{TRUE} then isolates from the network are removed (default is \code{TRUE}).
#' @return A list containing the following components:
#' \item{\code{A}}{ A sparse adjacency matrix of class 'dgCMatrix' representing the simulated network.}
#' \item{\code{Z}}{ A numeric \eqn{N \times K} cluster assignment matrix with rows representing the cluster an actor belongs to (i.e. indicated by a value of 1.0).}
#' \item{\code{U}}{ A numeric \eqn{N \times D} matrix with rows representing an actor's position in a \eqn{D}-dimensional social space. }
#' \item{\code{RE}}{ A numeric \eqn{N \times 1} matrix representing the actor specific random sociality effect (i.e., s) OR a \eqn{N \times 2} matrix representing the actor specific random sender and receiver effects (i.e., s and r, respectively).}
#' \item{\code{precision_R_effects}}{ The specific precision_R_effects used to simulate \code{RE}.}
#' \item{\code{model}}{ A character string representing the specific \code{model} used to simulate the network.}
#' @examples
#' \donttest{
#' mus <- matrix(c(-1,-1,1,-1,1,1), 
#'               nrow = 3,
#'               ncol = 2, 
#'               byrow = TRUE)
#' omegas <- array(c(diag(rep(7,2)),
#'                   diag(rep(7,2)), 
#'                   diag(rep(7,2))), 
#'                   dim = c(2,2,3))
#' p <- rep(1/3, 3)
#' beta0 <- 1.0
#' JANE::sim_A(N = 100L, 
#'             model = "NDH",
#'             mus = mus, 
#'             omegas = omegas, 
#'             p = p, 
#'             params_LR = list(beta0 = beta0),
#'             remove_isolates = TRUE)
#'}
#' @export

sim_A <- function(N,
                  mus, 
                  omegas, 
                  p, 
                  model = "NDH",
                  family = "bernoulli",
                  params_LR,
                  params_weights = NULL,
                  q_prob = 0.0,
                  mean_noise_weights,
                  precision_noise_weights,
                  remove_isolates = TRUE){
  
  if(!model %in% c("NDH", "RS", "RSR")){
    stop("Model needs to be one of the following: NDH, RS, or RSR")
  }
  
  if(!family %in% c("bernoulli", "lognormal", "poisson")){
    stop("family needs to be one of the following: 'bernoulli', 'lognormal', 'poisson'")
  }
  
  if(missing(params_LR) || is.null(params_LR$beta0)){
    stop("Please supply a named list for params_LR At minimum params_LR$beta0 needs to be specified")
  }
  
  if(!(0.0 <= q_prob & q_prob <= 1.0)){
    stop("Please supply a proportion in [0,1] for q_prob")
  }
  
  if( ( is.null(params_weights) & (family != "bernoulli"))  || ( (family != "bernoulli") && is.null(params_weights$beta0) ) ){
    stop("Please supply a named list for params_weights. At minimum, for family = 'poisson' a params_weights$beta0 needs to be specified and for family = 'lognromal' params_weights$beta0 and params_weights$precision_weights needs to be suppied")
  }
  
  if(missing(mean_noise_weights) & (family != "bernoulli") & (q_prob > 0.0)){
    stop("Please supply a numeric for mean_noise_weights")
  }
  
  if(missing(precision_noise_weights) & (family == "lognormal") & (q_prob > 0.0)){
    stop("Please supply a numeric >0 for precision_noise_weights")
  }
  
  if(is.null(params_weights$precision_weights) & (family == "lognormal")){
    stop("Please supply a numeric >0 for params_weights$precision_weights")
  }
  
  if( (!is.null(params_weights$precision_weights) & (family == "lognormal")) && params_weights$precision_weights <= 0.0){
    stop("Please supply a numeric >0 for params_weights$precision_weights")
  }
  
  if( (!missing(mean_noise_weights) & (family == "poisson")) && mean_noise_weights <= 0.0){
    stop("Please supply a numeric >0 for mean_noise_weights")
  }
  
  if( (!missing(precision_noise_weights) & (family == "lognormal")) && precision_noise_weights <= 0.0){
    stop("Please supply a numeric >0 for precision_noise_weights")
  }
  
  # draw Z
  Z <-  t(stats::rmultinom(n = N, size = 1, prob = p))
  
  # draw U
  U <- matrix(stats::rnorm(N*ncol(mus)), nrow = N, ncol = ncol(mus))
  n_k <- colSums(Z)
  for (k in 1:length(n_k)){
    U[Z[, k] == 1,] <- (U[Z[, k] == 1, ] %*% t(solve(chol(omegas[,,k])))) + ( rep(1, n_k[k]) %*% t(mus[k,]) )
  }
  
  if(model == "RS"){
    
    if(!is.null(params_LR$precision_R_effects)){
      if(!(length(params_LR$precision_R_effects) == 1 && params_LR$precision_R_effects > 0)){
        stop("For RS Model, with respect to params_LR, please supply a positive scalar precision value for omega_{1s}^2")
      }
    } else {
      params_LR$precision_R_effects <- stats::rgamma(n = 1, shape = 1, scale = 1)
    }
    
    # draw s
    s <- stats::rnorm(n = N, mean = 0, sd = sqrt(1/params_LR$precision_R_effects))
    
    # draw A
    A <- draw_A_RS_c(U = U, beta0 = params_LR$beta0, s = s)
    
    RE <- matrix(s, nrow = N, ncol = 1)
    colnames(RE) <- "s"
    
  } else if(model == "RSR"){
    
    if(!is.null(params_LR$precision_R_effects)){
      if(!(length(params_LR$precision_R_effects) == 4 && all(dim(params_LR$precision_R_effects) == c(2,2)) && all(eigen(params_LR$precision_R_effects)$values > 0))) {
        stop("For RSR Model, with respect to params_LR, please supply a 2x2 p.d. precision matrix")
      }
    } else {
      params_LR$precision_R_effects <- stats::rWishart(n = 1, df  = 2 + 1, Sigma = diag(2))[,,1]
    }
    
    # draw s and r
    RE <- matrix(stats::rnorm(n = N*2), ncol = 2) %*% t(solve(chol(params_LR$precision_R_effects)))
    
    # draw A
    A <- draw_A_RSR_c(U = U, beta0 = params_LR$beta0, s = RE[,1], r = RE[,2])
    colnames(RE) <- c("s", "r")
    
  } else {
    
    # draw A
    A <- draw_A_NDH_c(U = U, beta0 = params_LR$beta0)
    RE <- NULL
    params_LR$precision_R_effects <- NULL
    
  }
  
  if(remove_isolates){
    
    isolates <- which(rowSums(A)==0 & colSums(A)==0)
    
    if(length(isolates) == 0){
      
      message("No isolates to remove")
      
    } else {
      
      message(paste0(length(isolates), " isolate(s) removed"))
      A <- A[-isolates, -isolates]
      Z <- Z[-isolates, ]
      U <- U[-isolates, ]
      RE <- RE[-isolates, , drop = F]
      N <- nrow(A)
    }
    
  } 
  
  params_LR$RE <- RE
  
  # Extract edge indices
  W <- A
  edge_indices <- as.matrix(summary(W))
  
  if(model != "RSR"){
    
    edge_indices <- edge_indices[edge_indices[,"j"]>edge_indices[,"i"], ]
    
  }
  
  if(family == "bernoulli"){
    
    params_weights <- NULL
    
  } else {
  
    if(model == "RS"){
      
      if(!is.null(params_weights$precision_R_effects)){
        if(!(length(params_weights$precision_R_effects) == 1 && params_weights$precision_R_effects > 0)){
          stop("For RS Model, with respect to params_weights, please supply a positive scalar precision value for omega_{2s}^2")
        }
      } else {
        params_weights$precision_R_effects <- stats::rgamma(n = 1, shape = 1, scale = 1)
      }
      
      # draw s
      s <- stats::rnorm(n = N, mean = 0, sd = sqrt(1/params_weights$precision_R_effects))
      RE_W <- matrix(s, nrow = N, ncol = 1)
      colnames(RE_W) <- "s"
      
      compute_mean_edge_weight(edge_indices = edge_indices,
                               beta0 = params_weights$beta0,
                               RE = RE_W,
                               model = "RS")
      mean_edges <- edge_indices[, 3]
      
    } else if(model == "RSR"){
      
      if(!is.null(params_weights$precision_R_effects)){
        if(!(length(params_weights$precision_R_effects) == 4 && all(dim(params_weights$precision_R_effects) == c(2,2)) && all(eigen(params_weights$precision_R_effects)$values > 0))) {
          stop("For RSR Model, with respect to params_weights, please supply a 2x2 p.d. precision matrix")
        }
      } else {
        params_weights$precision_R_effects<- stats::rWishart(n = 1, df  = 2 + 1, Sigma = diag(2))[,,1]
      }
      
      RE_W <- matrix(stats::rnorm(n = N*2), ncol = 2) %*% t(solve(chol(params_weights$precision_R_effects)))
      colnames(RE_W) <- c("s", "r")
      
      compute_mean_edge_weight(edge_indices = edge_indices,
                               beta0 = params_weights$beta0,
                               RE = RE_W,
                               model = "RSR")
      mean_edges <- edge_indices[, 3]
      
    } else {
      
      RE_W <- NULL
      params_weights$precision_R_effects <- NULL
      mean_edges <- params_weights$beta0
      
    }
    
    params_weights$RE <- RE_W
    
    if(family == "lognormal"){
      
      true_weights <- unname(exp(mean_edges + stats::rnorm(n = nrow(edge_indices), mean = 0.0, sd = sqrt(1/params_weights$precision_weights))))
      edge_indices[, 3] <- true_weights
      W[as.matrix(edge_indices[, c("i", "j")])] <- true_weights
      
      if(model != "RSR"){
        
        W[as.matrix(edge_indices[, c("j", "i")])] <- true_weights
        
      }
      
    }
    
    if(family == "poisson"){
      
      true_weights <- extraDistr::rtpois(n = nrow(edge_indices), lambda = exp(mean_edges), a = 0.0)
      edge_indices[, 3] <- true_weights
      W[as.matrix(edge_indices[, c("i", "j")])] <- true_weights
      
      if(model != "RSR"){
        
        W[as.matrix(edge_indices[, c("j", "i")])] <- true_weights
        
      }
      
    }
    
  }
  
  if(q_prob > 0.0){
    
    # non-link indices
    temp_A <- (A - 1.0)*-1.0
    diag(temp_A) <- 0.0
    non_edge_indices <- as.matrix(summary(as(temp_A, "sparseMatrix")))
    
    if(model != "RSR"){
      
      non_edge_indices <- non_edge_indices[non_edge_indices[,"j"]>non_edge_indices[,"i"], ]
      
    }  
    
    # Draw cluster labels
    Z_W <-  apply(t(stats::rmultinom(n = nrow(non_edge_indices),
                                     size = 1.0, prob = c(1.0-q_prob, q_prob))),
                  1.0, which.max)
    
    noise_weights <- numeric(length = length(Z_W))
    
    if(family == "bernoulli"){
      
      noise_weights[which(Z_W==2)] <- 1.0
      mean_noise_weights <- NULL
      precision_noise_weights <- NULL
      
    } else if(family == "poisson"){
      
      noise_weights[which(Z_W==2)] <- extraDistr::rtpois(n = sum(Z_W==2), lambda = mean_noise_weights, a = 0.0)
      precision_noise_weights <- NULL
      
    } else {
      
      noise_weights[which(Z_W==2)] <- unname(exp(mean_noise_weights + stats::rnorm(n = sum(Z_W==2), mean = 0.0, sd = sqrt(1/precision_noise_weights))))
      
    }
    
    W[cbind(non_edge_indices[,1], non_edge_indices[,2])] <- noise_weights
    
    if(model != "RSR"){
      
      W[cbind(non_edge_indices[,2], non_edge_indices[,1])] <- noise_weights
      
    }
    
    non_edge_indices[, 3] <- noise_weights
    non_edge_indices <- cbind(non_edge_indices, Z_W)
    Z_W_out <- rbind(cbind(edge_indices, 1), 
                 non_edge_indices[non_edge_indices[, "Z_W"] == 2, ])
    colnames(Z_W_out) <- c("i", "j", "weight", "Z_W")
    Z_W_out <- Z_W_out[order(Z_W_out[, "j"],
                             Z_W_out[, "i"]), ]
                 
  } else {
    
    mean_noise_weights <- NULL
    precision_noise_weights <- NULL
    non_edge_indices <- NULL
    Z_W_out <- NULL
    
  }

  return(list(A = A,
              W = W,
              Z_U = Z,
              U = U,
              Z_W = Z_W_out,
              params_LR = params_LR,
              params_weights = params_weights,
              q_prob = q_prob,
              mean_noise_weights = mean_noise_weights,
              precision_noise_weights = precision_noise_weights,
              model = model))
  
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

#' @useDynLib JANE   
compute_mean_edge_weight <- function(edge_indices, beta0, RE, model) {
  invisible(.Call('_JANE_compute_mean_edge_weight', PACKAGE = 'JANE', edge_indices, beta0, RE, model))
}
