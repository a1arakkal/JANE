#' Simulate unweighted networks from latent space cluster models
#' @description Simulate an unweighted network from a \eqn{D}-dimensional latent space cluster model with \eqn{K} clusters and \eqn{N} actors. The \emph{squared} euclidean distance is used (i.e., \eqn{dist(U_i,U_j)^2}), where \eqn{U_i} and \eqn{U_j} are the respective actor's positions in an unobserved social space.
#' @param N An integer specifying the number of actors in the network.
#' @param mus A numeric \eqn{K \times D} matrix specifying the mean vectors of the multivariate normal distribution for the latent positions of the \eqn{K} clusters.
#' @param omegas A numeric \eqn{D \times D \times K} array specifying the precision matrices of the multivariate normal distribution for the latent positions of the \eqn{K} clusters.
#' @param p A numeric vector of length \eqn{K} specifying the mixture weights of the finite multivariate normal mixture distribution for the latent positions.
#' @param beta0 A numeric value specifying the intercept parameter for the logistic regression model.
#' @param model A character string to specify the model to simulate the network from:
#'  \itemize{
#'   \item{'NDH': generates an \strong{undirected} network with no degree heterogeneity}
#'   \item{'RS': generates an \strong{undirected} network with degree heterogeneity, specifically by including actor specific random sociality effects}
#'   \item{'RSR': generates a \strong{directed} network with degree heterogeneity, specifically by including actor specific random sender and receiver effects}
#'   }
#' @param precision_R_effects Precision parameters for random degree heterogeneity effects:
#'  \itemize{
#'   \item{'NDH': does not apply, can leave as missing}
#'   \item{'RS': a numeric value specifying the precision parameter of the normal distribution of the random sociality effect, if missing will generate from a gamma(shape = 1, rate = 1)}
#'   \item{'RSR': a numeric matrix specifying the precision matrix of the multivariate normal distribution of the random sender and receiver effects, if missing will generate from a Wishart(df = 3, Sigma = \eqn{I_2})}
#'   } 
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
#'             beta0 = beta0, 
#'             remove_isolates = TRUE)
#'}
#' @export

sim_A <- function(N, mus, omegas, p, beta0, 
                  model,
                  precision_R_effects,
                  remove_isolates = TRUE){
  
  if(!model %in% c("NDH", "RS", "RSR")){
    stop("Model needs to be one of the following: NDH, RS, or RSR")
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
