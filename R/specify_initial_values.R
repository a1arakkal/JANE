#' Specify starting values for EM algorithm
#' @description A function that allows the user to specify starting values for the EM algorithm in a structure accepted by \code{JANE}. 
#' @param A A square matrix or sparse matrix of class 'dgCMatrix' representing the adjacency matrix of the unweighted network of interest.
#' @param D An integer specifying the dimension of the latent positions.
#' @param K An integer specifying the total number of clusters.
#' @param model A character string specifying the model:
#'  \itemize{
#'   \item{'NDH': \strong{undirected} network with no degree heterogeneity}
#'   \item{'RS': \strong{undirected} network with degree heterogeneity}
#'   \item{'RSR': \strong{directed} network with degree heterogeneity}
#'   }
#' @param n_interior_knots An integer specifying the number of interior knots used in fitting a natural cubic spline for degree heterogeneity models (i.e., 'RS' and 'RSR' only; default is \code{NULL}).   
#' @param mus A numeric \eqn{K \times D} matrix specifying the mean vectors of the multivariate normal distribution for the latent positions of the \eqn{K} clusters.
#' @param omegas A numeric \eqn{D \times D \times K} array specifying the precision matrices of the multivariate normal distribution for the latent positions of the \eqn{K} clusters.
#' @param p A numeric vector of length \eqn{K} specifying the mixture weights of the finite multivariate normal mixture distribution for the latent positions.
#' @param beta A numeric vector specifying the regression coefficients for the logistic regression model. Specifically, a vector of length \code{1 + (model =="RS")*(n_interior_knots + 1) +  (model =="RSR")*2*(n_interior_knots + 1)}.
#' @param U A numeric \eqn{N \times D} matrix with rows specifying an actor's position in a \eqn{D}-dimensional social space.
#' @param Z A numeric \eqn{N \times K} matrix with rows representing the conditional probability that an actor belongs to the cluster \eqn{K = k} for \eqn{k = 1,\ldots,K}.
#' @details
#' To match \code{\link[JANE]{JANE}}, this function will remove isolates from the adjacency matrix A and determine the total number of actors after excluding isolates. If this is not done, errors with respect to incorrect dimensions in the starting values will be generated when executing \code{\link[JANE]{JANE}}.
#' 
#' Similarly to match \code{\link[JANE]{JANE}}, if an unsymmetric adjacency matrix A is supplied for \code{model %in% c('NDH', 'RS')} the user will be asked if they would like to proceed with converting A to a symmetric matrix (i.e., \code{A <- 1.0 * ( (A + t(A)) > 0.0 )}).
#' @return A list of starting values for the EM algorithm generated from the input values in a structure accepted by \code{JANE}.
#' @examples
#' \donttest{
#' # Simulate network
#' mus <- matrix(c(-1,-1,1,-1,1,1), 
#'               nrow = 3,
#'               ncol = 2, 
#'               byrow = TRUE)
#' omegas <- array(c(diag(rep(7,2)),
#'                   diag(rep(7,2)), 
#'                   diag(rep(7,2))), 
#'                 dim = c(2,2,3))
#' p <- rep(1/3, 3)
#' beta0 <- -1
#' sim_data <- JANE::sim_A(N = 100L, 
#'                         model = "RSR",
#'                         mus = mus, 
#'                         omegas = omegas, 
#'                         p = p, 
#'                         params_LR = list(beta0 = beta0),
#'                         remove_isolates = TRUE)
#' 
#' # Specify starting values
#' D <- 3L
#' K <- 5L
#' N <- nrow(sim_data$A)
#' n_interior_knots <- 5L
#' 
#' U <- matrix(stats::rnorm(N*D), nrow = N, ncol = D)
#' omegas <- stats::rWishart(n = K, df = D+1, Sigma = diag(D))
#' mus <- matrix(stats::rnorm(K*D), nrow = K, ncol = D)
#' p <- extraDistr::rdirichlet(n = 1, rep(3,K))[1,]
#' Z <-  extraDistr::rdirichlet(n = N, alpha = rep(1, K))
#' beta <- stats::rnorm(n = 1 + 2*(1 + n_interior_knots))
#' 
#' my_starting_values <- JANE::specify_initial_values(A = sim_data$A,
#'                                                    D = D,
#'                                                    K = K,
#'                                                    model = "RSR",
#'                                                    n_interior_knots = n_interior_knots,
#'                                                    U = U,
#'                                                    omegas = omegas, 
#'                                                    mus = mus, 
#'                                                    p = p, 
#'                                                    Z = Z,
#'                                                    beta = beta)         
#' 
#' # Run JANE using my_starting_values (no need to specify D and K as function will 
#' # determine those values from my_starting_values)
#' res <- JANE::JANE(A = sim_data$A,
#'                   initialization = my_starting_values,
#'                   model = "RSR")
#' }
#' @export
specify_initial_values <- function(A,
                                   D,
                                   K,
                                   model,
                                   n_interior_knots = NULL,
                                   U,
                                   omegas, 
                                   mus, 
                                   p, 
                                   Z,
                                   beta){ 
  
  # Stop if any argument is missing
  defined <- ls()
  passed <- c(names(as.list(match.call())[-1]), "n_interior_knots")
  
  if (any(!defined %in% passed)) {
    stop(paste("Please supply values for argument(s): ", paste(setdiff(defined, passed), collapse=", ")))
  }
  
  # Stop if model not supplied or length > 1
  if(missing(model) || !(length(model) == 1 & is.character(model))){
    stop("Argument 'model' missing or not a character vector of length 1, please supply a model (i.e., 'NDH', 'RS', or 'RSR')")
  }
  
  # Check model 
  if(!model %in% c("NDH", "RS", "RSR")){
    stop("Model needs to be one of the following: 'NDH', 'RS', or 'RSR'")
  }
  
  # Stop if n_interior_knots is NULL for RS and RSR
  if((model %in% c("RS", "RSR")) & is.null(n_interior_knots)){
    stop("Model 'RS' or 'RSR' requires an integer value for n_interior_knots")
  }
  
  # Stop if D or K not numeric
  if(!(is.numeric(D) & is.numeric(K)) | !(length(D) == 1 & length(K) == 1)){
    stop("Please supply scalar integer values for D and K")
  }
  
  # Stop if A not supplied
  if(missing(A) || !(inherits(A, "matrix") | inherits(A, "dgCMatrix"))){
    stop("Argument 'A' missing or is not a matrix, please supply an adjacency matrix")
  }
  
  # Stop if A is not square
  if(nrow(A) != ncol(A)){
    stop("Please supply a square adjacency matrix")
  }

  # Check for class of A
  if(!"dgCMatrix" %in% class(A)){
    A <- methods::as(A, "dgCMatrix")
  }
  
  # Stop if not unweighted network
  if(!all(A@x == 1.0)){
    stop("Please supply an unweighted network")
  }
  
  # Check for self loops 
  if(!all(diag(A) == 0.0)){
    stop("Self-loop(s) detected (i.e., diagonal element(s) of A are not 0)")
  }
  
  # If unsymmetric A provided for model = "NDH" or "RS" convert to symmetric A and warn
  if(!isSymmetric(A) & (model %in% c("NDH", "RS"))){
    input <- utils::menu(c("Yes", "No"), 
                         title = paste0("Unsymmetric A matrix supplied for model = ",
                                        model, ", do you want to convert A to a symmetric matrix?"))
    if(input == 1){
      A <- 1.0 * ( (A + t(A)) > 0.0 )
      message("Converting A to symmetric matrix")
    } else {
      stop("A needs to be symmetric for model = ", model)
    }
  }

  # if A is named either col or row use that as id, else give id as row number
  if (( length(unique(rownames(A))) == nrow(A) ) | ( length(unique(colnames(A))) == nrow(A) )){
    ids <- if(length(unique(rownames(A))) != nrow(A)){ colnames(A) } else { rownames(A) }
  } else {
    message(paste0("Rownames/colnames for A missing or not unique, generating node IDs 1:", nrow(A)))
    ids <- 1:nrow(A)
  }
  
  # Check for isolates and remove
  isolates <- which(rowSums(A)==0 & colSums(A)==0)
  
  if (length(isolates) > 0){
    A <- A[-isolates, -isolates]
    message(paste0(length(isolates), " isolate(s) removed. Specifically node(s): ", paste(ids[isolates], collapse = ", ")))
    ids <- ids[-isolates]
  }
  
  # Check if max K supplied >= N
  if(K >= nrow(A)){
    stop("Number of clusters is greater than or equal to number of actors, please supply smaller K(s)")
  }
  
  list_name <- list(U = U, 
                    omegas = omegas, 
                    mus = mus, 
                    p = p, 
                    prob_matrix = Z,
                    beta = beta)
  
  check_initial_values(list_name = list_name,
                       A = A,
                       K = K,
                       D = D,
                       n_interior_knots = n_interior_knots,
                       model = model)
  return(list_name)
  
}
