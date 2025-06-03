#' Specify starting values for EM algorithm
#' @description A function that allows the user to specify starting values for the EM algorithm in a structure accepted by \code{\link[JANE]{JANE}}. 
#' @param A A square matrix or sparse matrix of class 'dgCMatrix' representing the adjacency matrix of the network of interest.
#' @param D An integer specifying the dimension of the latent positions.
#' @param K An integer specifying the total number of clusters.
#' @param family A character string specifying the distribution of the edge weights.
#'  \itemize{
#'   \item{'bernoulli': Expects an \strong{unweighted} network; utilizes a Bernoulli distribution with a logit link (default)}
#'   \item{'lognormal': Expects a \strong{weighted} network with positive, non-zero, continuous edge weights; utilizes a log-normal distribution with an identity link to model the \emph{non-noise weights} and log-normal distribution with mean (on log-scale) specified as \code{guess_noise_weights} to model the \emph{noise weights}}
#'   \item{'exp_lognormal': Expects a \strong{weighted} network with positive, non-zero, continuous edge weights; utilizes a log-normal distribution with an identity link to model the \emph{non-noise weights} and exponential distribution with mean specified as \code{guess_noise_weights} to model the \emph{noise weights}}
#'   \item{'poisson': Expects a \strong{weighted} network with edge weights representing non-zero counts; utilizes a zero-truncated Poisson distribution with a log link to model the \emph{non-noise weights} and zero-truncated Poisson distribution with mean specified as \code{guess_noise_weights} to model the \emph{noise weights}}
#'   }
#' @param model A character string specifying the model:
#'  \itemize{
#'   \item{'NDH': \strong{undirected} network with no degree heterogeneity}
#'   \item{'RS': \strong{undirected} network with degree heterogeneity}
#'   \item{'RSR': \strong{directed} network with degree heterogeneity}
#'   }
#' @param noise_weights A logical; if TRUE then a Hurdle model is used to account for noise weights, if FALSE simply utilizes the supplied network (converted to an unweighted binary network if a weighted network is supplied, i.e., (A > 0.0)*1.0) and fits a latent space cluster model (default is FALSE).
#' @param n_interior_knots An integer specifying the number of interior knots used in fitting a natural cubic spline for degree heterogeneity models (i.e., 'RS' and 'RSR' only; default is \code{NULL}).   
#' @param mus A numeric \eqn{K \times D} matrix specifying the mean vectors of the \eqn{K} \eqn{D}-variate normal distributions for the latent positions.
#' @param omegas A numeric \eqn{D \times D \times K} array specifying the precision matrices of the \eqn{K} \eqn{D}-variate normal distributions for the latent positions.
#' @param p A numeric vector of length \eqn{K} specifying the mixture weights of the finite multivariate normal mixture distribution for the latent positions.
#' @param beta A numeric vector specifying the regression coefficients for the logistic regression model. Specifically, a vector of length \cr \code{1 + (model =="RS")*(n_interior_knots + 1) +} \cr \code{(model =="RSR")*2*(n_interior_knots + 1)}.
#' @param beta2 A numeric vector specifying the regression coefficients for the zero-truncated Poisson or log-normal GLM. Specifically, a vector of length \cr \code{1 + (model =="RS")*(n_interior_knots + 1) +} \cr \code{(model =="RSR")*2*(n_interior_knots + 1)}. \cr Only relevant when \code{noise_weights = TRUE & family != 'bernoulli'}.
#' @param precision_weights A positive numeric scalar specifying the precision (on the log scale) of the log-normal weight distribution. Only relevant when \code{noise_weights = TRUE & family \%in\% c('lognormal', 'exp_lognormal')}.
#' @param precision_noise_weights A positive numeric scalar specifying the precision (on the log scale) of the log-normal noise weight distribution. Only relevant when \code{noise_weights = TRUE & family = 'lognormal'}.
#' @param U A numeric \eqn{N \times D} matrix with rows specifying an actor's position in a \eqn{D}-dimensional social space.
#' @param Z A numeric \eqn{N \times K} matrix with rows representing the conditional probability that an actor belongs to the cluster \eqn{K = k} for \eqn{k = 1,\ldots,K}.
#' @details
#' To match \code{\link[JANE]{JANE}}, this function will remove isolates from the adjacency matrix A and determine the total number of actors after excluding isolates. If this is not done, errors with respect to incorrect dimensions in the starting values will be generated when executing \code{\link[JANE]{JANE}}.
#' 
#' Similarly to match \code{\link[JANE]{JANE}}, if an unsymmetric adjacency matrix A is supplied for \code{model %in% c('NDH', 'RS')} the user will be asked if they would like to proceed with converting A to a symmetric matrix (i.e., \code{A <- 1.0 * ( (A + t(A)) > 0.0 )}). Additionally, if a weighted network is supplied and \code{noise_weights = FALSE}, then the network will be converted to an unweighted binary network (i.e., (A > 0.0)*1.0).
#' @return A list of starting values for the EM algorithm generated from the input values in a structure accepted by \code{\link[JANE]{JANE}}.
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
                                   family = "bernoulli",
                                   noise_weights = FALSE,
                                   n_interior_knots = NULL,
                                   U,
                                   omegas, 
                                   mus, 
                                   p, 
                                   Z,
                                   beta, 
                                   beta2,
                                   precision_weights,
                                   precision_noise_weights){ 
  
  # Stop if any argument is missing
  defined <- ls()
  passed <- names(as.list(match.call())[-1])
  
  if(!noise_weights){
    required <- defined[!defined %in% c("noise_weights", "family", "n_interior_knots", "beta2", "precision_weights", "precision_noise_weights")]
  } else {
    if(family == "bernoulli"){
      required <- defined[!defined %in% c("noise_weights", "family", "n_interior_knots", "beta2", "precision_weights", "precision_noise_weights")]
    } else if(family == "poisson"){
      required <- defined[!defined %in% c("noise_weights", "family", "n_interior_knots", "precision_weights", "precision_noise_weights")]
    } else if(family == "exp_lognormal"){
      required <- defined[!defined %in% c("noise_weights", "family", "n_interior_knots", "precision_noise_weights")]
    } else {
      required <- defined[!defined %in% c("noise_weights", "family", "n_interior_knots")]
    }
  }
  
  if (any(!required %in% passed)) {
    stop(paste("Please supply values for argument(s): ", paste(setdiff(required, passed), collapse=", ")))
  }

  # Stop if everything but model, family, noise_weights is not numeric
  check_numeric <- sapply(required[!(required %in% c("A", "model", "family", "noise_weights"))],
                          function(x){is.numeric(eval(parse(text = x)))})
  
  if(any(!check_numeric)){
    stop(paste0("Please supply numeric values in the correct structure for: ", paste0(names(check_numeric[!check_numeric]), collapse=", ")))
  }
  
  # Stop if model length > 1
  if(!(length(model) == 1 & is.character(model))){
    stop("Argument 'model' is not a character vector of length 1, please supply a model (i.e., 'NDH', 'RS', or 'RSR')")
  }
  
  # Stop if family length > 1
  if(!(length(family) == 1 & is.character(family))){
    stop("Argument 'family' is not a character vector of length 1, please supply a family (i.e., 'bernoulli', 'poisson', 'exp_lognormal', or 'lognormal')")
  }
  
  # Check family 
  if(!family %in% c("bernoulli", "poisson", "exp_lognormal", "lognormal")){
    stop("family needs to be one of the following: 'bernoulli', 'poisson', 'exp_lognormal', or 'lognormal")
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
  if(!(inherits(A, "matrix") | inherits(A, "dgCMatrix"))){
    stop("Argument 'A' is not a matrix, please supply an adjacency matrix")
  }
  
  # Stop if A is not square
  if(nrow(A) != ncol(A)){
    stop("Please supply a square adjacency matrix")
  }

  # Check for class of A
  if(!"dgCMatrix" %in% class(A)){
    A <- methods::as(A, "dgCMatrix")
  }
  
  # Check for self loops 
  if(!all(diag(A) == 0.0)){
    stop("Self-loop(s) detected (i.e., diagonal element(s) of A are not 0)")
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
  
  # Check if K supplied >= N
  if(K >= nrow(A)){
    stop("Number of clusters is greater than or equal to number of actors, please supply smaller K(s)")
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
  
  # Check if edges weights are >0
  if(!all(A@x > 0.0)){
    stop("Negative edge weights detected, edge weights need to be >0. This error can also be generated if a sparse adjacency matrix is supplied containing values of 0 for its 'non-zero' values (e.g., if a 'non-zero' element in the sparse matrix is replaced by a 0).")
  }
  
  # Check if weighted network supplied for family = "bernoulli" and ask user if they want to convert to an unweighted network
  if(family == "bernoulli" & !all(A@x == 1.0)){
    
    input <- utils::menu(c("Yes", "No"), 
                         title = paste0("Weighted edges found in A matrix for family = ",
                                        family, ", do you want to convert the edges in the A matrix to unweighted edges?"))
    if(input == 1){
      A <- 1.0 * ( A > 0.0 )
      message("Converting edges in A to unweighted edges")
    } else {
      stop("A needs to have unweighted edges for family = ", family)
    }
    
  }
  
  # Check if unweighted network supplied for family %in% c("lognormal", "poisson")
  if(family %in% c("lognormal", "exp_lognormal", "poisson") & all(A@x == 1.0)){
    stop(paste0("Edges in the A matrix are unweighted with family = ", family, ". Please supply an A matrix with weighted edges for family %in% c('lognormal', 'exp_lognormal', 'poisson')"))
  }
  
  # Check if discrete data not supplied for family = "poisson"
  if(family == "poisson" & !all( (A@x %% 1.0) == 0.0 ) ){
    stop(paste0("Non-discrete edge weights detected with family = ", family, ". Please supply an A matrix with discrete edge weights for family == 'poisson'"))
  }
  
  # Check noise_weights input convert to unweighted network if noise_weights == FALSE & family %in% c("lognormal", 'exp_lognormal', "poisson")
  if(!noise_weights & family %in% c("lognormal", "exp_lognormal", "poisson")){
    A <- 1.0 * ( A > 0.0 )
    family <- "bernoulli"
    message("noise_weights == FALSE & family %in% c('lognormal', 'exp_lognormal', 'poisson'), converting A to unweighted matrix and fitting latent space cluster model assuming no noise weights and family = 'bernoulli'")
  }
  
  list_name <- list(U = U, 
                    omegas = omegas, 
                    mus = mus, 
                    p = p, 
                    prob_matrix = Z,
                    beta = beta)
  
  if(noise_weights){
    
    if(family != "bernoulli"){
      list_name$beta2 <- beta2
    }
    
    if(family %in% c("lognormal", "exp_lognormal")){
      
      list_name$precision_weights <- precision_weights
      
      if(family == "lognormal"){
        list_name$precision_noise_weights <- precision_noise_weights
      }
      
    }
    
  }
  
  check_initial_values(list_name = list_name,
                       A = A,
                       K = K,
                       D = D,
                       n_interior_knots = n_interior_knots,
                       model = model,
                       family = family,
                       noise_weights = noise_weights)
  return(list_name)
  
}
