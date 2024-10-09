#' Specify prior hyperparameters for EM algorithm
#' @description A function that allows the user to specify the prior hyperparameters for the EM algorithm in a structure accepted by \code{JANE}. 
#' @param D An integer specifying the dimension of the latent positions.
#' @param K An integer specifying the total number of clusters.
#' @param model A character string specifying the model:
#'  \itemize{
#'   \item{'NDH': \strong{undirected} network with no degree heterogeneity}
#'   \item{'RS': \strong{undirected} network with degree heterogeneity}
#'   \item{'RSR': \strong{directed} network with degree heterogeneity}
#'   }
#' @param n_interior_knots An integer specifying the number of interior knots used in fitting a natural cubic spline for degree heterogeneity models (i.e., 'RS' and 'RSR' only; default is \code{NULL}).   
#' @param a A numeric vector of length D...
#' @param b A numeric value...
#' @param c A numeric value...
#' @param G A numeric \eqn{N \times D} matrix...
#' @param nu A numeric vector of length K...
#' @param e A numeric vector. Specifically a vector of length \code{1 + (model =="RS")*(n_interior_knots + 1) +  (model =="RSR")*2*(n_interior_knots + 1)}...
#' @param f A numeric square matrix. Specifically dimension \code{1 + (model =="RS")*(n_interior_knots + 1) +  (model =="RSR")*2*(n_interior_knots + 1)}...
#' @details
#' For the current implementation we require that all elements of the nu vector be >= 1 to prevent against negative mixture probabilities for empty clusters.
#' @return A list of prior hyperparameters for the EM algorithm generated from the input values in a structure accepted by \code{JANE}.
#' @export
specify_priors <- function(D, 
                           K,
                           model,
                           n_interior_knots = NULL,
                           a, # prior on mean of mus (vector of length D)
                           b, # prior on precision of mus (scalar)
                           c, # prior on df for omegas (scalar)
                           G, # prior on sigma for omegas (DxD matrix)
                           nu, # prior for dirichlet vector of length k
                           e, # prior mean on beta (vector of length 1 + (model =="RS")*(n_interior_knots + 1) +  (model =="RSR")*2*(n_interior_knots + 1))
                           f){ # prior precision on beta (square matrix of dim 1 + (model =="RS")*(n_interior_knots + 1) +  (model =="RSR")*2*(n_interior_knots + 1))
  
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
  
  # Stop if everything but model is not numeric
  check_numeric <- sapply(defined[names(defined) != "model"], is.numeric) 
  if(any(!check_numeric)){
    stop(paste0("Please supply numeric values in the correct structure for: ", paste0(names(check_numeric[!check_numeric]), collapse=", ")))
  }
  
  # add check that nu >= 1
  if (any(nu < 1)){
    stop("For the current implementation we require that all elements of the nu vector be >= 1 to prevent against negative mixture probabilities for empty clusters")
  }
  
  priors <- list(
    a = t(a),
    b = b,
    c = c,
    G = G,
    nu = nu,
    e = e,
    f = f
  )
  
  check_priors(priors = priors,
               D = D,
               K = K,
               n_interior_knots = n_interior_knots,
               model = model)
  
  return(priors)
}