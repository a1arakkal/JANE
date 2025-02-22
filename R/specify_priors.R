#' Specify prior hyperparameters for EM algorithm
#' @description A function that allows the user to specify the prior hyperparameters for the EM algorithm in a structure accepted by \code{JANE}. 
#' @param D An integer specifying the dimension of the latent positions.
#' @param K An integer specifying the total number of clusters.
#' @param family A character string specifying the model:
#' @param model A character string specifying the model:
#'  \itemize{
#'   \item{'NDH': \strong{undirected} network with no degree heterogeneity}
#'   \item{'RS': \strong{undirected} network with degree heterogeneity}
#'   \item{'RSR': \strong{directed} network with degree heterogeneity}
#'   }
#' @param noise_weights A character string specifying the model:
#' @param n_interior_knots An integer specifying the number of interior knots used in fitting a natural cubic spline for degree heterogeneity models (i.e., 'RS' and 'RSR' only; default is \code{NULL}).   
#' @param a A numeric vector of length \eqn{D} specifying the mean of the multivariate normal prior on \eqn{\mu_k} for \eqn{k = 1,\ldots,K}, where \eqn{\mu_k} represents the mean of the multivariate normal distribution for the latent positions of the \eqn{k^{th}} cluster.
#' @param b A numeric value specifying the scaling factor on the precision of the multivariate normal prior on \eqn{\mu_k} for \eqn{k = 1,\ldots,K}, where \eqn{\mu_k} represents the mean of the multivariate normal distribution for the latent positions of the \eqn{k^{th}} cluster.
#' @param c A numeric value specifying the degrees of freedom of the Wishart prior on \eqn{\Omega_k} for \eqn{k = 1,\ldots,K}, where \eqn{\Omega_k} represents the precision of the multivariate normal distribution for the latent positions of the \eqn{k^{th}} cluster.
#' @param G A numeric \eqn{D \times D} matrix specifying the inverse of the scale matrix of the Wishart prior on \eqn{\Omega_k} for \eqn{k = 1,\ldots,K}, where \eqn{\Omega_k} represents the precision of the multivariate normal distribution for the latent positions of the \eqn{k^{th}} cluster.
#' @param nu A numeric vector of length \eqn{K} specifying the concentration parameters of the Dirichlet prior on \eqn{p}, where \eqn{p} represents the mixture weights of the finite multivariate normal mixture distribution for the latent positions.
#' @param e A numeric vector of length \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the mean of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @param f A numeric square matrix of dimension \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the precision of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @param h A numeric vector of length \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the mean of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @param l A numeric square matrix of dimension \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the precision of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @param e_2 A numeric vector of length \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the mean of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @param f_2 A numeric square matrix of dimension \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the precision of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @param m_1 A numeric vector of length \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the mean of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @param o_1 A numeric square matrix of dimension \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the precision of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @param m_2 A numeric vector of length \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the mean of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @param o_2 A numeric square matrix of dimension \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the precision of the multivariate normal prior on \eqn{\beta}, where \eqn{\beta} represents the coefficients of the logistic regression model.
#' @details
#' 
#' \strong{Prior on \eqn{\mu_k} and \eqn{\Omega_k}} (note: the same prior is used for \eqn{k = 1,\ldots,K}) :
#' 
#' \eqn{\pi(\mu_k, \Omega_k) = \pi(\mu_k | \Omega_k) \pi(\Omega_k)}, thus
#' \deqn{\mu_k | \Omega_k \sim MVN(a, (b\Omega_k)^{-1})}
#' \deqn{\Omega_k \sim Wishart(c, G^{-1})}
#' 
#' \strong{Prior on \eqn{p}}:
#' 
#' For the current implementation we require that all elements of the nu vector be >= 1 to prevent against negative mixture weights for empty clusters.
#' \deqn{p \sim Dirichlet(\nu_1 ,\ldots,\nu_K)}
#' 
#' \strong{Prior on \eqn{\beta}}:
#' \deqn{\beta \sim MVN(e, f^{-1})}
#' 
#' @return A list of prior hyperparameters for the EM algorithm generated from the input values in a structure accepted by \code{JANE}.
#' 
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
#'                   dim = c(2,2,3))
#' p <- rep(1/3, 3)
#' beta0 <- 1.0
#' sim_data <- JANE::sim_A(N = 100L, 
#'                         model = "RS",
#'                         mus = mus, 
#'                         omegas = omegas, 
#'                         p = p, 
#'                         params_LR = list(beta0 = beta0), 
#'                         remove_isolates = TRUE)
#'                         
#'                         
#' # Specify prior hyperparameters
#' D <- 3L
#' K <- 5L
#' n_interior_knots <- 5L
#' 
#' a <- rep(1, D)
#' b <- 3
#' c <- 4
#' G <- 10*diag(D)
#' nu <- rep(2, K)
#' e <- rep(0.5, 1 + (n_interior_knots + 1))
#' f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
#' 
#' my_prior_hyperparameters <- specify_priors(D = D,
#'                                            K = K,
#'                                            model = "RS",
#'                                            n_interior_knots = n_interior_knots,
#'                                            a = a,
#'                                            b = b,
#'                                            c = c,
#'                                            G = G,
#'                                            nu = nu,
#'                                            e = e,
#'                                            f = f)
#'                                            
#' # Run JANE on simulated data using supplied prior hyperparameters
#' res <- JANE::JANE(A = sim_data$A,
#'                   D = D,
#'                   K = K,
#'                   initialization = "GNN",
#'                   model = "RS",
#'                   case_control = FALSE,
#'                   DA_type = "none",
#'                   control = list(priors = my_prior_hyperparameters))
#'                                                          
#' }
#' @export
specify_priors <- function(D, 
                           K,
                           model,
                           family = "bernoulli",
                           noise_weights = FALSE,
                           n_interior_knots = NULL,
                           a, 
                           b, 
                           c,
                           G, 
                           nu,
                           e, 
                           f,
                           h,
                           l,
                           e_2,
                           f_2,
                           m_1,
                           o_1,
                           m_2,
                           o_2){
  
  # Stop if any argument is missing
  defined <- ls()
  passed <- names(as.list(match.call())[-1])
  if(!noise_weights){
    required <- defined[!defined %in% c("noise_weights", "family", "n_interior_knots", "h", "l", "e_2", "f_2", "m_1", "o_1", "m_2", "o_2")]
  } else {
    if(family == "bernoulli"){
      required <- defined[!defined %in% c("noise_weights", "family", "n_interior_knots", "e_2", "f_2", "m_1", "o_1", "m_2", "o_2")]
    } else if(family == "poisson"){
      required <- defined[!defined %in% c("noise_weights", "family", "n_interior_knots", "m_1", "o_1", "m_2", "o_2")]
    } else {
      required <- defined[!defined %in% c("noise_weights", "family", "n_interior_knots")]
    }
  }
  
  if (any(!required %in% passed)) {
    stop(paste("Please supply values for argument(s): ", paste(setdiff(required, passed), collapse=", ")))
  }
  
  # Stop if model length > 1
  if(!(length(model) == 1 & is.character(model))){
    stop("Argument 'model'is not a character vector of length 1, please supply a model (i.e., 'NDH', 'RS', or 'RSR')")
  }

  # Stop if family length > 1
  if(!(length(family) == 1 & is.character(family))){
    stop("Argument 'family' is not a character vector of length 1, please supply a family (i.e., 'bernoulli', 'poisson', or 'lognormal')")
  }

  # Check family 
  if(!family %in% c("bernoulli", "poisson", "lognormal")){
    stop("family needs to be one of the following: 'bernoulli', 'poisson', or 'lognormal")
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
  
  # Stop if everything but model, family, noise_weights is not numeric
  check_numeric <- sapply(required[!(required %in% c("model", "family", "noise_weights"))],
                          function(x){is.numeric(eval(parse(text = x)))})
  
  if(any(!check_numeric)){
    stop(paste0("Please supply numeric values in the correct structure for: ", paste0(names(check_numeric[!check_numeric]), collapse=", ")))
  }
  
  # add check that nu >= 1
  if (any(nu < 1)){
    stop("For the current implementation we require that all elements of the nu vector be >= 1 to prevent against negative mixture probabilities for empty clusters")
  }
  
  
  if (!noise_weights){
    priors <- list(
      a = t(a),
      b = b,
      c = c,
      G = G,
      nu = nu,
      e = e,
      f = f
    )
  } else {
    if(family == "bernoulli"){
      priors <- list(
        a = t(a),
        b = b,
        c = c,
        G = G,
        nu = nu,
        e = e,
        f = f,
        h = h, 
        l = l
      )
    } else if(family == "poisson"){
      priors <- list(
        a = t(a),
        b = b,
        c = c,
        G = G,
        nu = nu,
        e = e,
        f = f,
        h = h, 
        l = l,
        e_2 = e_2,
        f_2 = f_2
      )
    } else {
      priors <- list(
        a = t(a),
        b = b,
        c = c,
        G = G,
        nu = nu,
        e = e,
        f = f,
        h = h, 
        l = l,
        e_2 = e_2,
        f_2 = f_2,
        m_1 = m_1,
        o_1 = o_1,
        m_2 = m_2,
        o_2 = o_2
      )
    }
  }
  
  check_priors(priors = priors,
               D = D,
               K = K,
               n_interior_knots = n_interior_knots,
               model = model,
               family = family,
               noise_weights = noise_weights)
  
  return(priors)
}