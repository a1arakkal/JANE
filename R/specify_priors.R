#' Specify prior hyperparameters for EM algorithm
#' @description A function that allows the user to specify the prior hyperparameters for the EM algorithm in a structure accepted by \code{\link[JANE]{JANE}}. 
#' @param D An integer specifying the dimension of the latent positions (default is 2).
#' @param K An integer specifying the total number of clusters (default is 2).
#' @param family A character string specifying the distribution of the edge weights.
#'  \itemize{
#'   \item{'bernoulli': for \strong{unweighted} networks; utilizes a Bernoulli distribution with a logit link (default)}
#'   \item{'lognormal': for \strong{weighted} networks with positive, non-zero, continuous edge weights; utilizes a log-normal distribution with an identity link}
#'   \item{'poisson': for \strong{weighted} networks with edge weights representing non-zero counts; utilizes a zero-truncated Poisson distribution with a log link}
#'   }
#' @param model A character string specifying the model:
#' @param model A character string specifying the model:
#'  \itemize{
#'   \item{'NDH': \strong{undirected} network with no degree heterogeneity (or connection strength heterogeneity if working with weighted network)}
#'   \item{'RS': \strong{undirected} network with degree heterogeneity (and connection strength heterogeneity if working with weighted network)}
#'   \item{'RSR': \strong{directed} network with degree heterogeneity (and connection strength heterogeneity if working with weighted network)}
#'   }
#' @param noise_weights A logical; if TRUE then a Hurdle model is used to account for noise weights, if FALSE simply utilizes the supplied network (converted to an unweighted binary network if a weighted network is supplied, i.e., (A > 0.0)*1.0) and fits a latent space cluster model (default is FALSE).
#' @param n_interior_knots An integer specifying the number of interior knots used in fitting a natural cubic spline for degree heterogeneity (and connection strength heterogeneity if working with weighted network) models (i.e., 'RS' and 'RSR' only; default is \code{NULL}).   
#' @param a A numeric vector of length \eqn{D} specifying the mean of the multivariate normal prior on \eqn{\mu_k} for \eqn{k = 1,\ldots,K}, where \eqn{\mu_k} represents the mean of the multivariate normal distribution for the latent positions of the \eqn{k^{th}} cluster.
#' @param b A positive numeric scalar specifying the scaling factor on the precision of the multivariate normal prior on \eqn{\mu_k} for \eqn{k = 1,\ldots,K}, where \eqn{\mu_k} represents the mean of the multivariate normal distribution for the latent positions of the \eqn{k^{th}} cluster.
#' @param c A numeric scalar \eqn{\ge} \eqn{D} specifying the degrees of freedom of the Wishart prior on \eqn{\Omega_k} for \eqn{k = 1,\ldots,K}, where \eqn{\Omega_k} represents the precision of the multivariate normal distribution for the latent positions of the \eqn{k^{th}} cluster.
#' @param G A numeric \eqn{D \times D} matrix specifying the inverse of the scale matrix of the Wishart prior on \eqn{\Omega_k} for \eqn{k = 1,\ldots,K}, where \eqn{\Omega_k} represents the precision of the multivariate normal distribution for the latent positions of the \eqn{k^{th}} cluster.
#' @param nu A positive numeric vector of length \eqn{K} specifying the concentration parameters of the Dirichlet prior on \eqn{p}, where \eqn{p} represents the mixture weights of the finite multivariate normal mixture distribution for the latent positions.
#' @param e A numeric vector of length \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the mean of the multivariate normal prior on \eqn{\beta_{LR}}, where \eqn{\beta_{LR}} represents the coefficients of the logistic regression model.
#' @param f A numeric p.s.d square matrix of dimension \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the precision of the multivariate normal prior on \eqn{\beta_{LR}}, where \eqn{\beta_{LR}} represents the coefficients of the logistic regression model.
#' @param h A positive numeric scalar specifying the first shape parameter for the Beta prior on \eqn{q}, where \eqn{q} is the proportion of non-edges in the "true" underlying network converted to noise edges. Only relevant when \code{noise_weights = TRUE}.
#' @param l A positive numeric scalar specifying the second shape parameter for the Beta prior on \eqn{q}, where \eqn{q} is the proportion of non-edges in the "true" underlying network converted to noise edges. Only relevant when \code{noise_weights = TRUE}.
#' @param e_2 A numeric vector of length \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the mean of the multivariate normal prior on \eqn{\beta_{GLM}}, where \eqn{\beta_{GLM}} represents the coefficients of the zero-truncated Poisson or log-normal GLM. Only relevant when \code{noise_weights = TRUE & family != 'bernoulli'}.
#' @param f_2 A numeric p.s.d square matrix of dimension \code{1 + (model =='RS')*(n_interior_knots + 1) +  (model =='RSR')*2*(n_interior_knots + 1)} specifying the precision of the multivariate normal prior on \eqn{\beta_{GLM}}, where \eqn{\beta_{GLM}} represents the coefficients of the zero-truncated Poisson or log-normal GLM. Only relevant when \code{noise_weights = TRUE & family != 'bernoulli'}.
#' @param m_1 A positive numeric scalar specifying the shape parameter for the Gamma prior on \eqn{\tau^2_{weights}}, where \eqn{\tau^2_{weights}} is the precision (on the log scale) of the log-normal weight distribution. Note, this value is scaled by 0.5, see 'Details'. Only relevant when \code{noise_weights = TRUE & family = 'lognormal'}.
#' @param o_1 A positive numeric scalar specifying the rate parameter for the Gamma prior on \eqn{\tau^2_{weights}}, where \eqn{\tau^2_{weights}} is the precision (on the log scale) of the log-normal weight distribution. Note, this value is scaled by 0.5, see 'Details'. Only relevant when \code{noise_weights = TRUE & family = 'lognormal'}.
#' @param m_2 A positive numeric scalar specifying the shape parameter for the Gamma prior on \eqn{\tau^2_{noise \ weights}}, where \eqn{\tau^2_{noise \ weights}} is the precision (on the log scale) of the log-normal noise weight distribution. Note, this value is scaled by 0.5, see 'Details'. Only relevant when \code{noise_weights = TRUE & family = 'lognormal'}.
#' @param o_2 A positive numeric scalar specifying the rate parameter for the Gamma prior on \eqn{\tau^2_{noise \ weights}}, where \eqn{\tau^2_{noise \ weights}} is the precision (on the log scale) of the log-normal noise weight distribution. Note, this value is scaled by 0.5, see 'Details'. Only relevant when \code{noise_weights = TRUE & family = 'lognormal'}.
#' @details
#' 
#' \strong{Prior on \eqn{\boldsymbol{\mu}_k} and \eqn{\boldsymbol{\Omega}_k}} (note: the same prior is used for \eqn{k = 1,\ldots,K}) :
#' 
#' \deqn{\boldsymbol{\Omega}_k \sim Wishart(c, \boldsymbol{G}^{-1})}
#' \deqn{\boldsymbol{\mu}_k | \boldsymbol{\Omega}_k \sim MVN(\boldsymbol{a}, (b\boldsymbol{\Omega}_k)^{-1})}
#' 
#' \strong{Prior on \eqn{\boldsymbol{p}}}:
#' 
#' For the current implementation we require that all elements of the \code{nu} vector be \eqn{\ge 1} to prevent against negative mixture weights for empty clusters.
#' \deqn{\boldsymbol{p} \sim Dirichlet(\nu_1 ,\ldots,\nu_K)}
#' 
#' \strong{Prior on \eqn{\boldsymbol{\beta}_{LR}}}:
#' \deqn{\boldsymbol{\beta}_{LR} \sim MVN(\boldsymbol{e}, \boldsymbol{F}^{-1})}
#'
#' \strong{Prior on \eqn{q}}:
#' \deqn{q \sim Beta(h, l)}
#' 
#' \strong{\emph{Zero-truncated Poisson}}
#' 
#' \strong{Prior on \eqn{\boldsymbol{\beta}_{GLM}}}:
#' \deqn{\boldsymbol{\beta}_{GLM} \sim MVN(\boldsymbol{e}_{2}, \boldsymbol{F}_{2}^{-1})}
#' 
#' \strong{\emph{Log-normal}}
#' 
#' \strong{Prior on \eqn{\tau^2_{weights}}}:
#' \deqn{\tau^2_{weights} \sim Gamma(\frac{m_1}{2}, \frac{o_1}{2})}
#' 
#' \strong{Prior on \eqn{\boldsymbol{\beta}_{GLM}}}:
#' \deqn{\boldsymbol{\beta}_{GLM}|\tau^2_{weights} \sim MVN(\boldsymbol{e}_{2}, (\tau^2_{weights}\boldsymbol{F}_{2})^{-1})}
#' 
#' \strong{Prior on \eqn{\tau^2_{noise \ weights}}}:
#' \deqn{\tau^2_{noise \ weights} \sim Gamma(\frac{m_2}{2}, \frac{o_2}{2})}
#' 
#' Unevaluated calls can be supplied as values for specific hyperparameters. This is particularly useful when running \code{JANE} for multiple combinations of \code{K} and \code{D}. See 'examples' section below for implementation examples.
#' 
#' @return A list of S3 \code{\link{class}} "\code{JANE.priors}" representing prior hyperparameters for the EM algorithm, in a structure accepted by \code{JANE}.
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
#' # Specify prior hyperparameters as unevaluated calls
#' n_interior_knots <- 5L
#' e <- rep(0.5, 1 + (n_interior_knots + 1))
#' f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
#' 
#' my_prior_hyperparameters <- specify_priors(model = "RS",
#'                                            n_interior_knots = n_interior_knots,
#'                                            a = quote(rep(1, D)),
#'                                            b = b,
#'                                            c = quote(D + 1),
#'                                            G = quote(10*diag(D)),
#'                                            nu = quote(rep(2, K)),
#'                                            e = e,
#'                                            f = f)
#'                                            
#' # # Run JANE on simulated data using supplied prior hyperparameters (NOT RUN)
#' # future::plan(future::multisession, workers = 5)
#' # res <- JANE::JANE(A = sim_data$A,
#' #                    D = 2:5,
#' #                    K = 2:10,
#' #                    initialization = "GNN",
#' #                    model = "RS",
#' #                    case_control = FALSE,
#' #                    DA_type = "none",
#' #                    control = list(priors = my_prior_hyperparameters))
#' # future::plan(future::sequential)
#'                 
#'                                                          
#' }
#' @export
specify_priors <- function(D = 2,
                           K = 2,
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
    required <- defined[!defined %in% c("D", "K", "noise_weights", "family", "n_interior_knots", "h", "l", "e_2", "f_2", "m_1", "o_1", "m_2", "o_2")]
  } else {
    if(family == "bernoulli"){
      required <- defined[!defined %in% c("D", "K", "noise_weights", "family", "n_interior_knots", "e_2", "f_2", "m_1", "o_1", "m_2", "o_2")]
    } else if(family == "poisson"){
      required <- defined[!defined %in% c("D", "K", "noise_weights", "family", "n_interior_knots", "m_1", "o_1", "m_2", "o_2")]
    } else {
      required <- defined[!defined %in% c("D", "K", "noise_weights", "family", "n_interior_knots")]
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
                          function(x){is.numeric(eval(parse(text = x))) | is.call(eval(parse(text = x)))})
  
  if(any(!check_numeric)){
    stop(paste0("Please supply numeric values in the correct structure for: ", paste0(names(check_numeric[!check_numeric]), collapse=", ")))
  }
  
  # add check that nu >= 1
  if (any(if(is.call(nu)){eval(nu) < 1} else {nu < 1})){
    stop("For the current implementation we require that all elements of the nu vector be >= 1 to prevent against negative mixture probabilities for empty clusters")
  }
  
  if (xor(is.call(a), is.call(G))){
    stop("If either 'a' or 'G' is specified as a call (i.e., depending on D), then both must be specified as calls")
  }
  
  if (is.call(a) & is.call(G)){
    if(length(grep("\\bD\\b", as.character(a))) == 0 | length(grep("\\bD\\b", as.character(G))) == 0){
      stop("If 'a' and 'G' are specified as a call then both must be functions of 'D'")
    }
  }
  
  if (is.call(nu)){
    if(length(grep("\\bK\\b", as.character(nu))) == 0){
      stop("If 'nu' is specified as a call then it must be a function of 'K'")
    }
  }
  
  if(is.call(a)){
    test_a <- eval(a)
    if(is.null(dim(test_a)) || dim(test_a)[1] != 1){
      a <- bquote(t(.(a)))
    }
  } else {
    a <- t(a)
  }
  
  if (!noise_weights){
    priors <- list(
      a = a,
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
        a = a,
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
        a = a,
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
        a = a,
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
  
  test_priors <- lapply(priors, eval, envir = environment())
  
  check_priors(priors = test_priors,
               D = D,
               K = K,
               n_interior_knots = n_interior_knots,
               model = model,
               family = family,
               noise_weights = noise_weights)
  
  return(structure(priors, class = "JANE.priors"))
}
