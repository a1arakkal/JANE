#' Fit JANE
#' @description Fit a latent space cluster model, with or without noise edges, using an EM algorithm.
#' @param A A square matrix or sparse matrix of class 'dgCMatrix' representing the adjacency matrix of the network of interest.
#' @param D Integer (scalar or vector) specifying the dimension of the latent space (default is 2).
#' @param K Integer (scalar or vector) specifying the number of clusters to consider (default is 2).
#' @param family A character string specifying the distribution of the edge weights.
#'  \itemize{
#'   \item{'bernoulli': Expects an \strong{unweighted} network; utilizes a Bernoulli distribution with a logit link (default)}
#'   \item{'lognormal': Expects a \strong{weighted} network with positive, non-zero, continuous edge weights; utilizes a log-normal distribution with an identity link}
#'   \item{'poisson': Expects a \strong{weighted} network with edge weights representing non-zero counts; utilizes a zero-truncated Poisson distribution with a log link}
#'   }
#' @param noise_weights A logical; if \code{TRUE} then a Hurdle model is used to account for noise weights, if \code{FALSE} simply utilizes the supplied network (converted to an unweighted binary network if a weighted network is supplied, i.e., (A > 0.0)*1.0) and fits a latent space cluster model (default is \code{FALSE}). 
#' @param guess_noise_weights Only applicable if \code{noise_weights = TRUE}. A numeric value specifying the best guess for the mean of the noise weight distribution for \code{family \%in\% c('lognormal', 'poisson')} (mean is on the log-scale for lognormal) \strong{OR} proportion (i.e. in (0,1)) of all edges that are noise edges for \code{family = 'bernoulli'}. If \code{NULL} (i.e., default) and \code{noise_weights = TRUE} then the 1st percentile of the non-zero weights will be used for \code{family \%in\% c('lognormal', 'poisson')} and 1% will be used for \code{family = 'bernoulli'}.
#' @param model A character string specifying the model to fit:
#'  \itemize{
#'   \item{'NDH': \strong{undirected} network with no degree heterogeneity}
#'   \item{'RS': \strong{undirected} network with degree heterogeneity}
#'   \item{'RSR': \strong{directed} network with degree heterogeneity}
#'   }
#' @param initialization A character string or a list to specify the initial values for the EM algorithm:
#'  \itemize{
#'   \item{'GNN': uses a type of graphical neural network approach to generate initial values (default)}
#'   \item{'random': uses random initial values}
#'   \item{A user supplied list of initial values. See \code{\link[JANE]{specify_initial_values}} on how to specify initial values}
#'   }
#' @param case_control A logical; if \code{TRUE} then uses a case-control approximation approach (default is \code{FALSE}).
#' @param DA_type (Experimental) A character string to specify the type of deterministic annealing approach to use
#'  \itemize{
#'   \item{'none': does not employ a deterministic annealing approach (default)}
#'   \item{'cooling': (Experimental) employes a traditional deterministic annealing approach where temperature decreases}
#'   \item{'heating': (Experimental) employes a deterministic anti-annealing approach where temperature increases}
#'   \item{'hybrid': (Experimental) employes a combination of the 'cooling' and 'heating' approach}
#'   }
#' @param seed (optional) An integer value to specify the seed for reproducibility.
#' @param control A list of control parameters. See 'Details'.
#' @return A list of S3 \code{\link{class}} "\code{JANE}" containing the following components:
#' \item{\code{input_params}}{ A list containing the input parameters for \code{IC_selection}, \code{case_control}, \code{DA_type}, \code{model}, \code{family}, and \code{noise_weights} used in the function call.}
#' \item{\code{A}}{ The square sparse adjacency matrix of class 'dgCMatrix' used in fitting the latent space cluster model. This matrix can be different than the input A matrix as isolates are removed.}
#' \item{\code{IC_out}}{ A matrix containing the relevant information criteria for all combinations of \code{K}, \code{D}, and \code{n_start} considered. The 'selected' column indicates the chosen optimal fit.}
#' \item{\code{all_convergence_ind}}{ A matrix containing the convergence information (i.e., 1 = converged, 0 = did not converge) and number of iterations for all combinations of \code{K}, \code{D}, \code{n_start}, and \code{beta_temperature} considered.} 
#' \item{\code{optimal_res}}{ A list containing the estimated parameters of interest based on the optimal fit selected. It is recommended to use \code{summary()} to extract the parameters of interest. See \code{\link[JANE]{summary.JANE}} for more details.}
#' \item{\code{optimal_starting}}{ A list containing the starting parameters used in the EM algorithm that resulted in the optimal fit selected. It is recommended to use \code{summary()} to extract the parameters of interest. See \code{\link[JANE]{summary.JANE}} for more details.}
#' @details
#' 
#' Isolates are removed from the adjacency matrix A. If an unsymmetric adjacency matrix A is supplied for \code{model %in% c('NDH', 'RS')} the user will be asked if they would like to proceed with converting A to a symmetric matrix (i.e., \code{A <- 1.0 * ( (A + t(A)) > 0.0 )}). Additionally, if a weighted network is supplied and \code{noise_weights = FALSE}, then the network will be converted to an unweighted binary network (i.e., (A > 0.0)*1.0) and a latent space cluster model is fit.
#'  
#' \strong{\code{control}:}
#' 
#' The \code{control} argument is a named list that the user can supply containing the following components:
#' \describe{
#' \item{\code{verbose}}{A logical; if \code{TRUE} causes additional information to be printed out about the progress of the EM algorithm (default is \code{FALSE}).}
#' \item{\code{max_its}}{An integer specifying the maximum number of iterations for the EM algorithm (default is \code{1e3}).}
#' \item{\code{min_its}}{An integer specifying the minimum number of iterations for the EM algorithm (default is \code{10}).}
#' \item{\code{priors}}{A list of prior hyperparameters (default is \code{NULL}). See \code{\link[JANE]{specify_priors}} on how to specify the hyperparameters.}
#' \item{\code{n_interior_knots}}{(only relevant for \code{model \%in\% c('RS', 'RSR')}) An integer specifying the number of interior knots used in fitting a natural cubic spline for degree heterogeneity models (default is \code{5}).}
#' \item{\code{termination_rule}}{A character string to specify the termination rule to determine the convergence of the EM algorithm: \itemize{
#'                                 \item{\code{'prob_mat'}: uses change in the absolute difference in \eqn{\hat{Z}^{U}} (i.e., the \eqn{N \times K} cluster membership probability matrix) between subsequent iterations (default)}
#'                                 \item{\code{'Q'}: uses change in the absolute difference in the objective function of the E-step evaluated using parameters from subsequent iterations}
#'                                 \item{\code{'ARI'}: comparing the classifications between subsequent iterations using adjusted Rand index}
#'                                 \item{\code{'NMI'}: comparing the classifications between subsequent iterations using normalized mutual information}
#'                                 \item{\code{'CER'}: comparing the classifications between subsequent iterations using classification error rate}
#'                                 }}
#' \item{\code{tolerance}}{A numeric specifying the tolerance used for \code{termination_rule \%in\% c('Q', 'prob_mat')} (default is \code{1e-3}).}                             
#' \item{\code{tolerance_ARI}}{A numeric specifying the tolerance used for \code{termination_rule = 'ARI'} (default is \code{0.999}).}    
#' \item{\code{tolerance_NMI}}{A numeric specifying the tolerance used for \code{termination_rule = 'NMI'} (default is \code{0.999}).}    
#' \item{\code{tolerance_CER}}{A numeric specifying the tolerance used for \code{termination_rule = 'CER'} (default is \code{0.01}).}    
#' \item{\code{n_its_start_CA}}{An integer specifying what iteration to start computing the change in cumulative averages (note: the change in the cumulative average of \eqn{\hat{U}}, the latent position matrix, is not tracked when \code{termination_rule = 'Q'}) (default is \code{20}).}
#' \item{\code{tolerance_diff_CA}}{A numeric specifying the tolerance used for the change in cumulative average of \code{termination_rule} metric and \eqn{\hat{U}} (note: the change in the cumulative average of \eqn{\hat{U}} is not tracked when \code{termination_rule = 'Q'}) (default is \code{1e-3}).}
#' \item{\code{consecutive_diff_CA}}{An integer specifying the tolerance for the number of consecutive instances where the change in cumulative average is less than \code{tolerance_diff_CA} (default is \code{5}).}    
#' \item{\code{quantile_diff}}{A numeric in \code{[0,1]} specifying the quantile used in computing the change in the absolute difference of \eqn{\hat{Z}^{U}} and \eqn{\hat{U}} between subsequent iterations (default is \code{1}, i.e., max).}
#' \item{\code{beta_temp_schedule}}{(Experimental) A numeric vector specifying the temperature schedule for deterministic annealing (default is \code{1}, i.e., deterministic annealing not utilized).}       
#' \item{\code{n_control}}{An integer specifying the fixed number of controls (i.e., non-links) sampled for each actor; only relevant when \code{case_control = TRUE} (default is \code{100} when \code{case_control = TRUE} and \code{NULL} when \code{case_control = FALSE}).}  
#' \item{\code{n_start}}{An integer specifying the maximum number of starts for the EM algorithm (default is \code{5}).}    
#' \item{\code{max_retry}}{An integer specifying the maximum number of re-attempts if starting values cause issues with EM algorithm (default is \code{5}).}    
#' \item{\code{IC_selection}}{A character string to specify the information criteria used to select the optimal fit based on the combinations of \code{K}, \code{D}, and \code{n_start} considered: \itemize{
#'                                 \item{\code{'BIC_model'}: BIC computed from logistic regression or Hurdle model component}
#'                                 \item{\code{'BIC_mbc'}: BIC computed from model based clustering component}
#'                                 \item{\code{'ICL_mbc'}: ICL computed from model based clustering component}
#'                                 \item{\code{'Total_BIC'}: sum of \code{'BIC_model'} and \code{'BIC_mbc'}}
#'                                 \item{\code{'Total_ICL'}: sum of \code{'BIC_model'} and \code{'ICL_mbc'} (default)}
#'                                 }}
#' \item{\code{sd_random_U_GNN}}{(only relevant when \code{initialization = 'GNN'}) A positive numeric value specifying the standard deviation for the random draws from a normal distribution to initialize \eqn{U} (default is \code{1}).} 
#' \item{\code{max_retry_GNN}}{(only relevant when \code{initialization = 'GNN'}) An integer specifying the maximum number of re-attempts for the \code{GNN} approach before switching to random starting values (default is \code{10}).} 
#' \item{\code{n_its_GNN}}{(only relevant when \code{initialization = 'GNN'}) An integer specifying the maximum number of iterations for the \code{GNN} approach (default is \code{10}).} 
#' \item{\code{downsampling_GNN}}{(only relevant when \code{initialization = 'GNN'}) A logical; if \code{TRUE} employs downsampling s.t. the number of links and non-links are balanced for the \code{GNN} approach (default is \code{TRUE}).}
#'}
#'
#' \strong{Running \code{JANE} in parallel:}
#' 
#' \code{JANE} integrates the \pkg{future} and \pkg{future.apply} packages to fit the various combinations of \code{K}, \code{D}, and \code{n_start} in parallel. The 'Examples' section below provides an example of how to run \code{JANE} in parallel. See \code{\link[future]{plan}} and \code{\link[future.apply]{future.apply}} for more details.
#' 
#' \strong{Choosing the number of clusters:}
#' 
#' \code{JANE} allows for the following model selection criteria to choose the number of clusters (smaller values are favored):
#' \itemize{
#'      \item{\code{'BIC_model'}: BIC computed from logistic regression or Hurdle model component}
#'      \item{\code{'BIC_mbc'}: BIC computed from model based clustering component}
#'      \item{\code{'ICL_mbc'}: ICL (Biernacki et al. (2000)) computed from model based clustering component}
#'      \item{\code{'Total_BIC'}: Sum of \code{'BIC_model'} and \code{'BIC_mbc'}, this is the model selection criterion proposed by Handcock et al. (2007)}
#'      \item{\code{'Total_ICL'}: (default) sum of \code{'BIC_model'} and \code{'ICL_mbc'}, this criterion is similar to \code{'Total_BIC'}, but uses ICL for the model based clustering component, which tends to favor more well-separated clusters.}
#'}
#' 
#' Based on simulation studies, Biernacki et al. (2000) recommends that when the interest in mixture modeling is cluster analysis, instead of density estimation, the \eqn{ICL_{mbc}} criterion should be favored over the \eqn{BIC_{mbc}} criterion, as the \eqn{BIC_{mbc}} criterion tends to overestimate the number of clusters. The \eqn{BIC_{mbc}} criterion is designed to choose the number of components in a mixture model rather than the number of clusters.
#' 
#' \emph{Warning}: It is not certain whether it is appropriate to use the model selection criterion above to select \code{D}.
#' 
#' @references
#' 
#' Biernacki, C., Celeux, G., Govaert, G., 2000. Assessing a mixture model for clustering with the integrated completed likelihood. IEEE Transactions on Pattern Analysis and Machine Intelligence 22, 719–725.
#'
#' Handcock, M.S., Raftery, A.E., Tantrum, J.M., 2007. Model-based clustering for social networks. Journal of the Royal Statistical Society Series A: Statistics in Society 170, 301–354.
#' 
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
#'                         model = "NDH",
#'                         mus = mus, 
#'                         omegas = omegas, 
#'                         p = p, 
#'                         params_LR = list(beta0 = beta0),
#'                         remove_isolates = TRUE)
#'                         
#' # Run JANE on simulated data
#' res <- JANE::JANE(A = sim_data$A,
#'                   D = 2L,
#'                   K = 3L,
#'                   initialization = "GNN", 
#'                   model = "NDH",
#'                   case_control = FALSE,
#'                   DA_type = "none")
#'
#' # Run JANE on simulated data - consider multiple D and K
#' res <- JANE::JANE(A = sim_data$A,
#'                   D = 2:5,
#'                   K = 2:10,
#'                   initialization = "GNN", 
#'                   model = "NDH",
#'                   case_control = FALSE,
#'                   DA_type = "none")
#'                   
#' # Run JANE on simulated data - parallel with 5 cores
#' future::plan(future::multisession, workers = 5)
#' res <- JANE::JANE(A = sim_data$A,
#'                   D = 2L,
#'                   K = 3L,
#'                   initialization = "GNN", 
#'                   model = "NDH",
#'                   case_control = FALSE,
#'                   DA_type = "none")
#' future::plan(future::sequential)
#' 
#' # Run JANE on simulated data - case/control approach with 20 controls sampled for each actor
#' res <- JANE::JANE(A = sim_data$A,
#'                   D = 2L,
#'                   K = 3L,
#'                   initialization = "GNN", 
#'                   model = "NDH",
#'                   case_control = TRUE,
#'                   DA_type = "none",
#'                   control = list(n_control = 20))
#'                    
#' # Reproducibility
#' res1 <- JANE::JANE(A = sim_data$A,
#'                    D = 2L,
#'                    K = 3L,
#'                    initialization = "GNN", 
#'                    seed = 1234,
#'                    model = "NDH",
#'                    case_control = FALSE,
#'                    DA_type = "none")
#' 
#' res2 <- JANE::JANE(A = sim_data$A,
#'                    D = 2L,
#'                    K = 3L,
#'                    initialization = "GNN", 
#'                    seed = 1234,
#'                    model = "NDH",
#'                    case_control = FALSE,
#'                    DA_type = "none")  
#' 
#' ## Check if results match
#' all.equal(res1, res2)    
#' 
#' # Another reproducibility example where the seed was not set. 
#' # It is possible to replicate the results using the starting values due to 
#' # the nature of EM algorithms
#' res3 <- JANE::JANE(A = sim_data$A,
#'                    D = 2L,
#'                    K = 3L,
#'                    initialization = "GNN", 
#'                    model = "NDH",
#'                    case_control = FALSE,
#'                    DA_type = "none")
#' ## Extract starting values                    
#' start_vals <- res3$optimal_start  
#' 
#' ## Run JANE using extracted starting values, no need to specify D and K 
#' ## below as function will determine those values from start_vals
#' res4 <- JANE::JANE(A = sim_data$A,
#'                    initialization = start_vals, 
#'                    model = "NDH",
#'                    case_control = FALSE,
#'                    DA_type = "none")
#'                    
#' ## Check if optimal_res are identical
#' all.equal(res3$optimal_res, res4$optimal_res)                   
#' }                            
#' @export
JANE <- function(A,
                 D = 2,
                 K = 2,
                 family = "bernoulli",
                 noise_weights = FALSE,
                 guess_noise_weights = NULL,
                 model,
                 initialization = "GNN", 
                 case_control = FALSE,
                 DA_type = "none", 
                 seed = NULL, 
                 control = list()){
  
  con <- list(
    verbose = FALSE,
    max_its = 1e3, 
    min_its = 10,  
    priors = NULL,  
    n_interior_knots = 5, 
    termination_rule = "prob_mat",
    tolerance = 1e-3, 
    tolerance_ARI = 0.999, 
    tolerance_NMI = 0.999, 
    tolerance_CER = 0.01,  
    n_its_start_CA = 20, 
    tolerance_diff_CA = 1e-3, 
    consecutive_diff_CA = 5,
    quantile_diff = 1,
    beta_temp_schedule = 1, 
    n_control = NULL, 
    n_start = 5,
    max_retry = 5, 
    IC_selection = "Total_ICL", 
    sd_random_U_GNN = 1, 
    max_retry_GNN = 10, 
    n_its_GNN = 10,
    downsampling_GNN = TRUE 
  )
  
  # Stop if A not supplied
  if(missing(A) || !(inherits(A, "matrix") | inherits(A, "dgCMatrix"))){
    stop("Argument 'A' missing or is not a matrix, please supply an adjacency matrix")
  }
  
  # Stop if A is not square
  if(nrow(A) != ncol(A)){
    stop("Please supply a square adjacency matrix")
  }
  
  # Stop if model not supplied
  if(missing(model) || !(length(model) == 1 & is.character(model))){
    stop("Argument 'model' missing or not a character vector of length 1, please supply a model (i.e., 'NDH', 'RS', or 'RSR')")
  }
  
  cl <- match.call()
  cl$A <- quote(A)
  cl$case_control <- NULL
  cl$seed <- NULL
  cl$DA_type <- NULL
  cl$family <- eval(family)
  cl$noise_weights <- eval(noise_weights)
  
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
  
  # Check if max K supplied >= N
  if(max(K) >= nrow(A)){
    stop("Number of clusters is greater than or equal to number of actors, please supply smaller K(s)")
  }
  
  # Check model 
  if(!model %in% c("NDH", "RS", "RSR")){
    stop("Model needs to be one of the following: 'NDH', 'RS', or 'RSR'")
  } else{
    cl$model <- eval(model)
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
  
  # Check family input
  if(!family %in% c("bernoulli", "lognormal", "poisson")){
    stop("family needs to be one of the following: 'bernoulli', 'lognormal', 'poisson'")
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
  if(family %in% c("lognormal", "poisson") & all(A@x == 1.0)){
    stop(paste0("Edges in the A matrix are unweighted with family = ", family, ". Please supply an A matrix with weighted edges for family %in% c('lognormal', 'poisson')"))
  }
  
  # Check if discrete data not supplied for family = "poisson"
  if(family == "poisson" & !all( (A@x %% 1.0) == 0.0 ) ){
    stop(paste0("Non-discrete edge weights detected with family = ", family, ". Please supply an A matrix with discrete edge weights for family == 'poisson'"))
  }
  
  # Check noise_weights input convert to unweighted network if noise_weights == FALSE & family %in% c("lognormal", "poisson")
  if(!noise_weights & family %in% c("lognormal", "poisson")){
    A <- 1.0 * ( A > 0.0 )
    family <- "bernoulli"
    message("noise_weights == FALSE & family %in% c('lognormal', 'poisson'), converting A to unweighted matrix and fitting latent space cluster model assuming no noise weights and family = 'bernoulli'")
  }
  
  # Check guess_noise_weights
  if(!noise_weights){
    
    guess_noise_weights <- NULL
    q_prob <- NULL
    
  } else {
    
    if(is.null(guess_noise_weights)){
      
      if(family != 'bernoulli'){
        
        prob_noise <- 0.01
        
        if(family == "poisson"){
          guess_noise_weights <- unname(stats::quantile(x = A@x, 
                                                        probs = prob_noise))
          density_A <- sum(A@x > guess_noise_weights)/(nrow(A) * (nrow(A)-1.0))
        } else {
          guess_noise_weights <- unname(stats::quantile(x = log(A@x), 
                                                        probs = prob_noise))
          density_A <- sum(log(A@x) > guess_noise_weights)/(nrow(A) * (nrow(A)-1.0))
        }
        
        q_prob <- (prob_noise*density_A)/((1.0-density_A)*(1.0-prob_noise))
        
      } else {
        
        guess_noise_weights <- 0.01
        density_A <- (length(A@x) * (1.0-guess_noise_weights))/(nrow(A) * (nrow(A)-1.0))
        q_prob <- (guess_noise_weights*density_A)/((1.0-density_A)*(1.0-guess_noise_weights))
        
      }
      
    } else {
      
      if(family != 'bernoulli'){
        
        if (guess_noise_weights <= 0 & family == "poisson"){
          stop("Please supply a postive non-zero numeric value for guess_noise_weights with family = 'poisson'.")
        }
       
        if (family == "poisson"){
          prob_noise <- mean(A@x <= guess_noise_weights)
          density_A <- sum(A@x > guess_noise_weights)/(nrow(A) * (nrow(A)-1.0))
          q_prob <- (prob_noise*density_A)/((1.0-density_A)*(1.0-prob_noise))
        } else {
          prob_noise <- mean(log(A@x) <= guess_noise_weights)
          density_A <- sum(log(A@x) > guess_noise_weights)/(nrow(A) * (nrow(A)-1.0))
          q_prob <- (prob_noise*density_A)/((1.0-density_A)*(1.0-prob_noise))
        }
        
      } else {
        if(!(0 < guess_noise_weights & guess_noise_weights < 1)){
          stop("Please supply a vaild proportion in (0,1) for guess_noise_weights with family ='bernoulli'.")
        }

        density_A <- (length(A@x) * (1.0-guess_noise_weights))/(nrow(A) * (nrow(A)-1.0))
        q_prob <- (guess_noise_weights*density_A)/((1.0-density_A)*(1.0-guess_noise_weights))
        
      }
      
    }
  }

  cl$guess_noise_weights <- eval(guess_noise_weights)
  cl$q_prob <- eval(q_prob)

  # Extract edge indices and weights and store as prob_matrix_W
  prob_matrix_W <- cbind(as.matrix(summary(A)), 1.0, 0.0)
  colnames(prob_matrix_W) <- c("i", "j", "weight", "hat_zij1", "hat_zij2")  
  
  # Only store upper triangular for RS and NDH as it is symmetric
  if(model != "RSR"){
    
    prob_matrix_W <- prob_matrix_W[prob_matrix_W[,"j"]>prob_matrix_W[,"i"], , drop = FALSE]
    
  }
  
  cl$prob_matrix_W <- eval(prob_matrix_W)
  
  # Convert A into unweighted network (weights already extracted and stored in prob_matrix_W)
  if(!all(A@x == 1.0)){
    A <- 1.0 * ( A > 0.0 )
  }
  
  # Check initialization
  if(!is.list(initialization) && (!initialization %in% c("random", "GNN"))){
    stop("Please provide one of the following for initialization: 'random', 'GNN', or a named list with the necessary starting paramters")
  } else {
    cl$initialization <- eval(initialization)
  }
  
  # Check DA_type
  if(!DA_type %in% c("none", "cooling", "heating", "hybrid")){
    stop("Please provide one of the following for DA_type: 'none', 'cooling', 'heating', or 'hybrid'")
  }
  
  # Update con n_control if case_control is T
  if(case_control == T){
    con$n_control <- 100
  } else {
    control$n_control <- NULL
  }
  
  # Update con beta_temp_schedule by DA_type argument
  if(DA_type == "cooling"){
    con$beta_temp_schedule <- seq(0.5, 1, length.out = 10)
  } else if (DA_type == "heating"){
    con$beta_temp_schedule <- seq(1.5, 1, length.out = 10)
  } else if (DA_type == "hybrid"){
    con$beta_temp_schedule <- c(seq(0.5, 1.5, length.out = 5), seq(1.4, 1, length.out = 5))
  } else {
    control$beta_temp_schedule <- 1
  }
  
  # Updated con n_start if user supplied starting values
  if (is.list(initialization)){
    control$n_start <- 1
  }
  
  # Check if names of elements in list match control_default
  nmsC <- names(con)
  namc <- names(control)
  
  # Check for duplicated params in control
  if(sum(duplicated(namc)) > 0){
    stop("Duplicated parameter(s) in control: ", namc[duplicated(namc)])
  }
  
  noNms <- namc[!namc %in% nmsC]
  if (length(noNms) != 0) {
    potential_matches <- nmsC[stringdist::amatch(x = noNms, table = nmsC, method = "jaccard", maxDist = Inf)]
    stop("Unknown parameter(s) in control: ", paste(noNms, collapse = ", "),
         ". Closest match: ", paste(unique(potential_matches), collapse = ", "))
  }
  
  con[namc] <- control
  
  # Check value of IC_selection
  if(!con[["IC_selection"]] %in% c("BIC_model", "BIC_mbc", "ICL_mbc", "Total_BIC", "Total_ICL")) {
    stop("Please provide one of the following for IC_selection: 'BIC_model', 'BIC_mbc', 'ICL_mbc', 'Total_BIC', or 'Total_ICL'")
  } 
  
  # Check value of n_control supplied and compute n_control
  if(!is.null(con[["n_control"]]) && !(con[["n_control"]]<= min(nrow(A)-rowSums(A), nrow(A)-colSums(A)) & con[["n_control"]]>0)){
    stop("Please supply a n_control value in (0, min(nrow(A)-rowSums(A), nrow(A)-colSums(A))]")
  } 
  
  # check value of quantile_diff
  if(!(con[["quantile_diff"]]>=0 & con[["quantile_diff"]]<=1)){
    stop("Please supply a quantile_diff in [0,1]")
  }
  
  # Check termination rule supplied 
  if(!con[["termination_rule"]] %in% c("ARI", "NMI", "CER", "prob_mat", "Q")){
    stop("Please provide one of the following termination rules: 'ARI', 'NMI', 'CER', 'prob_mat', or 'Q'")
  }
  
  # Check termination rule supplied when noise_weights = TRUE
  if(noise_weights & (con[["termination_rule"]] %in% c("ARI", "NMI", "CER", "Q"))){
    stop("For the current implementation, when noise_weights = TRUE the only available termination_rule is 'prob_mat'")
  }
  
  # Check n_its_start_CA
  if (con$n_its_start_CA < 1){
    stop("Please supply a n_its_start_CA value >= 1")
  }
  
  cl$control <- eval(con)
  
  # Check initialization list if supplied
  if(is.list(initialization)){
    
    K <- length(initialization$p)
    D <- ncol(initialization$U)
    cl$K <- K
    cl$D <- D
    
    check_initial_values(list_name = initialization,
                         A = A,
                         K = K,
                         D = D,
                         n_interior_knots = con$n_interior_knots,
                         model = model,
                         family = family,
                         noise_weights = noise_weights)
    
    if(noise_weights){
      initialization$q_prob <- q_prob
      initialization$guess_noise_weights <- guess_noise_weights
      cl$initialization <- eval(initialization)
    }
    
  }
  
  # Check prior if supplied
  if(!is.null(con$priors)){
    check_priors(priors = con$priors,
                 D = D,
                 K = K,
                 n_interior_knots = con$n_interior_knots,
                 model = model,
                 family = family,
                 noise_weights = noise_weights)
  }
  
  combinations_2run <- as.matrix(expand.grid(K = unique(K), D = unique(D), n_start = 1:con$n_start))
  combinations_2run <- combinations_2run[order(combinations_2run[, "K"],
                                               combinations_2run[, "D"], 
                                               combinations_2run[, "n_start"]), , drop = F]
  rownames(combinations_2run) <- NULL
  cl$combinations_2run <- combinations_2run
  
  old_handlers <- progressr::handlers() # returns current set of progression handlers
  
  progressr::handlers(progressr::handler_progress(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                                  complete = "=",   # Completion bar character
                                                  incomplete = "-", # Incomplete bar character
                                                  current = ">",    # Current bar character
                                                  clear = FALSE,    # If TRUE, clears the bar when finish
                                                  width = 100))
  
  if(nrow(combinations_2run) > 1){
  
    progressr::with_progress({
      p <- progressr::progressor(steps = nrow(combinations_2run))
      parallel_res <- future.apply::future_lapply(X = 1:nrow(combinations_2run), 
                                                  FUN = function(x){
                                                    out <- inner_parallel(x = x,
                                                                          call_def = cl,
                                                                          A = A)
                                                    p()
                                                    return(out)
                                                  },
                                                  future.globals = FALSE,
                                                  future.packages = "JANE",
                                                  future.seed = ifelse(is.null(seed), TRUE, seed))
    })
    
  } else {
    
    progressr::with_progress({
      parallel_res <- future.apply::future_lapply(X = 1:nrow(combinations_2run), 
                                                  FUN = function(x){
                                                    out <- inner_parallel(x = x,
                                                                          call_def = cl,
                                                                          A = A)
                                                    return(out)
                                                  },
                                                  future.globals = FALSE,
                                                  future.packages = "JANE",
                                                  future.seed = ifelse(is.null(seed), TRUE, seed))
    })
    
  }
  
  on.exit(progressr::handlers(old_handlers), add = TRUE)
  
  IC_out <- do.call("rbind", lapply(parallel_res, function(x){x$EM_results$IC}))
  IC_out <- cbind(combinations_2run, IC_out)

  # optimal based on min BIC --------------------------------------------------- 
  selected <- rep(0, nrow(IC_out))
  optimal_pos <- which(IC_out[, con$IC_selection] == min(IC_out[, con$IC_selection]))
  selected[optimal_pos] <- 1
  IC_out <- cbind(IC_out, selected)
  optimal_pos <- ifelse(length(optimal_pos) > 1, optimal_pos[1], optimal_pos)
  
  optimal_res <- parallel_res[[optimal_pos]]$EM_results
  if (length(optimal_res) > 1){
    optimal_res[["cluster_labels"]] <- apply(optimal_res$prob_matrix, 1, which.max)
    names(optimal_res[["cluster_labels"]]) <- ids
    rownames(optimal_res$prob_matrix) <- ids
    rownames(optimal_res$U) <- ids
  } else {
    optimal_res <- NULL
  }
  
  optimal_starting <- parallel_res[[optimal_pos]]$starting_params
  if(!is.null(optimal_starting) & is.list(optimal_starting)){
    optimal_starting[["cluster_labels"]] <- apply(optimal_starting$prob_matrix, 1, which.max)
    names(optimal_starting[["cluster_labels"]]) <- ids
    rownames(optimal_starting$prob_matrix) <- ids
    rownames(optimal_starting$U) <- ids
  } else {
    optimal_starting <- NULL
  }

  if(is.list(initialization)){
    IC_out[, "n_start"] <- NA
    IC_out[, "selected"] <- NA
  }
  
  # Keep info on convergence across all combinations_2run ----------------------
  all_convergence_ind <- do.call("rbind", lapply(1:length(parallel_res), function(x){
    
    if(!is.null(parallel_res[[x]]$EM_results$convergence_ind)){
      tmp <- rep(1, nrow(parallel_res[[x]]$EM_results$convergence_ind)) %x% combinations_2run[x, c("K", "D", "n_start"), drop = F]
      colnames(tmp) <-  c("K", "D", "n_start")
      cbind(tmp,
            parallel_res[[x]]$EM_results$convergence_ind)
    } else {
      convergence_ind <- matrix(0, nrow = length(con$beta_temp_schedule), ncol = 3)
      colnames(convergence_ind) <- c("beta_temperature", "convergence_ind", "n_iterations")
      convergence_ind[, "beta_temperature"] <- con$beta_temp_schedule
      tmp <- rep(1, nrow(convergence_ind)) %x% combinations_2run[x, c("K", "D", "n_start"), drop = F]
      colnames(tmp) <-  c("K", "D", "n_start")
      cbind(tmp, convergence_ind)
    }
    
  }))

  
  return(structure(list(optimal_res = optimal_res[sort(names(optimal_res))],
                        optimal_starting = optimal_starting[sort(names(optimal_starting))],
                        input_params =  list(IC_selection = con$IC_selection,
                                             case_control = case_control,
                                             DA_type = DA_type,
                                             model = model,
                                             family = family,
                                             noise_weights = noise_weights),
                        IC_out = IC_out,
                        all_convergence_ind = all_convergence_ind,
                        A = A), class = "JANE"))
  
}


EM_inner <- function(A,
                     D,
                     K,
                     model,
                     noise_weights,
                     prob_matrix_W,
                     family,
                     starting_params,
                     control,
                     ...){
  
  extra_args <- list(...)
  
  # Run initialize function
  current <- suppressWarnings(initialize_fun(A = A,
                                             list_name = starting_params, 
                                             family = family,
                                             noise_weights = noise_weights,
                                             prob_matrix_W = prob_matrix_W,
                                             model = model, 
                                             n_interior_knots = control$n_interior_knots,
                                             n_control = control$n_control, 
                                             priors = control$priors, 
                                             K = K,
                                             D = D))
  
  current$termination_rule <- control$termination_rule
  
  current$termination_metric <- array(NA, dim = c(control$max_its, 
                                                  ncol = 6 + (!(control$termination_rule %in% c("ARI", "NMI", "CER", "Q"))) - (control$termination_rule == "Q") + 2*(current$noise_weights)*(control$termination_rule != "Q"),
                                                  length(control$beta_temp_schedule)))
  
  colnames(current$termination_metric) <- if(ncol(current$termination_metric) == 9){
    
    c("ARI", "NMI", "CER", 
      paste0(control$termination_rule, ifelse(control$termination_rule == "prob_mat", 
                                              paste0("_*",control$quantile_diff),
                                              "")), 
      paste0(control$termination_rule, "_W", ifelse(control$termination_rule == "prob_mat", 
                                              paste0("_*",control$quantile_diff),
                                              "")), 
      "abs_diff_U", "abs_diff_MA_metric", "abs_diff_MA_metric_W", "abs_diff_MA_U")
  
  } else if(ncol(current$termination_metric) == 8){
    
    c("ARI", "NMI", "CER", paste0(control$termination_rule, "_W"), 
      "abs_diff_U", "abs_diff_MA_metric", "abs_diff_MA_metric_W", "abs_diff_MA_U")
    
  } else if(ncol(current$termination_metric) == 7){
    
    c("ARI", "NMI", "CER", paste0(control$termination_rule, ifelse(control$termination_rule == "prob_mat", 
                                                                   paste0("_*",control$quantile_diff),
                                                                   "")), 
      "abs_diff_U", "abs_diff_MA_metric", "abs_diff_MA_U")
    
  } else if(ncol(current$termination_metric) == 5){
    
    c("ARI", "NMI", "CER", control$termination_rule, "abs_diff_MA_metric")
    
  } else { 
    
    c("ARI", "NMI", "CER", "abs_diff_U", "abs_diff_MA", "abs_diff_MA_U")
    
  }
  
  col_term_metric <- grep(control$termination_rule, 
                          colnames(current$termination_metric))
  
  if (control$termination_rule != "Q"){
    col_term_U <- grep("abs_diff_U", 
                       colnames(current$termination_metric))
  }
  
  current$convergence_ind <- matrix(0, nrow = length(control$beta_temp_schedule), ncol = 3)
  colnames(current$convergence_ind) <- c("beta_temperature", "convergence_ind", "n_iterations")
  current$convergence_ind[, "beta_temperature"] <- control$beta_temp_schedule
  current$convergence_ind[, "n_iterations"] <- rep(control$max_its, time = length(control$beta_temp_schedule))
  
  if(nrow(extra_args$combinations_2run) == 1){
    p <-  progressr::progressor(steps = (control$max_its*length(control$beta_temp_schedule)) / 50) # updated every 50
    tick <- 0
  }
  
  # Start loop
  for (beta_temp in 1:length(control$beta_temp_schedule)){
    
    for(n_its in 1:control$max_its){
      
      # update prob_matrix
      current$fun_list$update_prob_matrix(prob_matrix = current$prob_matrix, 
                                          mus = current$mus, omegas = current$omegas, 
                                          p = current$p, U = current$U,
                                          temp_beta = as.double(control$beta_temp_schedule[beta_temp]))
      
      if(noise_weights){
        
        # update prob_matrix_W
        current$fun_list$update_prob_matrix_W(prob_matrix_W = current$prob_matrix_W,
                                              model = current$model,
                                              family = current$family,
                                              beta = current$beta,
                                              beta2 = if(current$family == "bernoulli"){matrix(0.0, 1, 1)}else{current$beta2},
                                              precision_weights = ifelse(current$family != "lognormal", 0.0, current$precision_weights),
                                              precision_noise_weights = ifelse(current$family != "lognormal", 0.0, current$precision_noise_weights),
                                              guess_noise_weights = current$guess_noise_weights,
                                              U = current$U,
                                              X = current$X,
                                              X2 = if(current$family == "bernoulli"){matrix(0.0, 1, 1)}else{current$X2},
                                              q = current$q_prob,
                                              temp_beta = as.double(control$beta_temp_schedule[beta_temp]))
        
        # update A
        A[current$prob_matrix_W[, 1:2, drop = FALSE]] <- current$prob_matrix_W[, 4]
        
        if(current$model != "RSR"){
          A[current$prob_matrix_W[, 2:1, drop = FALSE]] <- current$prob_matrix_W[, 4]
        }
        
        # update q_prob
        current$fun_list$update_q_prob(q_prob = current$q_prob,
                                       prob_matrix_W = current$prob_matrix_W,
                                       N = nrow(A),
                                       model = current$model,
                                       h = current$priors$h,
                                       l = current$priors$l)
        
        if(family != "bernoulli"){
          
          # update beta2
          update_beta2(beta2 = current$beta2,
                       prob_matrix_W = current$prob_matrix_W,
                       X2 = current$X2, 
                       model = current$model, 
                       family = current$family,
                       f_2 = if(current$model == "NDH"){matrix(current$priors$f_2, 1, 1)}else{current$priors$f_2}, 
                       e_2 = if(current$model == "NDH"){matrix(current$priors$e_2, 1, 1)}else{current$priors$e_2})
          
          
        }
        
        if (family == "lognormal"){
          
          # update precision_weights
          update_precision_weights(precision_weights = current$precision_weights,
                                   beta2 = current$beta2,
                                   prob_matrix_W = current$prob_matrix_W,
                                   X2 = current$X2, 
                                   model = current$model, 
                                   f_2 = if(current$model == "NDH"){matrix(current$priors$f_2, 1, 1)}else{current$priors$f_2}, 
                                   e_2 = if(current$model == "NDH"){matrix(current$priors$e_2, 1, 1)}else{current$priors$e_2},
                                   m_1 = current$priors$m_1,
                                   o_1 = current$priors$o_1)
          
          # update precision_noise_weights
          update_precision_noise_weights(precision_noise_weights = current$precision_noise_weights,
                                         prob_matrix_W = current$prob_matrix_W,
                                         guess_noise_weights = current$guess_noise_weights,
                                         m_2 = current$priors$m_2,
                                         o_2 = current$priors$o_2)
          
        }
        
      }
      
      # update U
      current$fun_list$update_U(U = current$U, 
                                A = A, 
                                n_control = control$n_control,
                                mus = current$mus, 
                                omegas = current$omegas, 
                                prob_matrix = current$prob_matrix, 
                                beta = current$beta,
                                X =  current$X,
                                model = current$model)
      
      # update p
      current$fun_list$update_p(prob_matrix = current$prob_matrix, p = current$p, nu = current$priors$nu)
      
      
      # update mus and omegas
      current$fun_list$update_mus_omegas(prob_matrix = current$prob_matrix,
                                         U = current$U, b = current$priors$b, a = current$priors$a,
                                         c = current$priors$c, G = current$priors$G,
                                         mus = current$mus, omegas = current$omegas)
      
      # update beta
      current$fun_list$update_beta(A = A, 
                                   n_control = control$n_control,
                                   U = current$U,
                                   beta = current$beta, 
                                   f = current$priors$f, 
                                   e = current$priors$e,
                                   X =  current$X,
                                   model = current$model)
      
      
      # check termination
      check_convergence <- terminate_EM(A = A,
                                        current = current,  
                                        termination_rule = control$termination_rule,
                                        tolerance = control$tolerance,
                                        tolerance_ARI = control$tolerance_ARI,
                                        tolerance_NMI = control$tolerance_NMI,
                                        tolerance_CER = control$tolerance_CER,
                                        quantile_diff = control$quantile_diff,
                                        n_control = control$n_control)
      
      
      # Additional check for convergence using cumulative average (CA)
      # Here we check if the CA of the termination metric is changing by iterations. If there
      # has been at least x control$n_its_start_CA number of consecutive iterations where 
      # the abs diff in CA is less than control$tolerance_diff_CA.
      
      if(n_its <= control$n_its_start_CA){
        
        if(!noise_weights){
          
          MA <- ifelse(is.infinite(check_convergence$metric[col_term_metric]), 0,
                       check_convergence$metric[col_term_metric])
          diff_MA <- Inf
          counter <- 0
          
        } else {
          
          MA <- ifelse(is.infinite(check_convergence$metric[col_term_metric[1]]), 0,
                       check_convergence$metric[col_term_metric[1]])
          diff_MA <- Inf
          counter <- 0
          
          MA_W <- ifelse(is.infinite(check_convergence$metric[col_term_metric[2]]), 0,
                       check_convergence$metric[col_term_metric[2]])
          diff_MA_W <- Inf
          counter_W <- 0
          
        }
        
        if(control$termination_rule != "Q"){
          MA_U <- ifelse(is.infinite(check_convergence$metric[col_term_U]), 0,
                         check_convergence$metric[col_term_U])
          diff_MA_U <- Inf
          counter_U <- 0
        }
        
      } else {
        
        n_its_temp <- n_its-control$n_its_start_CA+1
        
        if(!noise_weights){
          
          diff_MA_prev <- diff_MA
          diff_MA <- abs( ((MA - check_convergence$metric[col_term_metric]))/(n_its_temp) )
          counter <- (counter + 1) * (diff_MA<control$tolerance_diff_CA)*(diff_MA_prev<control$tolerance_diff_CA)
          MA <- (check_convergence$metric[col_term_metric] + (n_its_temp-1)*MA)/(n_its_temp)
          
        } else {
         
          diff_MA_prev <- diff_MA
          diff_MA <- abs( ((MA - check_convergence$metric[col_term_metric[1]]))/(n_its_temp) )
          counter <- (counter + 1) * (diff_MA<control$tolerance_diff_CA)*(diff_MA_prev<control$tolerance_diff_CA)
          MA <- (check_convergence$metric[col_term_metric[1]] + (n_its_temp-1)*MA)/(n_its_temp)
          
          diff_MA_prev_W <- diff_MA_W
          diff_MA_W <- abs( ((MA_W - check_convergence$metric[col_term_metric[2]]))/(n_its_temp) )
          counter_W <- (counter_W + 1) * (diff_MA_W<control$tolerance_diff_CA)*(diff_MA_prev_W<control$tolerance_diff_CA)
          MA_W <- (check_convergence$metric[col_term_metric[2]] + (n_its_temp-1)*MA_W)/(n_its_temp)
          
        }
        
        if(control$termination_rule != "Q"){
          diff_MA_U_prev <- diff_MA_U
          diff_MA_U <- abs( ((MA_U - check_convergence$metric[col_term_U]))/(n_its_temp) )
          counter_U <- (counter_U + 1) * (diff_MA_U<control$tolerance_diff_CA)*(diff_MA_U_prev<control$tolerance_diff_CA)
          MA_U <- (check_convergence$metric[col_term_U] + (n_its_temp-1)*MA_U)/(n_its_temp)
        }
        
      }
      
      if(control$termination_rule != "Q"){
        
        if(!noise_weights){
          
          current$termination_metric[n_its,,beta_temp] <- c(check_convergence$metric, diff_MA, diff_MA_U)
          
          if((check_convergence$terminate | (counter >= control$consecutive_diff_CA & counter_U >= control$consecutive_diff_CA)) & (n_its > control$min_its)){
            current$convergence_ind[beta_temp, "convergence_ind"] <- 1L
            current$convergence_ind[beta_temp, "n_iterations"] <- n_its
            break
          }
          
        } else {
          
          current$termination_metric[n_its,,beta_temp] <- c(check_convergence$metric, diff_MA, diff_MA_W, diff_MA_U)
          
          if((check_convergence$terminate | (counter >= control$consecutive_diff_CA & counter_W >= control$consecutive_diff_CA & counter_U >= control$consecutive_diff_CA)) & (n_its > control$min_its)){
            current$convergence_ind[beta_temp, "convergence_ind"] <- 1L
            current$convergence_ind[beta_temp, "n_iterations"] <- n_its
            break
          }
          
        }

      } else {
        
        current$termination_metric[n_its,,beta_temp] <- c(check_convergence$metric, diff_MA)
        
        if((check_convergence$terminate | (counter >= control$consecutive_diff_CA)) & (n_its > control$min_its)){
          current$convergence_ind[beta_temp, "convergence_ind"] <- 1L
          current$convergence_ind[beta_temp, "n_iterations"] <- n_its
          break
        }
        
      }
      
      if(nrow(extra_args$combinations_2run) == 1){
        tick <- tick + 1
        if(tick %% 50 == 0){
          p() # updated every 50
        }
      }
      
    }
    
  }
  
  current$termination_metric <- do.call("rbind", apply(current$termination_metric, 3, 
                                                       function(x){x[apply(x,1, function(x){!any(is.na(x) & !is.nan(x))}),]},
                                                       simplify = F))
  

  dont_return_names <- c("log_Q", "previous_prob_mat", "previous_prob_mat_W", "previous_U", "fun_list")
  out <- mget(x = names(current)[!(names(current) %in% dont_return_names)], envir = current)
  
  # Get IC info
  out$IC <- unlist(BICL(A = A, object = out))
  
  return(out[!(names(out) %in% c("X", "X2", "model", "family", "noise_weights"))])
  
}


inner_parallel <- function(x, call_def, A){
  
  combinations_2run_x <- call_def$combinations_2run[x, ]
  call_def$K <- combinations_2run_x["K"]
  call_def$D <- combinations_2run_x["D"]

  retry <- T
  retry_counter <- 0
  
  while(retry & retry_counter <= call_def$control$max_retry){
    
    if(!is.list(call_def$initialization)){
      
      call_def[[1]] <- as.symbol("initialize_starting_values")
      call_def$random_start <- ifelse(call_def$initialization == "GNN", F, T)
      call_def$starting_params <- eval(call_def)
      
    } else {
      
      call_def$starting_params <- call_def$initialization
      retry_counter <- Inf
      
    }
    
    run_fun <- tryCatch(
      {
        call_def[[1]] <- as.symbol("EM_inner")
        eval(call_def)
      },
      
      error = function(e) {
        if(call_def$control$verbose){
          if(!is.list(call_def$initialization)){
            message("Issues with starting values. Retrying with new starting values.\n")
          } 
        }
        NA
      },
      
      warning = function(w) {
        if(call_def$control$verbose){
          if(!is.list(call_def$initialization)){
            message("Issues with starting values. Retrying with new starting values.\n")
          } 
        }
        NA
      }
    ) 
    
    if(length(run_fun)>1){
      
      return(list(EM_results = run_fun,
                  starting_params = call_def$starting_params))
      
    } else {
      
      retry_counter <- retry_counter + 1
      
    }
    
  }
  
  if(retry){
    
    if(!is.list(call_def$initialization)){
      warning("Max retry (i.e., max_retry) attempts reached. Issues with starting values. Returning Inf values. If this occurs often consider using alternative initialization.")
    } else {
      warning("Issues with starting values. Returning Inf values. Consider using alternative initialization.")
    }

    
    return(list(EM_results = list(IC = c(BIC_model = Inf,
                                         BIC_mbc = Inf,
                                         ICL_mbc = Inf,
                                         Total_BIC = Inf,
                                         Total_ICL = Inf)),
                starting_params = Inf))
    
  }
  
}

#' @useDynLib JANE  
update_U <- function(U, A, mus, omegas, prob_matrix, beta, X, n_control, model) {
  invisible(.Call('_JANE_update_U', PACKAGE = 'JANE', U, A, mus, omegas, prob_matrix, beta, X, n_control, model))
}

#' @useDynLib JANE  
update_U_CC <- function(U, n_control, A, mus, omegas, prob_matrix, beta, X, model) {
  invisible(.Call('_JANE_update_U_CC', PACKAGE = 'JANE', U, n_control, A, mus, omegas, prob_matrix, beta, X, model))
}

#' @useDynLib JANE  
update_U_RE <- function(U, A, mus, omegas, prob_matrix, beta, X, model, n_control) {
  invisible(.Call('_JANE_update_U_RE', PACKAGE = 'JANE', U, A, mus, omegas, prob_matrix, beta, X, model, n_control))
}

#' @useDynLib JANE  
update_U_RE_CC <- function(U, n_control, A, mus, omegas, prob_matrix, beta, X, model) {
  invisible(.Call('_JANE_update_U_RE_CC', PACKAGE = 'JANE', U, n_control, A, mus, omegas, prob_matrix, beta, X, model))
}

#' @useDynLib JANE  
update_beta <- function(beta, A, U, f, e, X, n_control, model) {
  invisible(.Call('_JANE_update_beta', PACKAGE = 'JANE', beta, A, U, f, e, X, n_control, model))
}

#' @useDynLib JANE  
update_beta_CC <- function(beta, A, n_control, U, f, e, X, model) {
  invisible(.Call('_JANE_update_beta_CC', PACKAGE = 'JANE', beta, A, n_control, U, f, e, X, model))
}

#' @useDynLib JANE  
update_beta_RE <- function(beta, A, U, f, e, X, model, n_control) {
  invisible(.Call('_JANE_update_beta_RE', PACKAGE = 'JANE', beta, A, U, f, e, X, model, n_control))
}

#' @useDynLib JANE  
update_beta_RE_CC <- function(beta, A, n_control, U, f, e, X, model) {
  invisible(.Call('_JANE_update_beta_RE_CC', PACKAGE = 'JANE', beta, A, n_control, U, f, e, X, model))
}

#' @useDynLib JANE  
update_beta2 <- function(beta2, prob_matrix_W, f_2, e_2, X2, model, family) {
  invisible(.Call('_JANE_update_beta2', PACKAGE = 'JANE', beta2, prob_matrix_W, f_2, e_2, X2, model, family))
}

#' @useDynLib JANE  
update_mus_omegas <- function(prob_matrix, U, b, a, c, G, mus, omegas) {
  invisible(.Call('_JANE_update_mus_omegas', PACKAGE = 'JANE', prob_matrix, U, b, a, c, G, mus, omegas))
}

#' @useDynLib JANE  
update_p <- function(prob_matrix, p, nu) {
  invisible(.Call('_JANE_update_p', PACKAGE = 'JANE', prob_matrix, p, nu))
}

#' @useDynLib JANE  
update_prob_matrix_DA <- function(prob_matrix, mus, omegas, p, U, temp_beta) {
  invisible(.Call('_JANE_update_prob_matrix_DA', PACKAGE = 'JANE', prob_matrix, mus, omegas, p, U, temp_beta))
}

#' @useDynLib JANE  
update_prob_matrix_W_DA <- function(prob_matrix_W, model, family, beta, beta2, precision_weights, precision_noise_weights, guess_noise_weights, U, X, X2, q, temp_beta) {
  invisible(.Call('_JANE_update_prob_matrix_W_DA', PACKAGE = 'JANE', prob_matrix_W, model, family, beta, beta2, precision_weights, precision_noise_weights, guess_noise_weights, U, X, X2, q, temp_beta))
}

#' @useDynLib JANE  
update_q_prob <- function(q_prob, prob_matrix_W, model, N, h, l) {
  invisible(.Call('_JANE_update_q_prob', PACKAGE = 'JANE', q_prob, prob_matrix_W, model, N, h, l))
}

#' @useDynLib JANE  
trunc_poisson_density <- function(w, mean, log) {
  .Call('_JANE_trunc_poisson_density', PACKAGE = 'JANE', w, mean, log)
}

#' @useDynLib JANE  
lognormal_density <- function(w, precision, mean, log) {
  .Call('_JANE_lognormal_density', PACKAGE = 'JANE', w, precision, mean, log)
}

#' @useDynLib JANE  
update_precision_weights <- function(precision_weights, prob_matrix_W, model, beta2, X2, m_1, o_1, f_2, e_2) {
  invisible(.Call('_JANE_update_precision_weights', PACKAGE = 'JANE', precision_weights, prob_matrix_W, model, beta2, X2, m_1, o_1, f_2, e_2))
}

#' @useDynLib JANE  
update_precision_noise_weights <- function(precision_noise_weights, prob_matrix_W, guess_noise_weights, m_2, o_2) {
  invisible(.Call('_JANE_update_precision_noise_weights', PACKAGE = 'JANE', precision_noise_weights, prob_matrix_W, guess_noise_weights, m_2, o_2))
}
