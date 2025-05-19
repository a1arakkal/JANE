#' Summarizing JANE fits
#' @description S3 summary method for object of class "\code{JANE}".
#' @param object An object of S3 \code{\link{class}} "\code{JANE}", a result of a call to \code{\link[JANE]{JANE}}. 
#' @param true_labels (optional) A numeric, character, or factor vector of known true cluster labels. Must have the same length as number of actors in the fitted network. Need to account for potential isolates removed (default is \code{NULL}). 
#' @param initial_values A logical; if \code{TRUE} then summarize fit using the starting parameters used in the EM algorithm (default is \code{FALSE}, i.e., the results after the EM algorithm is run are summarized).
#' @param ... Unused.
#' @return A list of S3 \code{\link{class}} "\code{summary.JANE}" containing the following components (Note: \eqn{N} is the number of actors in the network, \eqn{K} is the number of clusters, and \eqn{D} is the dimension of the latent space):
#' \item{\code{coefficients}}{ A list containing the estimated coefficients from the logistic regression model (i.e., 'beta_LR') and, if relevant, the estimated coefficients from the zero- truncated Poisson or log-normal GLM (i.e., 'beta_GLM').}
#' \item{\code{U}}{ A numeric \eqn{N \times D} matrix with rows containing an actor's estimated latent position in a \eqn{D}-dimensional social space.}
#' \item{\code{p}}{ A numeric vector of length \eqn{K} containing the estimated mixture weights of the finite multivariate normal mixture distribution for the latent positions.}
#' \item{\code{mus}}{ A numeric \eqn{K \times D} matrix containing the estimated mean vectors of the \eqn{K} \eqn{D}-variate normal distributions for the latent positions.}
#' \item{\code{omegas}}{ A numeric \eqn{D \times D \times K} array containing the estimated precision matrices of the \eqn{K} \eqn{D}-variate normal distributions for the latent positions.}
#' \item{\code{Z_U}}{ A numeric \eqn{N \times K} matrix with rows containing the estimated conditional probability that an actor belongs to the cluster \eqn{K = k} for \eqn{k = 1,\ldots,K}.}
#' \item{\code{uncertainty}}{ A numeric vector of length \eqn{N} containing the uncertainty of the \eqn{i^{th}} actor's classification, derived as 1 - \eqn{max_k \hat{Z}^{U}_{ik}}.}
#' \item{\code{cluster_labels}}{ A numeric vector of length \eqn{N} containing the cluster assignment of each actor based on a hard clustering rule of \eqn{\{h | \hat{Z}^{U}_{ih} = max_k \hat{Z}^{U}_{ik}\}}.}
#' \item{\code{Z_W}}{ A numeric \eqn{|E| \times 6} matrix, with \eqn{|E|} representing the total number of edges in the network (for undirected networks, only the upper diagonal edges are retained). The first two columns (i.e., 'i' and 'j') contains the specific indices of the edge between the \eqn{i^{th}} and \eqn{j^{th}} actors, the third column (i.e., 'weight') contains the specific edge weight, the fourth column (i.e., 'hat_zij1') contains the estimated conditional probability that the specific edge is a non-noise edge, the fifth column (i.e., 'hat_zij2') contains the estimated conditional probability that the specific edge is a noise edge, and the sixth column (i.e., 'noise_edge_cluster_labels') contains the noise-edge cluster assignment of each edge based on a hard clustering rule of \eqn{\{h | \hat{Z}^{W}_{eh} = max(\hat{Z}^{W}_{e1}, \hat{Z}^{W}_{e2})\}} for \eqn{e = 1,\ldots,|E|}, where \eqn{\hat{Z}^{W}_{e1}} and \eqn{\hat{Z}^{W}_{e2}} are the estimated conditional probabilities that the \eqn{e^{th}} edge is a non-noise and noise edge, respectively (labels defined as, 1: non-noise edge and 2: noise edge). Will be \code{NULL} if \code{noise_weights = FALSE} or \code{initial_values = TRUE}. }
#' \item{\code{q_prob}}{ A numeric scalar representing the estimated proportion of non-edges in the "true" unobserved network that were converted to noise edges.}
#' \item{\code{precision_weights}}{ A numeric scalar representing the estimated precision (on the log scale) of the log-normal weight distribution. Only relevant for \code{family = 'lognormal' & noise_weights = TRUE}. }
#' \item{\code{precision_noise_weights}}{ A numeric scalar representing the estimated precision (on the log scale) of the log-normal noise weight distribution. Only relevant for \code{family = 'lognormal' & noise_weights = TRUE}. }
#' \item{\code{IC}}{ Information criteria values of the optimal fit selected, including \itemize{
#'                                 \item{\code{'BIC_model'}: BIC computed from logistic regression or Hurdle model component}
#'                                 \item{\code{'BIC_mbc'}: BIC computed from model based clustering component}
#'                                 \item{\code{'ICL_mbc'}: ICL computed from model based clustering component}
#'                                 \item{\code{'Total_BIC'}: sum of \code{'BIC_model'} and \code{'BIC_mbc'}}
#'                                 \item{\code{'Total_ICL'}: sum of \code{'BIC_model'} and \code{'ICL_mbc'}}
#'                                 }}
#' \item{\code{input_params}}{ A list with the following components: \itemize{
#'                           \item{\code{model}: A character string containing the specific \code{model} used (i.e., 'NDH', 'RS', or 'RSR')}
#'                           \item{\code{family}: A character string containing the specific \code{family} used (i.e., 'bernoulli', 'poisson', or 'lognormal')}
#'                           \item{\code{noise_weights}: A logical; if \code{TRUE} then the approach utilizing a Hurdle model accounting for noise edges was utilized}
#'                           \item{\code{IC_selection}: A character string containing the specific information criteria used to select the optimal fit (i.e., 'BIC_model', 'BIC_mbc', 'ICL_mbc', 'Total_BIC', or 'Total_ICL')}
#'                           \item{\code{case_control}: A logical; if \code{TRUE} then the case/control approach was utilized}
#'                           \item{\code{DA_type}: A character string containing the specific deterministic annealing approach utilized (i.e., 'none', 'cooling', 'heating', or 'hybrid')}
#'                           \item{\code{priors}: A list of the prior hyperparameters used. See \code{\link[JANE]{specify_priors}} for definitions.}
#'                           }}
#' \item{\code{clustering_performance}}{ (only if \code{true_labels} is \code{!NULL}) A list with the following components: \itemize{
#'                           \item{\code{CER}: A list with two components: (i) \code{misclassified}: The indices of the misclassified actors in a minimum error mapping between the cluster labels and the known true cluster labels (i.e., \code{true_labels}) and (ii) \code{errorRate}: The error rate corresponding to a minimum error mapping between the cluster labels and the known true cluster labels (see \code{\link[mclust]{classError}} for details)}
#'                           \item{\code{ARI}: A numeric value containing the adjusted Rand index comparing the cluster labels and the known true cluster labels (see \code{\link[mclust]{adjustedRandIndex}} for details)}
#'                           \item{\code{NMI}: A numeric value containing the normalized mutual information comparing the cluster labels and the known true cluster labels (see \code{\link[aricode]{NMI}} for details)}
#'                           \item{\code{confusion_matrix}: A numeric table containing the confusion matrix comparing the cluster labels and the known true cluster labels.}
#'                           }}
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
#' # Summarize fit 
#' summary(res)
#' 
#' # Summarize fit and compare to true cluster labels
#' summary(res, true_labels = apply(sim_data$Z_U, 1, which.max))
#' 
#' # Summarize fit using starting values of EM algorithm
#' summary(res, initial_values = TRUE)
#' }
#' @method summary JANE
#' @exportS3Method summary JANE
summary.JANE <- function(object, true_labels = NULL, initial_values = FALSE, ...){
  
  if(!inherits(object, "JANE")){
    stop("Object is not of class JANE")
  }
  
  if(!(length(object$IC_out[, "selected"]) == 1) && all(object$IC_out[, "selected"] == 1)){
    stop("Unable to select the best K, D, and n_start using infomation criteria. See output matrix IC_out for infomation criteria values. Try different initialization approaches or only specify one K, D, and n_start.")
  }
  
  if(is.null(object$optimal_res) | is.null(object$optimal_starting)){
    stop("Unable to fit the model")
  }
  
  IC_selection <- object$input_params$IC_selection
  case_control <- object$input_params$case_control
  DA_type <- object$input_params$DA_type
  model <- object$input_params$model
  family <- object$input_params$family
  noise_weights <- object$input_params$noise_weights
  
  priors <- object$optimal_res$priors
  priors$a <- priors$a[1,]
  priors$c <- unname(priors$c)
  
  if(!initial_values){
    summary_data <- object$optimal_res
  } else {
    summary_data <- object$optimal_starting
  }

  if(!is.null(true_labels)){
    if(length(true_labels) != length(summary_data$cluster_labels)){
      stop("The length of the vector of true labels supplied does not match the number of actors in the fitted network")
    }
    CER <- mclust::classError(summary_data$cluster_labels, true_labels)
    ARI <- mclust::adjustedRandIndex(summary_data$cluster_labels, true_labels)
    NMI <- aricode::NMI(summary_data$cluster_labels, true_labels)
  }
  
  # return list
  out <- list(coefficients = list(beta_LR = summary_data$beta),
              p = summary_data$p,
              U = summary_data$U,
              mus = summary_data$mus,
              omegas = summary_data$omegas,
              Z_U = summary_data$prob_matrix,
              uncertainty = 1-apply(summary_data$prob_matrix, 1, max),
              cluster_labels = summary_data$cluster_labels)
  
  if(noise_weights){
    
    if(!initial_values){
      out$Z_W <- summary_data$prob_matrix_W
      out$Z_W <- cbind(out$Z_W, noise_edge_cluster_labels = unname(apply(out$Z_W[, 4:5], 1, which.max)))
    }
    
    out$coefficients$beta_GLM <- summary_data$beta2
    out$q_prob <- summary_data$q_prob
    
    if(family == "lognormal"){
      out$precision_weights <- summary_data$precision_weights
      out$precision_noise_weights <- summary_data$precision_noise_weights
    }
    
  }
  
  if(!initial_values){
    out$IC <- summary_data$IC
  }
  
  out$input_params <- list(model = model,
                           family = family,
                           noise_weights = noise_weights,
                           IC_selection = IC_selection,
                           case_control = case_control,
                           DA_type = DA_type,
                           priors = priors)
    
  if (!is.null(true_labels)){
    out$clustering_performance = list(CER = CER,
                                      ARI = ARI,
                                      NMI = NMI,
                                      confusion_matrix = table(summary_data$cluster_labels, true_labels,
                                                               dnn = list("JANE clusters", "Truth")))
  }
  
  return(structure(out,
                   class = "summary.JANE"))
  
}

#' @exportS3Method print summary.JANE
print.summary.JANE <- function(x, ...){
  
  if(!inherits(x, "summary.JANE")){
    stop("Object is not of class summary.JANE")
  }
  
  p <- x$p
  K <- length(p)
  D <- ncol(x$U)
  n_k <- rep(0, K)
  names(n_k) <- 1:K
  n_k[names(table(x$cluster_labels))] <- table(x$cluster_labels)

  cat("Input parameters:\n")
  cat("Model =", x$input_params$model, "\n")
  cat("Family =", x$input_params$family, "\n")
  cat("Noise weights =", x$input_params$noise_weights, "\n")
  cat("Information criteria used to select optimal model =", x$input_params$IC_selection, "\n")
  cat("Case-control approximation utilized =", x$input_params$case_control, "\n")
  cat("Type of deterministic annealing =", x$input_params$DA_type, "\n")
  
  cat("\nOptimal configuration selected:\n")
  cat("K =", K, "and D =", D, "\n")
  
  if(!is.null(x$IC)){
    cat("\nInformation criteria:\n")
    print(x$IC)
  }
  
  cat("\nClustering table:\n")
  print(n_k)
  
  if (!is.null(x$clustering_performance)){
    cat("\nClustering performance:\n")
    cat("ARI:", x$clustering_performance$ARI, "\n")
    cat("NMI:", x$clustering_performance$NMI, "\n")
    cat("CER:", x$clustering_performance$CER$errorRate, "\n")
    
    cat("\nConfusion matrix:\n")
    print(x$clustering_performance$confusion_matrix)
  }
  
  cat("\nRegression coefficient(s):\n")
  print(x$coefficients)
  
  cat("\nMixture weights:\n")
  names(p) <- 1:length(p)
  print(p)
  
  cat("\nMeans (KxD matrix):\n")
  print(x$mus)
  
  cat("\nPrecision Matrices (DxDxK array):\n")
  print(x$omegas)
  
}

#' @exportS3Method print JANE
print.JANE <- function(x, ...){
  
  if(!inherits(x, "JANE")){
    stop("Object is not of class JANE")
  }
  
  if(!(length(x$IC_out[, "selected"]) == 1) && all(x$IC_out[, "selected"] == 1)){
    stop("Unable to select the best K, D, and n_start using infomation criteria. See output matrix IC_out for infomation criteria values. Try different initialization approaches or only specify one K, D, and n_start.")
  }
  
  if(is.null(x$optimal_res) | is.null(x$optimal_starting)){
    stop("Unable to fit the model")
  }
  
  model <- x$input_params$model
  p <- x$optimal_res$p
  K <- length(p)
  D <- ncol(x$optimal_res$U)
  n_k <- rep(0, K)
  names(n_k) <- 1:K
  n_k[names(table(x$optimal_res$cluster_labels))] <- table(x$optimal_res$cluster_labels)
  
  cat(model, "latent space network clustering with", D, "dimensions and", K, "clusters of sizes", paste0(n_k, collapse = ", "))
  
  cat("\n\nAvailable components:\n")
  print(names(x))
  
  cat("\nAvailable components for optimal_res:\n")
  print(names(x$optimal_res))
  
}


#' Plot JANE fits
#' @description S3 plot method for object of class "\code{JANE}".
#' @param x An object of S3 \code{\link{class}} "\code{JANE}", a result of a call to \code{\link[JANE]{JANE}}.
#' @param type A character string to select the type of plot:
#'  \itemize{
#'   \item{'lsnc': plot the network using the estimated latent positions and color-code actors by cluster (default)}
#'   \item{'misclassified': (can only be used if \code{true_labels} is \code{!NULL}) similar to 'lsnc', but will color misclassified actors in black}
#'   \item{'uncertainty': similar to 'lsnc', but here the color gradient applied represents the actor-specific classification uncertainty}
#'   \item{'trace_plot': presents various trace plots across the iterations of the EM algorithm}
#'   }
#' @param true_labels (optional) A numeric, character, or factor vector of known true cluster labels. Must have the same length as number of actors in the fitted network. Need to account for potential isolates removed. 
#' @param initial_values A logical; if \code{TRUE} then plots fit using the starting parameters used in the EM algorithm (default is \code{FALSE}, i.e., the results after the EM algorithm is run are plotted).
#' @param density_type 	Choose from one of the following three options: 'contour' (default), 'hdr', 'image', and 'persp' indicating the density plot type.
#' @param swap_axes A logical; if \code{TRUE} will swap the x and y axes (default is \code{FALSE}).
#' @param rotation_angle A numeric value that rotates the estimated latent positions and contours of the multivariate normal distributions clockwise (or counterclockwise if \code{swap_axes = TRUE}) through the specified angle about the origin (default is 0 degrees). Only relevant when \code{D} (i.e., dimension of the latent space) \code{>= 2} and \code{type != 'trace_plot'}.
#' @param zoom A numeric value > 0 that controls the % magnification of the plot (default is 100%).
#' @param alpha_edge A numeric value in \code{[0,1]} that controls the transparency of the network edges (default is 0.1).
#' @param alpha_node A numeric value in \code{[0,1]} that controls the transparency of the actors in the network (default is 1).
#' @param main An optional overall title for the plot.
#' @param xlab An optional title for the x axis.
#' @param ylab An optional title for the y axis.
#' @param cluster_cols An optional vector of colors for the clusters. Must have a length of at least \eqn{K}.
#' @param remove_noise_edges (only applicable if \code{JANE} was run with \code{noise_weights = TRUE}) A logical; if \code{TRUE} will remove noise edges based on \eqn{\{h | \hat{Z}^{W}_{eh} = max(\hat{Z}^{W}_{e1}, \hat{Z}^{W}_{e2})\}} for \eqn{e = 1,\ldots,|E|}, where \eqn{\hat{Z}^{W}_{e1}} and \eqn{\hat{Z}^{W}_{e2}} are the estimated conditional probabilities that the \eqn{e^{th}} edge is a non-noise and noise edge, respectively (default is \code{FALSE}).
#' @param ... Unused.
#' @details
#'The classification of actors into specific clusters is based on a hard clustering rule of \eqn{\{h | \hat{Z}^{U}_{ih} = max_k \hat{Z}^{U}_{ik}\}}. Additionally, the actor-specific classification uncertainty is derived as 1 - \eqn{max_k \hat{Z}^{U}_{ik}}.
#'
#'The trace plot contains up to five unique plots tracking various metrics across the iterations of the EM algorithm, depending on the \code{\link[JANE]{JANE}} control parameter \code{termination_rule}:
#'  \itemize{
#'   \item{\code{termination_rule = 'prob_mat'}: Five plots will be presented. Specifically, in the top panel, the plot on the left presents the change in the absolute difference in \eqn{{\hat{Z}^U}} (i.e., the \eqn{N \times K} cluster membership probability matrix) between subsequent iterations and,  if \code{noise_weights = TRUE}, the change in the absolute difference in \eqn{{\hat{Z}^W}} (i.e., the \eqn{|E| \times 2} edge weight cluster membership probability matrix) between subsequent iterations. The exact quantile of the absolute difference plotted are presented in parentheses and determined by the \code{\link[JANE]{JANE}} control parameter \code{quantile_diff}. For example, the default control parameter \code{quantile_diff} = 1, so the values being plotted are the max absolute difference in \eqn{{\hat{Z}^U}} (and potentially \eqn{{\hat{Z}^W}})  between subsequent iterations. The plot on the right of the top panel presents the absolute difference in the cumulative average of the absolute change in \eqn{{\hat{Z}^U}} (and potentially \eqn{{\hat{Z}^W}})  and \eqn{\hat{U}} (i.e., the \eqn{N \times D} matrix of latent positions) across subsequent iterations (absolute change in \eqn{{\hat{Z}^U}}, \eqn{{\hat{Z}^W}}, and \eqn{\hat{U}} are computed in an identical manner as described previously). This metric is only tracked beginning at an iteration determined by the \code{n_its_start_CA} control parameter in \code{\link[JANE]{JANE}}. Note, this plot may be empty if the EM algorithm converges before the \code{n_its_start_CA}-th iteration. Finally, the bottom panel presents ARI, NMI, and CER values comparing the classifications between subsequent iterations, respectively. Specifically, at a given iteration we determine the classification of actors in clusters based on a hard clustering rule of \eqn{\{h | \hat{Z}^{U}_{ih} = max_k \hat{Z}^{U}_{ik}\}} and given these labels from two subsequent iterations, we compute and plot the ARI, NMI and CER.}
#'   \item{\code{termination_rule = 'Q'}: Plots generated are similar to those described in the previous bullet point. However, instead of tracking the change in \eqn{{\hat{Z}^U}} (and potentially \eqn{{\hat{Z}^W}}) over iterations, here the absolute difference in the objective function of the E-step evaluated using parameters from subsequent iterations is tracked. Furthermore, the cumulative average of the absolute change in \eqn{\hat{U}} is no longer tracked.}
#'  \item{\code{termination_rule \%in\% c('ARI', 'NMI', 'CER')}: Four plots will be presented. Specifically, the top left panel presents a plot of the absolute difference in the cumulative average of the absolute change in the specific \code{termination_rule} employed and \eqn{\hat{U}} across iterations. As previously mentioned, if the EM algorithm converges before the \code{n_its_start_CA}-th iteration then this will be an empty plot. Furthermore, the other three plots present ARI, NMI, and CER values comparing the classifications between subsequent iterations, respectively.}
#'   }
#'   
#' @seealso \code{\link[mclust]{surfacePlot}}, \code{\link[mclust]{adjustedRandIndex}}, \code{\link[mclust]{classError}},  \code{\link[aricode]{NMI}}  
#' 
#' @return A plot of the network or trace plot of the EM run.
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
#' # plot trace plot
#' plot(res, type = "trace_plot")
#'                    
#' # plot network
#' plot(res)
#' 
#' # plot network - misclassified
#' plot(res, type = "misclassified", true_labels = apply(sim_data$Z_U, 1, which.max))
#' 
#' # plot network - uncertainty and swap axes
#' plot(res, type = "uncertainty", swap_axes = TRUE)
#' 
#' # plot network - but only show contours of MVNs
#' plot(res, swap_axes = TRUE, alpha_edge = 0, alpha_node = 0)
#' 
#' # plot using starting values of EM algorithm
#' plot(res, initial_values = TRUE)
#'}
#' @method plot JANE
#' @exportS3Method plot JANE
plot.JANE <- function(x, type = "lsnc", true_labels, initial_values = FALSE,
                      zoom = 100, density_type = "contour", rotation_angle = 0,
                      alpha_edge = 0.1, alpha_node = 1, swap_axes = FALSE,
                      main, xlab, ylab, cluster_cols, remove_noise_edges = FALSE, ...){
  
  if(!inherits(x, "JANE")){
    stop("Object is not of class JANE")
  }
  
  if(!(length(x$IC_out[, "selected"]) == 1) && all(x$IC_out[, "selected"] == 1)){
    stop("Unable to select the best K, D, and n_start using infomation criteria. See output matrix IC_out for infomation criteria values. Try different initialization approaches or only specify one K, D, and n_start.")
  }
  
  if(is.null(x$optimal_res) | is.null(x$optimal_starting)){
    stop("Unable to fit the model")
  }
  
  if(!type %in% c("lsnc", "trace_plot", "misclassified", "uncertainty")){
    stop("Please provide one of the following for type: 'lsnc', 'trace_plot', 'misclassified', or 'uncertainty'")
  }

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  
  trace_plot <- FALSE
  uncertainty <- FALSE
  misclassified <- NULL
  
  if(type == "trace_plot"){
    
    trace_plot <- TRUE
    plot_data <- x$optimal_res
    
  } else {
    
    temp_A <- x$A
    
    if(!initial_values){
      
      plot_data <- x$optimal_res
      
      if(remove_noise_edges){
        if(!x$input_params$noise_weights){
          stop("Can only remove noise edges if JANE was run with noise_weights = TRUE")
        } else {
          Z_W <- plot_data$prob_matrix_W
          Z_W_labels <- apply(Z_W[, 4:5], 1, which.max)
          noise <- Z_W[Z_W_labels == 2, 1:2]
          temp_A[noise] <- 0L
          temp_A <- methods::as(temp_A, "dgCMatrix")
        }
      }
      
    } else {
      
      if(remove_noise_edges){
        message("All edges are assumed to be non-noise when initializing model parameters")
      }
      
      plot_data <- x$optimal_starting
    }
    
  }
  
  plot_data$model <- x$input_params$model
    
  if(!missing(cluster_cols)){
    if(length(unique(cluster_cols)) < length(plot_data$p)){
      stop("The number of unique cluster colors supplied is less than the number of clusters")
    } else {
      cluster_cols <- cluster_cols[1:length(plot_data$p)]
    }
  }
  
  if(type == "misclassified"){
    if(!missing(true_labels)){
      if(length(true_labels) != length(plot_data$cluster_labels)){
        stop("The length of the vector of true labels supplied does not match the number of actors in the fitted network")
      } else {
        misclassified <- mclust::classError(plot_data$cluster_labels, true_labels)$misclassified
      }
    } else {
      stop("Please supply vector of true labels using the true_labels argument when specifying type = 'misclassified'")
    } 
  }
  
  if(type == "uncertainty"){
    uncertainty <- TRUE
  }
  
  D <- ncol(plot_data$U)
  K <- length(plot_data$p)
  
  if(!(0<=alpha_edge & alpha_edge<=1)){
    stop("alpha_edge not in [0,1]")
  }
  
  if(!(0<=alpha_node & alpha_node<=1)){
    stop("alpha_node not in [0,1]")
  }
  
  if(!(0<=zoom)){
    stop("zoom needs to be >= 0")
  }
  
  if(density_type == "persp"){
    alpha_edge <- 0
    alpha_node <- 0
    zoom <- 100
    uncertainty <- FALSE
    misclassified <- NULL
  }
  
  if(missing(main)){
    main <- ifelse(!is.null(misclassified),
                   "Latent Space Network Clustering - Misclassified Actors",
                   ifelse(!uncertainty,
                          "Latent Space Network Clustering",
                          "Latent Space Network Clustering - Actor-specific Clustering Uncertainty"))
  }
  
  if(trace_plot){
    
    trace_plot(plot_data)
    
  } else {
    
    if(D == 1){
      
      means <- plot_data$mus[,1]
      vars <- apply(plot_data$omegas, 3, function(x){1.0/x})
      xlim <- c(min(plot_data$U[,1]), max(plot_data$U[,1])) + (100/zoom)*c(-1,1)
      ylim <- c(0, max(sapply(1:length(means), 
                              function(x){max(stats::dnorm(x = seq(from = xlim[1], to = xlim[2], by = 0.1),
                                                           mean = means[x], 
                                                           sd = sqrt(vars[x])))}))*1.1)
      if(missing(cluster_cols)){
        colors <- grDevices::rainbow(n = K)
      } else {
        colors <- cluster_cols
      }

      color_actors <- colors[plot_data$cluster_labels]
      
      dnorm_fun <- function(x, i){
        stats::dnorm(x = x, mean = means[i], sd = sqrt(vars[i]))
      }
      
      if(uncertainty & is.null(misclassified)){
        uncer <- round(1-apply(plot_data$prob_matrix, 1, max), 2)
        nf <- graphics::layout(
          matrix(c(1,2), ncol=2, byrow=TRUE), 
          widths = c(1,0.25)
        )
        graphics::par(mar=c(4, 4, 2, 0.25), oma=c(0,0,1,0), las=1)
      }
      
      if(missing(xlab)){
        xlab <- "Dim 1"
      }
      
      if(missing(ylab)){
        ylab <- ""
      }
      
      plot(1,
           type = "n",
           xlab = xlab, 
           ylab = ylab, 
           ylim = ylim, 
           xlim = xlim,
           main = main,
           cex.main = ifelse(!is.null(misclassified), 1.0, ifelse(!uncertainty, 1.0, 0.8)))
      
      if(!is.null(misclassified)){
        
        for (i in 1:K){
          graphics::curve(dnorm_fun(x, i = i),
                          from = xlim[1], to = xlim[2],
                          n = 1000,
                          add = T,
                          col = "black")
        }
        
        graphics::points(cbind(plot_data$U[,1], 0), pch = "|", 
                         cex = 1.2, 
                         col = scales::alpha(ifelse(1:nrow(plot_data$U) %in% misclassified == T, "black", "white"),
                                             alpha_node))
        graphics::legend("topright", legend= "Misclassified actor", pch = "|",
                         cex = 0.8) 
        
      } else {
        if(!uncertainty){
          
          for (i in 1:K){
            graphics::curve(dnorm_fun(x, i = i),
                            from = xlim[1], to = xlim[2],
                            n = 1000,
                            add = T,
                            col = colors[i])
          }
          
          graphics::points(cbind(plot_data$U[,1], 0), pch = "|", 
                           cex = 1, 
                           col = scales::alpha(color_actors, alpha_node))
          graphics::legend("topright", legend= paste0("Cluster ", 1:length(colors)), lty = 1,
                           col = colors,
                           cex = 0.8) 
          
        } else {
          
          for (i in 1:K){
            graphics::curve(dnorm_fun(x, i = i),
                            from = xlim[1], to = xlim[2],
                            n = 1000,
                            add = T,
                            col = "black")
          }
          
          if (length(unique(uncer)) > 1){
            break_points <- cut(uncer, breaks = seq(min(uncer) - 1e-6, max(uncer), length.out = 11))
          } else {
            break_points <- as.factor(uncer)
          }
          
          cols <- grDevices::heat.colors(length(levels(break_points)), alpha_node, rev = TRUE)
          graphics::points(cbind(plot_data$U[,1], 0), pch = "|", 
                           cex = 1, col = cols[break_points])
          graphics::par(mar = c(5, 0, 5, 5.5))
          graphics::image(1, 1:length(levels(break_points)), t(seq_along(levels(break_points))), 
                          col = cols, axes = FALSE, xlab = "")
          labels <- strsplit(levels(break_points), ",")
          labels <-  unlist(lapply(labels, function(x){
            p1 <- as.numeric(sub(pattern = "(\\()", x = x[1] , replacement = ""))
            p2 <- as.numeric(sub(pattern = "(\\])", x = x[2] , replacement = ""))
            p1 <- ifelse(p1<0, 0, p1)
            if(is.na(p2)){
              paste0(format(round(p1, 2), nsmall = 2))
            } else {
              paste0("(",format(round(p1, 2), nsmall = 2),", ", format(round(p2, 2), nsmall = 2), "]")
            }
          }))
          graphics::axis(4, at = 1:length(labels), labels = labels, cex.axis=0.70)
         
        }
      }
      
      
    } else if(D == 2){
      
      if(missing(xlab)){
        xlab <- ifelse(!swap_axes, "Dim 1", "Dim 2")
      }
      
      if(missing(ylab)){
        ylab <- ifelse(!swap_axes, "Dim 2", "Dim 1")
      }
      
      cluster_cols <-  if(missing(cluster_cols)){
        NULL
      } else {
        cluster_cols
      }
      
      plot_data(A = temp_A,
                data = plot_data, 
                misclassified = misclassified,
                zoom = zoom, 
                type = density_type,
                swap_axes = swap_axes,
                rotation_angle = rotation_angle,
                alpha_edge = alpha_edge, 
                alpha_node = alpha_node,
                uncertainty = uncertainty,
                cluster_cols = cluster_cols,
                main = main,
                xlab = xlab,
                ylab = ylab)
    
      
    } else {
      
      total_n_plots <- as.matrix(expand.grid(D_1 = 1:D, D_2 = 1:D))
      total_n_plots <- total_n_plots[total_n_plots[, "D_2"]>total_n_plots[, "D_1"],]
      total_n_plots <- total_n_plots[order(total_n_plots[, "D_1"]), ]
      
      for (i in 1:nrow(total_n_plots)){
        
        plot_data_temp <- plot_data
        plot_data_temp$U <- plot_data_temp$U[, unname(total_n_plots[i, ])]
        plot_data_temp$omegas <- array(apply(plot_data_temp$omega, 3, 
                                             function(y){
                                               y <- chol2inv(chol(y)) # this is the covariance matrix
                                               sigma <- y[total_n_plots[i,], total_n_plots[i,]] # remove irrelevant variables
                                               return(chol2inv(chol(sigma))) # convert back to precision matrix as that is what plot_data uses
                                             }, simplify = T),
                                       dim = c(2,2,K))
        plot_data_temp$mus <- plot_data_temp$mus[, unname(total_n_plots[i, ])]
        
        xlab_def <- ifelse(!swap_axes,
                       paste0("Dim ", unname(total_n_plots[i,1])), 
                       paste0("Dim ", unname(total_n_plots[i,2])))
        
        ylab_def <- ifelse(!swap_axes,
                       paste0("Dim ", unname(total_n_plots[i,2])), 
                       paste0("Dim ", unname(total_n_plots[i,1])))
        
        xlab_supp <- if(missing(xlab)){
          xlab_def
        } else {
          xlab
        }
        
        ylab_supp <- if(missing(ylab)){
          ylab_def
        } else {
          ylab
        }
        
        cluster_cols <-  if(missing(cluster_cols)){
          NULL
        } else {
          cluster_cols
        }
        
        plot_data(A = temp_A,
                  data = plot_data_temp, 
                  misclassified = misclassified,
                  zoom = zoom, 
                  type = density_type,
                  swap_axes = swap_axes,
                  rotation_angle = rotation_angle,
                  alpha_edge = alpha_edge, 
                  alpha_node = alpha_node,
                  uncertainty = uncertainty,
                  cluster_cols = cluster_cols,
                  main = main,
                  xlab = xlab_supp,
                  ylab = ylab_supp)

        
        if (i != nrow(total_n_plots)){
          readline(paste0("Hit <Return> for plot ", i+1, " of ", nrow(total_n_plots),":"))  
        }
        
      }
      
    }
  }
  
}

