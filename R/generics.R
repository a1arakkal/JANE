#' @exportS3Method summary JANE
#' @export summary.JANE
summary.JANE <- function(object, true_labels = NULL, initial_values = FALSE, ...){
  
  if(!inherits(object, "JANE")){
    stop("Object is not of class JANE")
  }
  
  IC_selection <- object$input_params$IC_selection
  case_control <- object$input_params$case_control
  DA_type <- object$input_params$DA_type

  if(!initial_values){
    summary_data <- object$optimal_res
  } else {
    summary_data <- object$optimal_starting
  }

  model <- summary_data$model
  
  if(!is.null(true_labels)){
    if(length(true_labels) != length(summary_data$cluster_labels)){
      stop("The length of the vector of true labels supplied does not match the number of actors in the fitted network")
    }
    CER <- mclust::classError(summary_data$cluster_labels, true_labels)
    ARI <- mclust::adjustedRandIndex(summary_data$cluster_labels, true_labels)
    NMI <- aricode::NMI(summary_data$cluster_labels, true_labels)
  }
  
  # return list
  out <- list(coefficients = summary_data$beta,
              p_k = summary_data$p,
              U = summary_data$U,
              mus = summary_data$mus,
              omegas = summary_data$omegas,
              Z = summary_data$prob_matrix,
              cluster_labels = summary_data$cluster_labels)
  
  if(!initial_values){
    out$IC <- summary_data$IC
  }
  
  out$input_params <- list(model = model,
                           IC_selection = IC_selection,
                           case_control = case_control,
                           DA_type = DA_type)
    
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
  
  p <- x$p_k
  K <- length(p)
  D <- ncol(x$U)
  n_k <- table(x$cluster_labels)

  cat("Input paramters:\n")
  cat("Model =", x$input_params$model, "\n")
  cat("Information criteria used to select optimal model =", x$input_params$IC_selection, "\n")
  cat("Case control =", x$input_params$case_control, "\n")
  cat("Type of deterministic annealing =", x$input_params$DA_type, "\n")
  
  cat("\nOptimal configuration selected:\n")
  cat("K =", K, "and D =", D, "\n")
  
  if(!is.null(x$IC)){
    cat("\nInformation criteria:\n")
    print(x$IC)
  }
  
  cat("\nClustering table:")
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
  
  cat("\nMixture probabilities:\n")
  names(p) <- 1:length(p)
  print(p)
  
  cat("\nMeans (KxD matrix):\n")
  print(x$mus)
  
  cat("\nPrecision Matrices (DxDxK array):\n")
  print(x$omegas)
  
}

#' @exportS3Method print JANE
print.JANE <- function(x, ...){
  
  model <- x$optimal_res$model
  p <- x$optimal_res$p
  K <- length(p)
  D <- ncol(x$optimal_res$U)
  n_k <- table(x$optimal_res$cluster_labels)
  
  if(!inherits(x, "JANE")){
    stop("Object is not of class JANE")
  }
  
  cat(model, "latent space network clustering with", D, "dimensions and", K, "clusters of sizes", paste0(n_k, collapse = ", "))
  
  cat("\n\nAvailable components:\n")
  print(names(x))
  
  cat("\nAvailable components for optimal_res:\n")
  print(names(x$optimal_res))
  
}


#' @exportS3Method plot JANE
#' @export plot.JANE
plot.JANE <- function(x, type = "lsnc", true_labels, initial_values = FALSE,
                      zoom = 100, density_type = "contour", 
                      alpha_edge = 0.1, alpha_node = 1, swap_axes = FALSE, ...){
  
  if(!inherits(x, "JANE")){
    stop("Object is not of class JANE")
  }
  
  if(!type %in% c("lsnc", "trace_plot", "misclassified", "uncertainty")){
    stop("Please provide one of the following for type: 'lsnc', 'trace_plot', 'misclassified', or 'uncertainty'")
  }
  
  trace_plot <- FALSE
  uncertainty <- FALSE
  misclassified <- NULL
  
  if(type == "trace_plot"){
    trace_plot <- TRUE
    plot_data <- x$optimal_res
  } else {
    if(!initial_values){
      plot_data <- x$optimal_res
    } else {
      plot_data <- x$optimal_starting
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
      colors <- grDevices::rainbow(n = K)
      color_actors <- colors[plot_data$cluster_labels]
      
      dnorm_fun <- function(x, i){
        stats::dnorm(x = x, mean = means[i], sd = sqrt(vars[i]))
      }
      
      if(uncertainty & is.null(misclassified)){
        uncer <- round(1-apply(plot_data$prob_matrix, 1, max), 2)
        opar <- graphics::par(no.readonly=TRUE)
        nf <- graphics::layout(
          matrix(c(1,2), ncol=2, byrow=TRUE), 
          widths = c(3,0.5)
        )
        graphics::par(mar=c(4, 4, 2, 0.5), oma=c(1,1,1,1), las=1)
      }
      
      plot(1,
           type = "n",
           xlab = "Dim 1", 
           ylab = "", 
           ylim = ylim, 
           xlim = xlim,
           main = ifelse(!is.null(misclassified),
                         "Latent Space Network Clustering - Misclassified Actors",
                         ifelse(!uncertainty,
                                "Latent Space Network Clustering",
                                "Latent Space Network Clustering - Actor-specific Clustering Uncertainty")),
           cex.main = ifelse(!is.null(misclassified), 1.0, ifelse(!uncertainty, 1.0, 0.8)))
      
      for (i in 1:K){
        graphics::curve(dnorm_fun(x, i = i),
                        from = xlim[1], to = xlim[2],
                        n = 1000,
                        add = T,
                        col = colors[i])
      }
      
      if(!is.null(misclassified)){
        graphics::points(cbind(plot_data$U[,1], 0), pch = "|", 
                         cex = 1, 
                         col = scales::alpha(ifelse(1:nrow(plot_data$U) %in% misclassified == T, "black", "white"),
                                             alpha_node))
      } else {
        if(!uncertainty){
          graphics::points(cbind(plot_data$U[,1], 0), pch = "|", 
                           cex = 1, 
                           col = scales::alpha(color_actors, alpha_node))
        } else {
          
          if (length(unique(uncer)) > 1){
            break_points <- cut(uncer, breaks = seq(min(uncer) - 1e-6, max(uncer), length.out = 11))
          } else {
            break_points <- as.factor(uncer)
          }
          
          cols <- grDevices::heat.colors(length(levels(break_points)), alpha_node, rev = TRUE)
          graphics::points(cbind(plot_data$U[,1], 0), pch = "|", 
                           cex = 1, col = cols[break_points])
          graphics::par(mar = c(5, 0, 5, 5))
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
          graphics::axis(4, at = 1:length(labels), labels = labels)
          on.exit(graphics::par(opar), add = TRUE)
          
        }
      }
      
      
    } else if(D == 2){
      
      plot_data(A = x$A,
                data = plot_data, 
                misclassified = misclassified,
                zoom = zoom, 
                type = density_type,
                swap_axes = swap_axes,
                alpha_edge = alpha_edge, 
                alpha_node = alpha_node,
                uncertainty = uncertainty,
                main = ifelse(!is.null(misclassified),
                               "Latent Space Network Clustering - Misclassified Actors",
                               ifelse(!uncertainty,
                                      "Latent Space Network Clustering",
                                      "Latent Space Network Clustering - Actor-specific Clustering Uncertainty")),
                xlab = ifelse(!swap_axes, "Dim 1", "Dim 2"),
                ylab = ifelse(!swap_axes, "Dim 2", "Dim 1"))
    
      
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
        
        plot_data(A = x$A,
                  data = plot_data_temp, 
                  misclassified = misclassified,
                  zoom = zoom, 
                  type = density_type,
                  swap_axes = swap_axes,
                  alpha_edge = alpha_edge, 
                  alpha_node = alpha_node,
                  uncertainty = uncertainty,
                  main = ifelse(!is.null(misclassified),
                                "Latent Space Network Clustering - Misclassified Actors",
                                ifelse(!uncertainty,
                                       "Latent Space Network Clustering",
                                       "Latent Space Network Clustering - Actor-specific Clustering Uncertainty")),
                  xlab = ifelse(!swap_axes,
                                paste0("Dim ", unname(total_n_plots[i,1])), 
                                paste0("Dim ", unname(total_n_plots[i,2]))),
                  ylab = ifelse(!swap_axes,
                                paste0("Dim ", unname(total_n_plots[i,2])), 
                                paste0("Dim ", unname(total_n_plots[i,1]))))

        
        if (i != nrow(total_n_plots)){
          readline(paste0("Hit enter for plot ", i+1, " of ", nrow(total_n_plots),":"))  
        }
        
      }
      
    }
  }
}
