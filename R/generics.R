#' @exportS3Method summary JANE
#' @export summary.JANE
summary.JANE <- function(object, true_labels = NULL, ...){
  
  if(!inherits(object, "JANE")){
    stop("Object is not of class JANE")
  }
  
  model <- object$optimal_res$model
  p <- object$optimal_res$p
  K <- length(p)
  D <- ncol(object$optimal_res$U)
  n_k <- table(object$optimal_res$cluster_labels)
  IC <- object$optimal_res$IC
  if(!is.null(true_labels)){
    if(length(true_labels) != length(object$optimal_res$cluster_labels)){
      stop("The length of the vector of true labels supplied does not match the number of actors in the fitted network")
    }
    CER <- mclust::classError(object$optimal_res$cluster_labels, true_labels)
    ARI <- mclust::adjustedRandIndex(object$optimal_res$cluster_labels, true_labels)
    NMI <- aricode::NMI(object$optimal_res$cluster_labels, true_labels)
  }
  
  
  cat("Optimal ", model, " model using ", object$IC_selection, " criteria:\n", sep = "")
  cat("K =", K, "and D =", D, "\n")
  
  cat("\nInformation criteria:\n")
  print(IC)
  
  cat("\nClustering table:")
  print(n_k)
  
  if (!is.null(true_labels)){
    cat("\nClustering performance:\n")
    cat("ARI:", ARI, "\n")
    cat("NMI:", NMI, "\n")
    cat("CER:", CER$errorRate, "\n")
    
    cat("\nConfusion matrix:\n")
    print(table(object$optimal_res$cluster_labels, true_labels,
                dnn = list("JANE clusters", "Truth")))
  }
  
  
  cat("\nMixing probabilities:\n")
  names(p) <- 1:length(p)
  print(p)
  
  cat("\nMeans (KxD matrix):\n")
  print(object$optimal_res$mus)
  
  cat("\nPrecision Matrices (DxDxK array):\n")
  print(object$optimal_res$omegas)
  
}


#' @exportS3Method print JANE
#' @export print.JANE
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
plot.JANE <- function(x, trace_plot = FALSE, true_labels = NULL,
                      zoom = 100, misclassified = NULL, type = "contour",
                      alpha_edge = 0.1, alpha_node = 1, swap_axes = FALSE, ...){
  
  D <- ncol(x$optimal_res$U)
  K <- length(x$optimal_res$p)
  if(!inherits(x, "JANE")){
    stop("Object is not of class JANE")
  }
  
  if(!is.null(true_labels)){
    if(length(true_labels) != length(x$optimal_res$cluster_labels)){
      stop("The length of the vector of true labels supplied does not match the number of actors in the fitted network")
    }
    CER <- mclust::classError(x$optimal_res$cluster_labels, true_labels)
  }
  
  if(!(0<=alpha_edge & alpha_edge<=1)){
    stop("alpha_edge not in [0,1]")
  }
  
  if(!(0<=alpha_node & alpha_node<=1)){
    stop("alpha_node not in [0,1]")
  }
  
  if(!(0<=zoom & zoom <=100)){
    stop("zoom not in [0,100]")
  }
  
  if(type == "persp"){
    alpha_edge <- 0
    alpha_node <- 0
    zoom <- 100
  }
  
  if(trace_plot){
    
    trace_plot(x$optimal_res)
    
  } else {
    
    if(D == 1){
      
      means <- x$optimal_res$mus[,1]
      vars <- apply(x$optimal_res$omegas, 3, function(x){1.0/x})
      xlim <- c(min(x$optimal_res$U[,1]), max(x$optimal_res$U[,1])) + (100/zoom)*c(-1,1)
      ylim <- c(0, max(sapply(1:length(means), 
                              function(x){max(stats::dnorm(x = seq(from = xlim[1], to = xlim[2], by = 0.1),
                                                           mean = means[x], 
                                                           sd = sqrt(vars[x])))}))*1.1)
      colors <- grDevices::rainbow(n = K)
      color_actors <- colors[x$optimal_res$cluster_labels]
      
      dnorm_fun <- function(x, i){
        stats::dnorm(x = x, mean = means[i], sd = sqrt(vars[i]))
      }
      
      plot(1,
           type ="n",
           xlab ="Dim 1", 
           ylab ="", 
           ylim = ylim, 
           xlim = xlim)
      
      for (i in 1:K){
        graphics::curve(dnorm_fun(x, i = i),
                        from = xlim[1], to = xlim[2],
                        n = 1000,
                        add = T,
                        col = colors[i])
      }
      
      if(!is.null(true_labels)){
        graphics::points(cbind(x$optimal_res$U[,1], 0), pch = "|", 
                         cex = 1, 
                         col = scales::alpha(ifelse(1:nrow(x$optimal_res$U) %in% CER$misclassified == T, "black", "white"),
                                             alpha_node))
        title(main = "Results of Latent Space Network Clustering using JANE - Misclassified Actors")
      } else {
        graphics::points(cbind(x$optimal_res$U[,1], 0), pch = "|", 
                         cex = 1, 
                         col = scales::alpha(color_actors, alpha_node))
        title(main = "Results of Latent Space Network Clustering using JANE")
      }
      
      
    } else if(D == 2){
      
      if(!is.null(true_labels)){
        plot_data(A = x$A,
                  data = x$optimal_res, 
                  misclassified = CER$misclassified,
                  zoom = zoom, 
                  type = type,
                  swap_axes = swap_axes,
                  alpha_edge = alpha_edge, 
                  alpha_node = alpha_node,
                  title = "Results of Latent Space Network Clustering using JANE - Misclassified Actors")
        if (!swap_axes){
          graphics::title(ylab = "Dim 2",
                          xlab = "Dim 1")
        } else {
          graphics::title(ylab = "Dim 1",
                          xlab = "Dim 2")
        }
        
      } else {
        plot_data(A = x$A,
                  data = x$optimal_res, 
                  zoom = zoom, 
                  type = type,
                  swap_axes = swap_axes,
                  alpha_edge = alpha_edge, 
                  alpha_node = alpha_node,
                  title = "Results of Latent Space Network Clustering using JANE")
        if (!swap_axes){
          graphics::title(ylab = "Dim 2",
                          xlab = "Dim 1")
        } else {
          graphics::title(ylab = "Dim 1",
                          xlab = "Dim 2")
        }
      }
      
    } else {
      
      total_n_plots <- as.matrix(expand.grid(D_1 = 1:D, D_2 = 1:D))
      total_n_plots <- total_n_plots[total_n_plots[, "D_2"]>total_n_plots[, "D_1"],]
      total_n_plots <- total_n_plots[order(total_n_plots[, "D_1"]), ]
      
      for (i in 1:nrow(total_n_plots)){
        
        plot_data_temp <- x$optimal_res
        plot_data_temp$U <- plot_data_temp$U[, unname(total_n_plots[i, ])]
        plot_data_temp$omegas <- array(apply(plot_data_temp$omega, 3, 
                                             function(y){
                                               y <- chol2inv(chol(y)) # this is the covariance matrix
                                               sigma <- y[total_n_plots[i,], total_n_plots[i,]] # remove irrelevant variables
                                               return(chol2inv(chol(sigma))) # convert back to precision matrix as that is what plot_data uses
                                             }, simplify = T),
                                       dim = c(2,2,K))
        plot_data_temp$mus <- plot_data_temp$mus[, unname(total_n_plots[i, ])]
        
        if(!is.null(true_labels)){
          plot_data(A = x$A,
                    data = plot_data_temp, 
                    misclassified = CER$misclassified,
                    zoom = zoom, 
                    type = type,
                    swap_axes = swap_axes,
                    alpha_edge = alpha_edge, 
                    alpha_node = alpha_node,
                    title = "Results of Latent Space Network Clustering using JANE - Misclassified Actors")
          if (!swap_axes){
            graphics::title(ylab = paste0("Dim ", unname(total_n_plots[i,2])),
                            xlab = paste0("Dim ", unname(total_n_plots[i,1])))
          } else {
            graphics::title(ylab = paste0("Dim ", unname(total_n_plots[i,1])),
                            xlab = paste0("Dim ", unname(total_n_plots[i,2])))
          }
          
          
        } else {
          plot_data(A = x$A,
                    data = plot_data_temp,
                    zoom = zoom, 
                    type = type,
                    swap_axes = swap_axes,
                    alpha_edge = alpha_edge, 
                    alpha_node = alpha_node,
                    title = "Results of Latent Space Network Clustering using JANE")
          if (!swap_axes){
            graphics::title(ylab = paste0("Dim ", unname(total_n_plots[i,2])),
                            xlab = paste0("Dim ", unname(total_n_plots[i,1])))
          } else {
            graphics::title(ylab = paste0("Dim ", unname(total_n_plots[i,1])),
                            xlab = paste0("Dim ", unname(total_n_plots[i,2])))
          }
        }
        
        if (i != nrow(total_n_plots)){
          readline(paste0("Hit enter for plot ", i+1, " of ", nrow(total_n_plots),":"))  
        }
        
      }
      
    }
  }
}
