
#' @export
plot_data <- function(A, data, zoom = 1, misclassified = NULL, type = "contour",
                          alpha_edge = 0.1, alpha_node = 1, title){
  
  U <- data$U
  Z <- data$cluster_labels
  mus <- data$mus
  omegas <- data$omegas
  model <- data$model
  undirected <- ifelse(model != "RSR", T, F)
  A_ig <- igraph::graph_from_adjacency_matrix(A,  mode = ifelse(undirected, "undirected", "directed"))
  k_elist <- igraph::as_edgelist(A_ig,  names= F)
  
  par <- list()
  par$pro <- rep(1, nrow(mus))
  par$mean <- t(mus)
  par$variance$sigma <- array(apply(omegas, 3, function(x){chol2inv(chol(x))}),
                              dim  = dim(omegas))
  
  mclust::surfacePlot(U, 
                      what = "density",
                      transformation = "none",
                      type = type,
                      parameters = par,
                      ylim = c(min(U[,2] ), max(U[,2])) + (1.0/zoom)*c(-1,1),
                      xlim = c(min(U[,1] ), max(U[,1])) + (1.0/zoom)*c(-1,1))
  
  if(undirected){
    
    graphics::segments(x0 = U[k_elist[,1],1],
                       x1 = U[k_elist[,2],1],
                       y0 = U[k_elist[,1],2],
                       y1 = U[k_elist[,2],2],
                       col= grDevices::gray(0.5, alpha_edge))
    
  } else {
    
    # get each arrow's length by converting x and y coords to inches
    units <- par(c('usr', 'pin'))
    x_to_inches <- with(units, pin[1L]/diff(usr[1:2])) # scale for x values to convert to inches
    y_to_inches <- with(units, pin[2L]/diff(usr[3:4])) # scale for y values to convert to inches
    
    distances <- matrix(data = 0.0, nrow = nrow(k_elist), ncol = 1)
    compute_dist(U = U %*% diag(c(x_to_inches, y_to_inches)), 
                 distances = distances, 
                 model = "NDH", 
                 X = matrix(0), 
                 indices = k_elist - 1,
                 downsampling = T) # compute L2 norm squared of rescaled U_i-U_j
    
    # find too short arrows causing warning (i.e. less than 1/1000 of an inch)
    idx_short_arrows <- which(sqrt(distances[,1])<0.001) # square root to get L2 norm
    
    # remove problem arrows
    if(length(idx_short_arrows)>0){
      k_elist <- k_elist[-idx_short_arrows, ]
    } 
    
    graphics::arrows(x0 = U[k_elist[,1],1],
           x1 = U[k_elist[,2],1],
           y0 = U[k_elist[,1],2],
           y1 = U[k_elist[,2],2],
           col= grDevices::gray(0.5, alpha_edge),
           length = 0.1,
           angle = 10) 
    
    
  }
  
  if(is.null(misclassified)){
    graphics::points(U,pch=16, cex = 0.8, col = alpha(Z, alpha_node))
  } else {
    graphics::points(U,pch=16, cex = 0.8, 
           col = scales::alpha(ifelse(1:nrow(A) %in% misclassified == T, "black", "white"),
                       alpha_node))
  }
  title(title)
}

