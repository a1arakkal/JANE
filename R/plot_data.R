
plot_data <- function(A, data, zoom = 100, misclassified = NULL, type = "contour",  rotation_angle = 0,
                      alpha_edge = 0.1, alpha_node = 1, swap_axes = FALSE, uncertainty = FALSE, 
                      main = NULL, xlab = NULL, ylab = NULL){
  
  rotation_radians <- (pi/180)*rotation_angle
  rot_mat <- matrix(c(cos(rotation_radians),
                      sin(rotation_radians), 
                      -sin(rotation_radians),
                      cos(rotation_radians)), nrow = 2)
  
  U <- data$U %*% rot_mat
  Z <- data$cluster_labels
  mus <- data$mus
  omegas <- data$omegas
  model <- data$model
  undirected <- ifelse(model != "RSR", T, F)
  A_ig <- igraph::graph_from_adjacency_matrix(A,  mode = ifelse(undirected, "undirected", "directed"))
  k_elist <- igraph::as_edgelist(A_ig,  names= F)
  

  par <- list()
  par$pro <- rep(1, nrow(mus))
  par$mean <- t(mus %*% rot_mat)
  par$variance$sigma <- array(apply(omegas, 3, function(x){t(rot_mat) %*% chol2inv(chol(x)) %*% rot_mat}),
                              dim  = dim(omegas))
  

  
  if (swap_axes){
    par$mean <- par$mean[2:1, ]
    par$variance$sigma <- par$variance$sigma[2:1, 2:1, ]
    U <- U[, 2:1]
  } 
  
  if(uncertainty & is.null(misclassified)){
    uncer <- round(1-apply(data$prob_matrix, 1, max), 2)
    opar <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(opar), add = TRUE)
    nf <- graphics::layout(
      matrix(c(1,2), ncol=2, byrow=TRUE), 
      widths = c(1,0.25)
    )
    graphics::par(mar=c(4, 4, 2, 0.25), oma=c(0,0,1,0), las=1)
  }
  
  mclust::surfacePlot(data = U, 
                      what = "density",
                      transformation = "none",
                      type = type,
                      parameters = par,
                      swapAxes = FALSE,
                      ylim = c(min(U[,2]), max(U[,2])) + (100/zoom)*c(-1,1),
                      xlim = c(min(U[,1]), max(U[,1])) + (100/zoom)*c(-1,1), 
                      xlab = xlab,
                      ylab = ylab)
  graphics::title(main = main, 
                  cex.main = ifelse(!is.null(misclassified), 1.0, ifelse(!uncertainty, 1.0, 0.8)))
  
  
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
    if(!uncertainty){
      
      graphics::points(U,pch = 16, cex = 0.8, col = scales::alpha(Z, alpha_node))
      
    } else {
      
      if (length(unique(uncer)) > 1){
        break_points <- cut(uncer, breaks = seq(min(uncer) - 1e-6, max(uncer), length.out = 11))
      } else {
        break_points <- as.factor(uncer)
      }
      
      cols <- grDevices::heat.colors(length(levels(break_points)), alpha_node, rev = TRUE)
      graphics::points(U, pch = 16, cex = 0.8, col = cols[break_points])
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
  } else {
    graphics::points(U, pch = 16, cex = 0.8, 
                     col = scales::alpha(ifelse(1:nrow(A) %in% misclassified == T, "black", "white"),
                                         alpha_node))
  }
  
}

