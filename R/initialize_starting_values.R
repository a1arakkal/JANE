
#' @export
initialize_starting_values <- function(A, 
                                       K,
                                       D, 
                                       model,
                                       random_start,
                                       control,
                                       ...){
  
  N <- nrow(A)
  
  # enter loop indicator
  retry <- T
  retry_count <- 0
  
  if (!random_start){
    
    # Construct a row normalized matrix s.t. for the ith row their neighbors 
    # have a value 1/(deg_out_i).
    
    out_degree <- rowSums(A)
    A_matrix_row_norm <- Matrix::Diagonal(x = 1.0 / ifelse(out_degree>0, out_degree, 1.0)) %*% A
    
    # Construct a row normalized matrix s.t. for the ith row their neighbors 
    # have a value 1/(deg_out_i).
    in_degree <-  colSums(A)
    A_matrix_col_norm <- Matrix::Diagonal(x = 1.0 / ifelse(in_degree>0, in_degree, 1.0)) %*% t(A)
    
    while (retry & retry_count <= control$max_retry_GNN){
      
      # Construct a matrix of random values for starting U
      starting_U <- matrix(stats::rnorm(N*D, sd = control$sd_random_U_GNN), nrow = N, ncol = D)
      
      for (i in 1:control$n_its_GNN){
        
        # For each iteration the ith row of starting_U will be the average of the neighbors
        # values of the current starting_U 
        starting_U <- ((A_matrix_row_norm %*% starting_U) + (A_matrix_col_norm %*% starting_U))/2.0
        
      }
      
      starting_U <- as.matrix(starting_U)
      
      # Rescale U
      if(model == "NDH"){
        
        if (control$downsampling_GNN){
          
          # find row col indices of edges from upper diag
          indices <- summary(A)
          edges_indices <- indices[indices$j > indices$i, colnames(indices) != "x"] # will have sum(A)*0.5 rows 
          edges_indices[, "index"] <- ( (edges_indices$j -1) * (nrow(A)) ) + edges_indices$i # compute to column based-element-wise indices ( (col-1)*(N_total_row) ) + row
          edges_indices <- as.matrix(edges_indices)
          rownames(edges_indices) <- NULL
          colnames(edges_indices) <- NULL
          
          down_sample_nonLink_indices <- matrix(NA, nrow = nrow(edges_indices), ncol = ncol(edges_indices)) # create matrix for results
          
          while(any(is.na(down_sample_nonLink_indices))){
            
            rows_NA <- rowSums(is.na(down_sample_nonLink_indices)) > 0.0 # logical for whether or not rows need to be filled
            rsample_i <- round(stats::runif(n = sum(rows_NA)*10.0, min = 1.0, max = nrow(A))) # random draw of i
            rsample_j <- round(stats::runif(n = sum(rows_NA)*10.0, min = 1.0, max = nrow(A))) # random draw of j
            temp <- cbind(rsample_i[rsample_j > rsample_i], rsample_j[rsample_j > rsample_i]) # select where j > i
            index <- ( (temp[,2] -1) * (nrow(A)) ) + temp[,1] # compute the column based-element-wise indices
            # only retain index that is not edge index (i.e., setdiff(index, edges_indices$index) different than setdiff(edges_indices$index,index)
            # where the edge indices that are not in index is retained)
            temp <- cbind(temp[index %in% setdiff(index, edges_indices[,3]), , drop = F], index[index %in% setdiff(index, edges_indices[,3])])
            temp <- temp[!duplicated(temp[,3]), , drop = F] # remove duplicates
            # only retain index that is not down_sample_nonLink_indices so that we dont have duplicate values
            new_vals <- setdiff(temp[,3], down_sample_nonLink_indices[,3][!is.na(down_sample_nonLink_indices[,3])])
            
            if(length(new_vals)>0){ # only update if there are new_values to do so
              
              temp <- temp[temp[,3] %in% new_vals, , drop = F] # select the rows that will be added to down_sample_nonLink_indices
              
              if(sum(rows_NA) >= nrow(temp)){
                down_sample_nonLink_indices[which(rows_NA)[1:nrow(temp)], ] <- temp # if rows to update is more than nrow temp then only nrow temp of down_sample_nonLink_indices
              } else {
                down_sample_nonLink_indices[rows_NA, ] <- temp[1:sum(rows_NA), ] # else update all na rows
              }
              
            }
            
          }
          
          distances <- matrix(data = 0.0, nrow = nrow(edges_indices)*2, ncol = 1) # balanced design so nrow edges_indices*2
          compute_dist(U = starting_U, 
                       distances = distances, model = "NDH", 
                       X = matrix(0), 
                       indices = rbind(edges_indices[,-3]-1, down_sample_nonLink_indices[,-3]-1),
                       downsampling = T)
          
          res <- stats::optim(par = rep(0, 2),
                              fn = log_like_C,
                              gr = gradient_C,
                              lower = c(-Inf, -Inf),
                              upper = c(Inf, 0),
                              method = "L-BFGS-B",
                              X = cbind(rep(1, nrow(distances)), distances), 
                              y = A[c(edges_indices[,3], down_sample_nonLink_indices[,3])])
          
        } else {
          
          distances <- matrix(data = 0.0, nrow = (N-1)*N*0.5, ncol = 1) # for undirected network we only need ((N-1)*N)/2
          compute_dist(U = starting_U, distances = distances, model = "NDH", X = matrix(0), indices = matrix(0), downsampling = F)
          
          res <- stats::optim(par = rep(0, 2),
                              fn = log_like_C,
                              gr = gradient_C,
                              lower = c(-Inf, -Inf),
                              upper = c(Inf, 0),
                              method = "L-BFGS-B",
                              X = cbind(rep(1, nrow(distances)), distances), 
                              y = A[lower.tri(A, diag = F)])
        }
        
        scale_U <- ifelse(res$par[length(res$par)] == 0.0, 1.0, sqrt(-1.0*res$par[length(res$par)]))
        starting_U <- starting_U*scale_U
        starting_beta <- res$par[-length(res$par)]
        
      } else if (model == "RS"){
        
        # generate NS basis matrix for degree
        X_basis <- splines::ns(x = rowSums(A), df = control$n_interior_knots + 1, intercept = F)
        
        if (control$downsampling_GNN){
          
          # find row col indices of edges from upper diag
          indices <- summary(A)
          edges_indices <- indices[indices$j > indices$i, colnames(indices) != "x"] # will have sum(A)*0.5 rows 
          edges_indices[, "index"] <- ( (edges_indices$j -1) * (nrow(A)) ) + edges_indices$i # compute to column based-element-wise indices ( (col-1)*(N_total_row) ) + row
          edges_indices <- as.matrix(edges_indices)
          rownames(edges_indices) <- NULL
          colnames(edges_indices) <- NULL
          
          down_sample_nonLink_indices <- matrix(NA, nrow = nrow(edges_indices), ncol = ncol(edges_indices)) # create matrix for results
          
          while(any(is.na(down_sample_nonLink_indices))){
            
            rows_NA <- rowSums(is.na(down_sample_nonLink_indices)) > 0.0 # logical for whether or not rows need to be filled
            rsample_i <- round(stats::runif(n = sum(rows_NA)*10.0, min = 1.0, max = nrow(A))) # random draw of i
            rsample_j <- round(stats::runif(n = sum(rows_NA)*10.0, min = 1.0, max = nrow(A))) # random draw of j
            temp <- cbind(rsample_i[rsample_j > rsample_i], rsample_j[rsample_j > rsample_i]) # select where j > i
            index <- ( (temp[,2] -1) * (nrow(A)) ) + temp[,1] # compute the column based-element-wise indices
            # only retain index that is not edge index (i.e., setdiff(index, edges_indices$index) different than setdiff(edges_indices$index,index)
            # where the edge indices that are not in index is retained)
            temp <- cbind(temp[index %in% setdiff(index, edges_indices[,3]), , drop = F], index[index %in% setdiff(index, edges_indices[,3])])
            temp <- temp[!duplicated(temp[,3]), , drop = F] # remove duplicates
            # only retain index that is not down_sample_nonLink_indices so that we dont have duplicate values
            new_vals <- setdiff(temp[,3], down_sample_nonLink_indices[,3][!is.na(down_sample_nonLink_indices[,3])])
            
            if(length(new_vals)>0){ # only update if there are new_values to do so
              
              temp <- temp[temp[,3] %in% new_vals, , drop = F] # select the rows that will be added to down_sample_nonLink_indices
              
              if(sum(rows_NA) >= nrow(temp)){
                down_sample_nonLink_indices[which(rows_NA)[1:nrow(temp)], ] <- temp # if rows to update is more than nrow temp then only nrow temp of down_sample_nonLink_indices
              } else {
                down_sample_nonLink_indices[rows_NA, ] <- temp[1:sum(rows_NA), ] # else update all na rows
              }
              
            }
            
          }
          
          distances <- matrix(data = 0.0, nrow = nrow(edges_indices)*2, ncol = 1 + (control$n_interior_knots + 1)) # balanced design so edges_indices*2
          compute_dist(U = starting_U, 
                       distances = distances,
                       model = "RS", 
                       X = X_basis,
                       indices = rbind(edges_indices[,-3]-1, down_sample_nonLink_indices[,-3]-1),
                       downsampling = T)
          
          res <- stats::optim(par = rep(0, 1 + 1 + (control$n_interior_knots + 1)),
                              fn = log_like_C,
                              gr = gradient_C,
                              lower = rep(-Inf, 1 + 1 + (control$n_interior_knots + 1)),
                              upper = c(rep(Inf, 1 + (control$n_interior_knots + 1)), 0),
                              method = "L-BFGS-B",
                              X = cbind(rep(1, nrow(distances)), distances), 
                              y = A[c(edges_indices[,3], down_sample_nonLink_indices[,3])])
          
        } else {
          
          distances <- matrix(data = 0.0, nrow = (N-1)*N*0.5, ncol = 1 + (control$n_interior_knots + 1)) # for undirected network we only need ((N-1)*N)/2
          compute_dist(U = starting_U, distances = distances, model = "RS", indices = matrix(0), X = X_basis, downsampling = F)
          
          res <- stats::optim(par = rep(0, 1 + 1 + (control$n_interior_knots + 1)),
                              fn = log_like_C,
                              gr = gradient_C,
                              lower = rep(-Inf, 1 + 1 + (control$n_interior_knots + 1)),
                              upper = c(rep(Inf, 1 + (control$n_interior_knots + 1)), 0),
                              method = "L-BFGS-B",
                              X = cbind(rep(1, nrow(distances)), distances), 
                              y = A[lower.tri(A, diag = F)])
          
        }
        
        scale_U <- ifelse(res$par[length(res$par)] == 0.0, 1.0, sqrt(-1.0*res$par[length(res$par)]))
        starting_U <- starting_U*scale_U
        starting_beta <- res$par[-length(res$par)]
        
      } else {
        
        # generate NS basis matrix for in and out degree
        X_basis <- cbind(splines::ns(x = rowSums(A), df = control$n_interior_knots + 1, intercept = F),
                         splines::ns(x = colSums(A), df = control$n_interior_knots + 1, intercept = F))
        
        if (control$downsampling_GNN){
          
          # find row col indices of edges from upper diag
          indices <- summary(A)
          edges_indices <- indices[indices$j != indices$i, colnames(indices) != "x"] # will have sum(A) rows 
          edges_indices[, "index"] <- ( (edges_indices$j -1) * (nrow(A)) ) + edges_indices$i # compute to column based-element-wise indices ( (col-1)*(N_total_row) ) + row
          edges_indices <- as.matrix(edges_indices)
          rownames(edges_indices) <- NULL
          colnames(edges_indices) <- NULL
          
          down_sample_nonLink_indices <- matrix(NA, nrow = nrow(edges_indices), ncol = ncol(edges_indices)) # create matrix for results
          
          while(any(is.na(down_sample_nonLink_indices))){
            
            rows_NA <- rowSums(is.na(down_sample_nonLink_indices)) > 0.0 # logical for whether or not rows need to be filled
            rsample_i <- round(stats::runif(n = sum(rows_NA)*10.0, min = 1.0, max = nrow(A))) # random draw of i
            rsample_j <- round(stats::runif(n = sum(rows_NA)*10.0, min = 1.0, max = nrow(A))) # random draw of j
            temp <- cbind(rsample_i[rsample_j != rsample_i], rsample_j[rsample_j != rsample_i]) # select where j != i
            index <- ( (temp[,2] -1) * (nrow(A)) ) + temp[,1] # compute the column based-element-wise indices
            # only retain index that is not edge index (i.e., setdiff(index, edges_indices$index) different than setdiff(edges_indices$index,index)
            # where the edge indices that are not in index is retained)
            temp <- cbind(temp[index %in% setdiff(index, edges_indices[,3]), , drop = F], index[index %in% setdiff(index, edges_indices[,3])])
            temp <- temp[!duplicated(temp[,3]), , drop = F] # remove duplicates
            # only retain index that is not down_sample_nonLink_indices so that we dont have duplicate values
            new_vals <- setdiff(temp[,3], down_sample_nonLink_indices[,3][!is.na(down_sample_nonLink_indices[,3])])
            
            if(length(new_vals)>0){ # only update if there are new_values to do so
              
              temp <- temp[temp[,3] %in% new_vals, , drop = F] # select the rows that will be added to down_sample_nonLink_indices
              
              if(sum(rows_NA) >= nrow(temp)){
                down_sample_nonLink_indices[which(rows_NA)[1:nrow(temp)], ] <- temp # if rows to update is more than nrow temp then only nrow temp of down_sample_nonLink_indices
              } else {
                down_sample_nonLink_indices[rows_NA, ] <- temp[1:sum(rows_NA), ] # else update all na rows
              }
              
            }
            
          }
          
          distances <- matrix(data = 0.0, nrow = nrow(edges_indices)*2, ncol = 1 + 2*(control$n_interior_knots + 1)) # balanced design so edges_indices*2
          compute_dist(U = starting_U, 
                       distances = distances,
                       model = "RSR", 
                       X = X_basis,
                       indices = rbind(edges_indices[,-3]-1, down_sample_nonLink_indices[,-3]-1),
                       downsampling = T)
          
          res <- stats::optim(par = rep(0, 1 + 1 + 2*(control$n_interior_knots + 1)),
                              fn = log_like_C,
                              gr = gradient_C,
                              lower = rep(-Inf, 1 + 1 + 2*(control$n_interior_knots + 1)),
                              upper = c(rep(Inf, 1 + 2*(control$n_interior_knots + 1)), 0),
                              method = "L-BFGS-B",
                              X = cbind(rep(1, nrow(distances)), distances), 
                              y = A[c(edges_indices[,3], down_sample_nonLink_indices[,3])])
          
        } else {
          
          distances <- matrix(data = 0.0, nrow = (N-1)*N, ncol = 1 + 2*(control$n_interior_knots + 1)) # for directed network we need ((N-1)*N)
          compute_dist(U = starting_U, distances = distances, model = "RSR", indices = matrix(0), X = X_basis, downsampling = F)
          
          temp_y <- as.matrix(A) 
          diag(temp_y) <- NA
          temp_y <- c(t(temp_y))
          temp_y <- temp_y[!is.na(temp_y)]
          
          res <- stats::optim(par = rep(0, 1 + 1 + 2*(control$n_interior_knots + 1)),
                              fn = log_like_C,
                              gr = gradient_C,
                              lower = rep(-Inf, 1 + 1 + 2*(control$n_interior_knots + 1)),
                              upper = c(rep(Inf, 1 + 2*(control$n_interior_knots + 1)), 0),
                              method = "L-BFGS-B",
                              X = cbind(rep(1, nrow(distances)), distances), 
                              y = temp_y)
        }
        
        scale_U <- ifelse(res$par[length(res$par)] == 0.0, 1.0, sqrt(-1.0*res$par[length(res$par)]))
        starting_U <- starting_U*scale_U
        starting_beta <- res$par[-length(res$par)]
        
      }
      
      
      # Step4: Run K-means algo for GMM based on starting U
      starting_params <- stats::kmeans(x = starting_U,
                                       centers = K,
                                       iter.max = 10,
                                       nstart = 5)
      
      starting_mus <- starting_params$centers
      clust_assignments <- starting_params$cluster
      n_k <- tabulate(clust_assignments)
      
      if(length(n_k) == K & all(n_k>1)){
        check_omegas <- tryCatch(
          {
            sapply(1:K, function(x){all(eigen(chol2inv(chol(stats::var(starting_U[clust_assignments == x, ]))))$values > 0.0)})
          },
          error = function(e) {
            F
          },
          warning = function(w) {
            F
          }
        ) 
      } else {
        check_omegas <- F
      }
      
      if(all(check_omegas)){
        
        starting_omegas <- array(data = NA, dim = c(D,D,K))
        
        for (i in 1:K){
          starting_omegas[,,i] <- chol2inv(chol(stats::var(starting_U[clust_assignments == i, ])))
        }
        
        # for starting omegas we wont use sample precision of each cluster as
        # observations within a cluster will be very similar and thus will have a high precision
        # so we will just start with the identity matrix
        # starting_omegas <- array(data = diag(D), dim = c(D,D,K))
        
        starting_prob_matrix <- matrix(0.0, nrow = N, ncol = K)
        starting_prob_matrix[cbind(1:N, clust_assignments)] <- 1.0
        
        starting_params <- list(
          
          U = starting_U,
          omegas = starting_omegas,
          mus = starting_mus,
          p = colSums(starting_prob_matrix)/sum(colSums(starting_prob_matrix)),
          prob_matrix = starting_prob_matrix,
          beta = starting_beta
          
        )
        
        retry <- F
        
      } else {
        
        retry_count <- retry_count + 1
        message("Intial omega(s) are not positive definite, trying again with new random U matrix.\n")
        
      }
      
    }
    
  }
  
  if(retry | random_start){
    
    if(retry & !random_start){
      message("Reached max re-try iterrations, switching to random values.\n")
    } 
    
    starting_params <- list(
      
      U = matrix(stats::rnorm(nrow(A)*D), nrow = nrow(A), ncol = D),
      omegas = stats::rWishart(n = K, df = D+1, Sigma = diag(D)),
      mus = matrix(stats::rnorm(K*D), nrow = K, ncol = D),
      p = extraDistr::rdirichlet(n = 1, rep(3,K))[1,],
      prob_matrix = extraDistr::rdirichlet(n = nrow(A), alpha = rep(1, K))
      
    )
    
    if (model == "NDH"){
      
      starting_params$beta <- stats::rnorm(n = 1)
      
    } else if (model == "RS"){
      
      starting_params$beta <- stats::rnorm(n = 1 + (1 + control$n_interior_knots))
      
    } else {
      
      starting_params$beta <- stats::rnorm(n = 1 + 2*(1 + control$n_interior_knots))
      
    }
    
  }
  
  return(starting_params)
  
}

#' @useDynLib JANE   
log_like_C <- function(par, X, y) {
  .Call('_JANE_log_like_C', PACKAGE = 'JANE', par, X, y)
}

#' @useDynLib JANE   
gradient_C <- function(par, X, y) {
  .Call('_JANE_gradient_C', PACKAGE = 'JANE', par, X, y)
}

#' @useDynLib JANE   
compute_dist <- function(U, distances, model, X, indices, downsampling) {
  invisible(.Call('_JANE_compute_dist', PACKAGE = 'JANE', U, distances, model, X, indices, downsampling))
}

