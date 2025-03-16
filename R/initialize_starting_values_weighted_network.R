initialize_starting_values_weighted_network <- function(A,
                                                        model,
                                                        random_start,
                                                        control,
                                                        family,
                                                        prob_matrix_W,
                                                        noise_weights,
                                                        guess_noise_weights,
                                                        q_prob){
  
  # enter loop indicator
  retry <- T
  retry_count <- 0
  
  if (!random_start){
    
    starting_params <- tryCatch(
      
      {
        
        if(family == "lognormal"){
          
          prob_matrix_W_temp <- prob_matrix_W[log(prob_matrix_W[, 3]) > guess_noise_weights, , drop = FALSE]
          w <- prob_matrix_W_temp[,3]
          precision_noise_weights <- 1.0/stats::var(log(prob_matrix_W[, 3])[log(prob_matrix_W[, 3]) <= guess_noise_weights])
          precision_noise_weights <- ifelse(is.infinite(precision_noise_weights) | is.nan(precision_noise_weights) | is.na(precision_noise_weights),
                                            stats::rgamma(n = 1, 1, 1), precision_noise_weights)
          
          if(model == "NDH"){
            
            beta2 <- mean(log(w))
            precision_weights <- 1/stats::var(log(w))
            
          } else if(model == "RS"){
            
            # generate NS basis matrix for node strength
            degree <- tapply(rep(1, length(w)*2), rbind(prob_matrix_W_temp[, c(1,2,3), drop = FALSE],
                                                        prob_matrix_W_temp[, c(2,1,3), drop = FALSE])[,1], FUN = sum)
            
            node_strength <- tapply(c(w,w), rbind(prob_matrix_W_temp[, c(1,2,3), drop = FALSE],
                                                  prob_matrix_W_temp[, c(2,1,3), drop = FALSE])[,1], FUN = sum)
            
            temp_scaled_node_strength <- node_strength/degree
            
            scaled_node_strength <- numeric(nrow(A))
            scaled_node_strength[as.numeric(names(temp_scaled_node_strength))] <- as.numeric(temp_scaled_node_strength)

            X_basis <- splines::ns(x = scaled_node_strength, df = control$n_interior_knots + 1, intercept = F)
            
            distances <- matrix(0, nrow = nrow(prob_matrix_W_temp), ncol =  1 + (control$n_interior_knots + 1))
            
            compute_dist(U = matrix(0, nrow = nrow(prob_matrix_W_temp), ncol = 2),
                         distances = distances,
                         model = "RS", 
                         indices = cbind(prob_matrix_W_temp[,1]-1, prob_matrix_W_temp[,2]-1),
                         X = X_basis, 
                         downsampling = TRUE)
            fit <- suppressWarnings(stats::lm(log(w)~distances[, -ncol(distances)]))
            beta2 <- unname(fit$coefficients)
            precision_weights <- 1/(summary(fit)$sigma)^2
            
          } else {
            
            # generate NS basis matrix for node strength
            degree_out <- tapply(rep(1, length(w)), prob_matrix_W_temp[, 1], FUN = sum)
            
            degree_in <- tapply(rep(1, length(w)), prob_matrix_W_temp[, 2], FUN = sum)
            
            node_strength_out <- tapply(w, prob_matrix_W_temp[, 1], FUN = sum)
            
            node_strength_in <- tapply(w, prob_matrix_W_temp[, 2], FUN = sum)
            
            temp_scaled_node_strength_out <- node_strength_out/degree_out
            
            temp_scaled_node_strength_in <- node_strength_in/degree_in
            
            scaled_node_strength_out <- numeric(nrow(A))
            scaled_node_strength_out[as.numeric(names(temp_scaled_node_strength_out))] <- as.numeric(temp_scaled_node_strength_out)
            
            scaled_node_strength_in <- numeric(nrow(A))
            scaled_node_strength_in[as.numeric(names(temp_scaled_node_strength_in))] <- as.numeric(temp_scaled_node_strength_in)
            
            X_basis <- cbind(splines::ns(x = scaled_node_strength_out, df = control$n_interior_knots + 1, intercept = F),
                             splines::ns(x = scaled_node_strength_in, df = control$n_interior_knots + 1, intercept = F))
            
            
            distances <- matrix(0, nrow = nrow(prob_matrix_W_temp), ncol =  1 + 2*(control$n_interior_knots + 1))
            
            compute_dist(U = matrix(0, nrow = nrow(prob_matrix_W_temp), ncol = 2),
                         distances = distances,
                         model = "RSR", 
                         indices = cbind(prob_matrix_W_temp[,1]-1, prob_matrix_W_temp[,2]-1),
                         X = X_basis, 
                         downsampling = TRUE)
            fit <- suppressWarnings(stats::lm(log(w)~distances[, -ncol(distances)]))
            beta2 <- unname(fit$coefficients)
            precision_weights <- 1/(summary(fit)$sigma)^2
            
          }
          
          # Check for NAs in beta2
          if(any(is.na(beta2))){
            stop()
          }
          
          list(
            q_prob = q_prob,
            guess_noise_weights = guess_noise_weights,
            beta2 = beta2,
            precision_weights = precision_weights,
            precision_noise_weights = precision_noise_weights
          )
          
        } else if(family == "poisson"){
          
          prob_matrix_W_temp <- prob_matrix_W[prob_matrix_W[, 3] > guess_noise_weights, , drop = FALSE]
          w <- prob_matrix_W_temp[,3]
          
          if(model == "NDH"){
            
            beta2 <- log(mean(w))
            
          } else if(model == "RS"){
            
            # generate NS basis matrix for node strength
            degree <- tapply(rep(1, length(w)*2), rbind(prob_matrix_W_temp[, c(1,2,3), drop = FALSE],
                                                        prob_matrix_W_temp[, c(2,1,3), drop = FALSE])[,1], FUN = sum)
            
            node_strength <- tapply(c(w,w), rbind(prob_matrix_W_temp[, c(1,2,3), drop = FALSE],
                                                  prob_matrix_W_temp[, c(2,1,3), drop = FALSE])[,1], FUN = sum)
            
            temp_scaled_node_strength <- node_strength/degree
            
            scaled_node_strength <- numeric(nrow(A))
            scaled_node_strength[as.numeric(names(temp_scaled_node_strength))] <- as.numeric(temp_scaled_node_strength)

            X_basis <- splines::ns(x = scaled_node_strength, df = control$n_interior_knots + 1, intercept = F)
            
            distances <- matrix(0, nrow = nrow(prob_matrix_W_temp), ncol =  1 + (control$n_interior_knots + 1))
            
            compute_dist(U = matrix(0, nrow = nrow(prob_matrix_W_temp), ncol = 2),
                         distances = distances,
                         model = "RS", 
                         indices = cbind(prob_matrix_W_temp[,1]-1, prob_matrix_W_temp[,2]-1),
                         X = X_basis, 
                         downsampling = TRUE)
            fit <- suppressWarnings(stats::glm(w~distances[, -ncol(distances)], family = "poisson"))
            beta2 <- unname(fit$coefficients)
            
          } else {
            
            # generate NS basis matrix for node strength
            degree_out <- tapply(rep(1, length(w)), prob_matrix_W_temp[, 1], FUN = sum)
            
            degree_in <- tapply(rep(1, length(w)), prob_matrix_W_temp[, 2], FUN = sum)
            
            node_strength_out <- tapply(w, prob_matrix_W_temp[, 1], FUN = sum)
            
            node_strength_in <- tapply(w, prob_matrix_W_temp[, 2], FUN = sum)
            
            temp_scaled_node_strength_out <- node_strength_out/degree_out
            
            temp_scaled_node_strength_in <- node_strength_in/degree_in
            
            scaled_node_strength_out <- numeric(nrow(A))
            scaled_node_strength_out[as.numeric(names(temp_scaled_node_strength_out))] <- as.numeric(temp_scaled_node_strength_out)
            
            scaled_node_strength_in <- numeric(nrow(A))
            scaled_node_strength_in[as.numeric(names(temp_scaled_node_strength_in))] <- as.numeric(temp_scaled_node_strength_in)
            
            X_basis <- cbind(splines::ns(x = scaled_node_strength_out, df = control$n_interior_knots + 1, intercept = F),
                             splines::ns(x = scaled_node_strength_in, df = control$n_interior_knots + 1, intercept = F))
            
            
            distances <- matrix(0, nrow = nrow(prob_matrix_W_temp), ncol =  1 + 2*(control$n_interior_knots + 1))
            
            compute_dist(U = matrix(0, nrow = nrow(prob_matrix_W_temp), ncol = 2),
                         distances = distances,
                         model = "RSR", 
                         indices = cbind(prob_matrix_W_temp[,1]-1, prob_matrix_W_temp[,2]-1),
                         X = X_basis, 
                         downsampling = TRUE)
            
            fit <- suppressWarnings(stats::glm(w~distances[, -ncol(distances)], family = "poisson"))
            beta2 <- unname(fit$coefficients)
            
          }
          
          # Check for NAs in beta2
          if(any(is.na(beta2))){
            stop()
          }
          
          list(
            q_prob = q_prob,
            guess_noise_weights = guess_noise_weights,
            beta2 = beta2
          )
          
        } else {
          
          list(
            q_prob = q_prob,
            guess_noise_weights = guess_noise_weights
          )
          
        }
        
        
      },
      error = function(e) {
        NULL
      }
    )
    
    if(is.null(starting_params)){
      retry_count <- retry_count + 1
      if(control$verbose){
        message("Issues generating starting values, trying again.\n")
      }
    } else {
      retry <- F
    }
    
  }
  
  
  if(retry | random_start){
    
    if(retry & !random_start){
      if(control$verbose){
        message("Reached max GNN reattempts, switching to random values.\n")
      }
    } 
    
    starting_params <- list(
      q_prob = q_prob,
      guess_noise_weights = guess_noise_weights
    )
    
    if(family != "bernoulli"){
      if (model == "NDH"){
        
        starting_params$beta2 <- stats::rnorm(n = 1)
        
      } else if (model == "RS"){
        
        starting_params$beta2 <- stats::rnorm(n = 1 + (1 + control$n_interior_knots))
        
      } else {
        
        starting_params$beta2 <- stats::rnorm(n = 1 + 2*(1 + control$n_interior_knots))
        
      }
    }
    
    if(family == "lognormal"){
      starting_params$precision_weights <- stats::rgamma(n = 1, 1, 1)
      starting_params$precision_noise_weights <- stats::rgamma(n = 1, 1, 1)
    }
    
  }
  
  return(starting_params)
  
}

