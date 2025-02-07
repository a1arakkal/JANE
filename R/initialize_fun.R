
initialize_fun <- function(A, family, noise_weights, prob_matrix_W, priors, list_name, model, n_interior_knots, n_control, K, D){
  
  # Create new environment -----------------------------------------------------
  
  current <- invisible(rlang::new_environment(rlang::duplicate(list_name))) # need to duplicate as C++ modify in place will also change list_name which we don't want
  current$model <- model
  current$noise_weights <- noise_weights
  current$family <- family
  current$prob_matrix_W <- prob_matrix_W
  
  # Update new environment based on model --------------------------------------
  
  if(model == "NDH"){
    
    current$X <- matrix(0, 1, 1)
    
    current$fun_list <- list(update_prob_matrix = update_prob_matrix_DA,
                             update_p = update_p,
                             update_U = update_U,
                             update_mus_omegas = update_mus_omegas,
                             update_beta = update_beta,
                             log_Q = log_Q)
    
    
    if(!is.null(n_control)){
      current$fun_list$update_U <- update_U_CC
      current$fun_list$update_beta <- update_beta_CC
    }
    
    # priors
    if(is.null(priors)){
      
      current$priors <- list(
        a = matrix(rep(0,D), nrow = 1), # prior on mean of mus
        b = 1, # prior on precision of mus
        c = D + 1, # prior on df for omegas
        G = diag(D), # prior on sigma for omegas
        nu = rep(3,K), # prior for dirichlet
        e = 0, # prior mean on beta
        f = 1/100 # prior precision on beta
      )
      
    } else {
      
      current$priors <- priors
      
    }
    
  } else if(model == "RS"){
    
    # generate NS basis matrix 
    current$X <- splines::ns(x = (1.0/(nrow(A)-1.0))*rowSums(A), df = n_interior_knots + 1, intercept = F)
    
    current$fun_list <- list(update_prob_matrix = update_prob_matrix_DA,
                             update_p = update_p,
                             update_U = update_U_RE,
                             update_mus_omegas = update_mus_omegas,
                             update_beta = update_beta_RE,
                             log_Q = log_Q_RE)
    
    
    if(!is.null(n_control)){
      current$fun_list$update_U <- update_U_RE_CC
      current$fun_list$update_beta <- update_beta_RE_CC
    }
    
    # priors
    if(is.null(priors)){
      
      current$priors <- list(
        a = matrix(rep(0,D), nrow = 1), # prior on mean of mus
        b = 1, # prior on precision of mus
        c = D + 1, # prior on df for omegas
        G = diag(D), # prior on sigma for omegas
        nu = rep(3,K), # prior for dirichlet
        e = rep(0, 1 + ncol(current$X)), # prior mean on beta
        f = diag(c(1/100, rep(1/(2.5^2), ncol(current$X)))) # prior precision on beta
      )
      
    } else {
      
      current$priors <- priors
      
    }
    
  } else {
    
    # generate NS basis matrix 
    current$X <- cbind(splines::ns(x = (1.0/(nrow(A)-1.0))*rowSums(A), df = n_interior_knots + 1, intercept = F),
                       splines::ns(x = (1.0/(nrow(A)-1.0))*colSums(A), df = n_interior_knots + 1, intercept = F))
    
    current$fun_list <- list(update_prob_matrix = update_prob_matrix_DA,
                             update_p = update_p,
                             update_U = update_U_RE,
                             update_mus_omegas = update_mus_omegas,
                             update_beta = update_beta_RE,
                             log_Q = log_Q_RE)
    
    
    if(!is.null(n_control)){
      current$fun_list$update_U <- update_U_RE_CC
      current$fun_list$update_beta <- update_beta_RE_CC
    }
    
    # priors
    if(is.null(priors)){
      
      current$priors <- list(
        a = matrix(rep(0,D), nrow = 1), # prior on mean of mus
        b = 1, # prior on precision of mus
        c = D + 1, # prior on df for omegas
        G = diag(D), # prior on sigma for omegas
        nu = rep(3,K), # prior for dirichlet
        e = rep(0, 1 + ncol(current$X)), # prior mean on beta
        f = diag(c(1/100, rep(1/(2.5^2), ncol(current$X)))) # prior precision on beta
      )
      
    } else {
      
      current$priors <- priors
      
    }
    
  }
  
  if(noise_weights){
    
    current$fun_list$update_prob_matrix_W <- update_prob_matrix_W_DA
    current$fun_list$update_q_prob <- update_q_prob
    current$previous_prob_mat_W <- current$prob_matrix_W * 1.0
    
    if(is.null(current$priors$h) | is.null(current$priors$l)){
      current$priors$h <-  1 # prior parameter 1 for q
      current$priors$l <-  1 # prior parameter 2 for q
    }
    
    if(family != "bernoulli"){
      
      if(model == "NDH"){
        
        current$X2 <- matrix(0, 1, 1)
        
      } else if(model == "RS"){
        
        w <- prob_matrix_W[,3]
        
        # generate NS basis matrix for node strength
        degree <- tapply(rep(1, length(w)*2), rbind(prob_matrix_W[, c(1,2,3)],
                                                    prob_matrix_W[, c(2,1,3)])[,1], FUN = sum)
        
        node_strength <- tapply(c(w,w), rbind(prob_matrix_W[, c(1,2,3)],
                                              prob_matrix_W[, c(2,1,3)])[,1], FUN = sum)
        
        temp_scaled_node_strength <- node_strength/degree
        
        scaled_node_strength <- numeric(nrow(A))
        scaled_node_strength[as.numeric(names(temp_scaled_node_strength))] <- as.numeric(temp_scaled_node_strength)
        
        current$X2 <- splines::ns(x = scaled_node_strength, df = n_interior_knots + 1, intercept = F)
        
      } else {
        
        w <- prob_matrix_W[,3]
        
        # generate NS basis matrix for node strength
        degree_out <- tapply(rep(1, length(w)), prob_matrix_W[, 1], FUN = sum)
        
        degree_in <- tapply(rep(1, length(w)), prob_matrix_W[, 2], FUN = sum)
        
        node_strength_out <- tapply(w, prob_matrix_W[, 1], FUN = sum)
        
        node_strength_in <- tapply(w, prob_matrix_W[, 2], FUN = sum)
        
        temp_scaled_node_strength_out <- node_strength_out/degree_out
        
        temp_scaled_node_strength_in <- node_strength_in/degree_in
        
        scaled_node_strength_out <- numeric(nrow(A))
        scaled_node_strength_out[as.numeric(names(temp_scaled_node_strength_out))] <- as.numeric(temp_scaled_node_strength_out)
        
        scaled_node_strength_in <- numeric(nrow(A))
        scaled_node_strength_in[as.numeric(names(temp_scaled_node_strength_in))] <- as.numeric(temp_scaled_node_strength_in)
        
        current$X2  <- cbind(splines::ns(x = scaled_node_strength_out, df = n_interior_knots + 1, intercept = F),
                             splines::ns(x = scaled_node_strength_in, df = n_interior_knots + 1, intercept = F))
        
      }
      
      if(is.null(current$priors$e_1) | is.null(current$priors$f_2)){
        
        current$priors$e_1 <- rep(0, 1 + ncol(current$X2)) # prior mean on beta2
        current$priors$f_2 <- diag(c(1/100, rep(1/(2.5^2), ncol(current$X2)))) # prior precision on beta2
        
      } 
      
    }
    
    if(family == "lognormal"){
      
      if(is.null(current$priors$m_1) | is.null(current$priors$o_1) | is.null(current$priors$m_2) | is.null(current$priors$o_2) ){
        
        current$priors$m_1 <- 2 # prior parameter 1 for tau
        current$priors$o_1 <- 2 # prior parameter 2 for tau
        current$priors$m_2 <- 2 # prior parameter 1 for tau_noise
        current$priors$o_2 <- 2 # prior parameter 2 for tau_noise
        
      } 
      
    }
    
  }
  
  current$log_Q <- Inf
  current$previous_prob_mat <- current$prob_matrix * 1.0
  current$previous_U <- current$U * 1.0
  
  return(current)
  
}