
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
    current$X <- splines::ns(x = rowSums(A), df = n_interior_knots + 1, intercept = F)
    
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
    current$X <- cbind(splines::ns(x = rowSums(A), df = n_interior_knots + 1, intercept = F),
                       splines::ns(x = colSums(A), df = n_interior_knots + 1, intercept = F))
    
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
    current$priors$h <-  1 # prior parameter 1 for q
    current$priors$l <-  1 # prior parameter 2 for q
  }
  
  current$log_Q <- Inf
  current$previous_prob_mat <- current$prob_matrix * 1.0
  current$previous_U <- current$U * 1.0
  
  return(current)
  
}