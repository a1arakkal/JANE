
check_priors <- function(priors,
                         D,
                         K,
                         n_interior_knots,
                         model,
                         family,
                         noise_weights){
  
  # Check input priors list has all hyperparameters needed ---------------------
  
  if (!noise_weights){
    hyperparms <- c("a", "b", "c", "G", "nu", "e", "f")
  } else {
    if(family == "bernoulli"){
      hyperparms <- c("a", "b", "c", "G", "nu", "e", "f", "h", "l")
    } else if(family == "poisson"){
      hyperparms <- c("a", "b", "c", "G", "nu", "e", "f", "h", "l", "e_2", "f_2")
    } else if(family == "exp_lognormal"){
      hyperparms <- c("a", "b", "c", "G", "nu", "e", "f", "h", "l", "e_2", "f_2", "m_1", "o_1")
    } else {
      hyperparms <- c("a", "b", "c", "G", "nu", "e", "f", "h", "l", "e_2", "f_2", "m_1", "o_1", "m_2", "o_2")
    }
  }
 
  check_hyperparms <- hyperparms %in% names(priors)
  if(!all(check_hyperparms)){
    stop(paste0("Missing hyperparameter(s) for: ", paste0(hyperparms[!check_hyperparms], collapse = ", ")))
  } 
  
  priors$a <- priors$a[1,] 
  
  if(model == "NDH"){
    
    # input check
    if (!noise_weights){
      correct_dim <- list(a = D,
                          b = 1,
                          c = 1,
                          G = c(D, D),
                          nu = K,
                          e = 1,
                          f = 1)
    } else {
      if(family == "bernoulli"){
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1,
                            f = 1,
                            h = 1,
                            l = 1)
      } else if(family == "poisson"){
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1,
                            f = 1,
                            h = 1,
                            l = 1,
                            e_2 = 1,
                            f_2 = 1)
      } else if(family == "exp_lognormal"){
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1,
                            f = 1,
                            h = 1,
                            l = 1,
                            e_2 = 1,
                            f_2 = 1,
                            m_1 = 1,
                            o_1 = 1)
      } else {
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1,
                            f = 1,
                            h = 1,
                            l = 1,
                            e_2 = 1,
                            f_2 = 1,
                            m_1 = 1,
                            o_1 = 1,
                            m_2 = 1,
                            o_2 = 1)
      }
    }

    dim_priors <- lapply(priors, function(x){if("array" %in% class(x)){dim(x)} else{length(x)}})
    check_dim <- sapply(names(dim_priors), function(x){all(correct_dim[[x]] == dim_priors[[x]])})
    
    if(!all(check_dim)){
      stop(paste0("Dimension of hyperparameters incorrect for: ", 
                  paste0(names(check_dim)[!check_dim], collapse = ", "), ". \n  ",
                  "Dimension needs to be => ", 
                  paste0(unname(sapply(names(check_dim)[!check_dim], function(x){ifelse(length(correct_dim[[x]]) >1,
                                                                                        paste0(x, ": array of dim (",paste0(correct_dim[[x]], collapse = ", "), ")"),
                                                                                        paste0(x, ": vector of length ", correct_dim[[x]]))})), 
                         collapse = ", "), 
                  ".")
      )
    } 
    
  } else if (model == "RS"){
    
    # input check
    if (!noise_weights){
      correct_dim <- list(a = D,
                          b = 1,
                          c = 1,
                          G = c(D, D),
                          nu = K,
                          e = 1 + n_interior_knots + 1,
                          f = c(1 + n_interior_knots + 1, 1 + n_interior_knots + 1))
    } else {
      if(family == "bernoulli"){
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1 + n_interior_knots + 1,
                            f = c(1 + n_interior_knots + 1, 1 + n_interior_knots + 1),
                            h = 1,
                            l = 1)
      } else if(family == "poisson"){
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1 + n_interior_knots + 1,
                            f = c(1 + n_interior_knots + 1, 1 + n_interior_knots + 1),
                            h = 1,
                            l = 1,
                            e_2 = 1 + n_interior_knots + 1,
                            f_2 = c(1 + n_interior_knots + 1, 1 + n_interior_knots + 1))
      } else if(family == "exp_lognormal"){
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1 + n_interior_knots + 1,
                            f = c(1 + n_interior_knots + 1, 1 + n_interior_knots + 1),
                            h = 1,
                            l = 1,
                            e_2 = 1 + n_interior_knots + 1,
                            f_2 = c(1 + n_interior_knots + 1, 1 + n_interior_knots + 1),
                            m_1 = 1,
                            o_1 = 1)
      } else {
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1 + n_interior_knots + 1,
                            f = c(1 + n_interior_knots + 1, 1 + n_interior_knots + 1),
                            h = 1,
                            l = 1,
                            e_2 = 1 + n_interior_knots + 1,
                            f_2 = c(1 + n_interior_knots + 1, 1 + n_interior_knots + 1),
                            m_1 = 1,
                            o_1 = 1,
                            m_2 = 1,
                            o_2 = 1)
      }
    }
    
    dim_priors <- lapply(priors, function(x){if("array" %in% class(x)){dim(x)} else{length(x)}})
    check_dim <- sapply(names(dim_priors), function(x){all(correct_dim[[x]] == dim_priors[[x]])})
    
    if(!all(check_dim)){
      stop(paste0("Dimension of hyperparameters incorrect for: ", 
                  paste0(names(check_dim)[!check_dim], collapse = ", "), ". \n  ",
                  "Dimension needs to be => ", 
                  paste0(unname(sapply(names(check_dim)[!check_dim], function(x){ifelse(length(correct_dim[[x]]) >1,
                                                                                        paste0(x, ": array of dim (",paste0(correct_dim[[x]], collapse = ", "), ")"),
                                                                                        paste0(x, ": vector of length ", correct_dim[[x]]))})), 
                         collapse = ", "), 
                  ".")
      )
    }
    
    
  } else {
    
    # input check
    if (!noise_weights){
      correct_dim <- list(a = D,
                          b = 1,
                          c = 1,
                          G = c(D, D),
                          nu = K,
                          e = 1 + 2*(n_interior_knots + 1),
                          f = c(1 + 2*(n_interior_knots + 1), 1 + 2*(n_interior_knots + 1)))
      
    } else {
      if(family == "bernoulli"){
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1 + 2*(n_interior_knots + 1),
                            f = c(1 + 2*(n_interior_knots + 1), 1 + 2*(n_interior_knots + 1)),
                            h = 1,
                            l = 1)
      } else if(family == "poisson"){
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1 + 2*(n_interior_knots + 1),
                            f = c(1 + 2*(n_interior_knots + 1), 1 + 2*(n_interior_knots + 1)),
                            h = 1,
                            l = 1,
                            e_2 = 1 + 2*(n_interior_knots + 1),
                            f_2 = c(1 + 2*(n_interior_knots + 1), 1 + 2*(n_interior_knots + 1)))
      } else if(family == "exp_lognormal"){
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1 + 2*(n_interior_knots + 1),
                            f = c(1 + 2*(n_interior_knots + 1), 1 + 2*(n_interior_knots + 1)),
                            h = 1,
                            l = 1,
                            e_2 = 1 + 2*(n_interior_knots + 1),
                            f_2 = c(1 + 2*(n_interior_knots + 1), 1 + 2*(n_interior_knots + 1)),
                            m_1 = 1,
                            o_1 = 1)
      } else {
        correct_dim <- list(a = D,
                            b = 1,
                            c = 1,
                            G = c(D, D),
                            nu = K,
                            e = 1 + 2*(n_interior_knots + 1),
                            f = c(1 + 2*(n_interior_knots + 1), 1 + 2*(n_interior_knots + 1)),
                            h = 1,
                            l = 1,
                            e_2 = 1 + 2*(n_interior_knots + 1),
                            f_2 = c(1 + 2*(n_interior_knots + 1), 1 + 2*(n_interior_knots + 1)),
                            m_1 = 1,
                            o_1 = 1,
                            m_2 = 1,
                            o_2 = 1)
      }
    }
  
    dim_priors <- lapply(priors, function(x){if("array" %in% class(x)){dim(x)} else{length(x)}})
    check_dim <- sapply(names(dim_priors), function(x){all(correct_dim[[x]] == dim_priors[[x]])})
    
    if(!all(check_dim)){
      stop(paste0("Dimension of hyperparameters incorrect for: ", 
                  paste0(names(check_dim)[!check_dim], collapse = ", "), ". \n   ",
                  "Dimension needs to be => ", 
                  paste0(unname(sapply(names(check_dim)[!check_dim], function(x){ifelse(length(correct_dim[[x]]) >1,
                                                                                        paste0(x, ": array of dim (",paste0(correct_dim[[x]], collapse = ", "), ")"),
                                                                                        paste0(x, ": vector of length ", correct_dim[[x]]))})), 
                         collapse = ", "), 
                  ".")
      )
    } 
    
  } 
  
  
}
