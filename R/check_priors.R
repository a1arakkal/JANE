
#' @export
check_priors <- function(priors,
                         D,
                         K,
                         n_interior_knots,
                         model){
  
  # Check input priors list has all hyperparameters needed ---------------------
  
  hyperparms <- c("a", "b", "c", "G", "nu", "e", "f")
  
  check_hyperparms <- hyperparms %in% names(priors)
  if(!all(check_hyperparms)){
    stop(paste0("Missing hyperparameter(s) for: ", paste0(hyperparms[!check_hyperparms], collapse = ", ")))
  } 
  
  priors$a <- priors$a[1,] 
  
  if(model == "NDH"){
    
    # input check
    correct_dim <- list(a = D,
                        b = 1,
                        c = 1,
                        G = c(D, D),
                        nu = K,
                        e = 1,
                        f = 1)
    
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
    correct_dim <- list(a = D,
                        b = 1,
                        c = 1,
                        G = c(D, D),
                        nu = K,
                        e = 1 + n_interior_knots + 1,
                        f = c(1 + n_interior_knots + 1, 1 + n_interior_knots + 1))
    
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
    correct_dim <- list(a = D,
                        b = 1,
                        c = 1,
                        G = c(D, D),
                        nu = K,
                        e = 1 + 2*(n_interior_knots + 1),
                        f = c(1 + 2*(n_interior_knots + 1), 1 + 2*(n_interior_knots + 1)))
    
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
