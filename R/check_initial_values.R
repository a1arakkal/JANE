
check_initial_values <- function(list_name,
                                 A,
                                 K,
                                 D,
                                 n_interior_knots,
                                 model,
                                 family,
                                 noise_weights){
  
  # Initial parameters dimension input check -----------------------------------
  params <- c("U", "omegas",  "mus", "p", "prob_matrix", "beta")
  correct_dim <- list(U = c(nrow(A),D),
                      omegas = c(D,D,K),
                      mus = c(K, D),
                      p = K,
                      prob_matrix = c(nrow(A), K))
  
  dim_init <- lapply(list_name[params[-length(params)]], function(x){if("array" %in% class(x)){dim(x)} else{length(x)}})
  check_dim <- sapply(names(dim_init), function(x){all(correct_dim[[x]] == dim_init[[x]])})
  
  if(!all(check_dim)){
    stop(paste0("Dimension of initial parameters incorrect for: ", 
                paste0(names(check_dim)[!check_dim], collapse = ", "), ". \n  ",
                "Dimension needs to be => ", 
                paste0(unname(sapply(names(check_dim)[!check_dim], function(x){ifelse(length(correct_dim[[x]]) >1,
                                                                                      paste0(x, ": array of dim (",paste0(correct_dim[[x]], collapse = ", "), ")"),
                                                                                      paste0(x, ": vector of length ", correct_dim[[x]]))})), 
                       collapse = ", "), 
                ".")
    )
  } 
  
  # Initial parameters values check ---------------------------------------------
  
  # check omegas
  check_omegas <- tryCatch(
    {
      apply(list_name[["omegas"]], 3, function(x){all(eigen(x)$values > 0)})
    },
    error = function(e) {
      rep(F,K)
    },
    warning = function(w) {
      rep(F,K)
    }
  ) 
  if(!all(check_omegas)){
    stop(paste0("Intial omega(s) are not positive definite. Check the following entries of the array: ",
                paste0("[,,", which(!check_omegas), "]", collapse = ", ")))
  }
  
  # check p
  check_p <- abs(sum(list_name[["p"]]) - 1) < 1e-12
  
  if(!check_p){
    stop("Intial vector p does not sum to 1")
  }
  
  # check prob_matrix
  check_prob_matrix <- abs(rowSums(list_name[["prob_matrix"]]) - 1) < 1e-12
  
  if(!all(check_prob_matrix)){
    stop("Rows of intial prob_matix do not sum to 1")
  }
  
  
  if(model == "NDH"){
    
    check_dim_beta <- length(list_name[["beta"]]) == 1
    
    if(!check_dim_beta){
      stop("Length of beta vector incorrect, needs to be 1")
    } 
    
  } else if (model == "RS"){
    
    check_dim_beta <- length(list_name[["beta"]]) == (1 + 1 + n_interior_knots)
    
    if(!check_dim_beta){
      stop(paste0("Length of beta vector incorrect, needs to be ", 1 + 1 + n_interior_knots))
    } 
    
    
  } else {
    
    check_dim_beta <- length(list_name[["beta"]]) == (1 + 2*(1 + n_interior_knots))
    
    if(!check_dim_beta){
      stop(paste0("Length of beta vector incorrect, needs to be ", 1 + 2*(1 + n_interior_knots)))
    } 
    
  } 
  
  if(noise_weights){

    if(family != "bernoulli"){
      
      if(model == "NDH"){
        
        check_dim_beta <- length(list_name[["beta2"]]) == 1
        
        if(!check_dim_beta){
          stop("Length of beta2 vector incorrect, needs to be 1")
        } 
        
      } else if (model == "RS"){
        
        check_dim_beta <- length(list_name[["beta2"]]) == (1 + 1 + n_interior_knots)
        
        if(!check_dim_beta){
          stop(paste0("Length of beta2 vector incorrect, needs to be ", 1 + 1 + n_interior_knots))
        } 
        
        
      } else {
        
        check_dim_beta <- length(list_name[["beta2"]]) == (1 + 2*(1 + n_interior_knots))
        
        if(!check_dim_beta){
          stop(paste0("Length of beta2 vector incorrect, needs to be ", 1 + 2*(1 + n_interior_knots)))
        } 
        
      } 
      
    }
    
    if(family %in% c("lognormal", "exp_lognormal")){
      
       check_dim_precision_weights <- length(list_name[["precision_weights"]]) == 1 & list_name[["precision_weights"]] > 0
      
      if(!check_dim_precision_weights){
        stop(paste0("precision_weights needs to be a postive scalar"))
      } 
       
       if (family == "lognormal"){
         
         check_dim_precision_noise_weights <- length(list_name[["precision_noise_weights"]]) == 1 & list_name[["precision_noise_weights"]] > 0
         
         if(!check_dim_precision_noise_weights){
           stop(paste0("precision_noise_weights needs to be a postive scalar"))
         } 
         
       }
       
    }
    
  }
  
}