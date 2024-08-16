
#' @export
EM <- function(A,
               D = 2,
               K = 2,
               model,
               initialization = "GNN", # random, GNN, or user supplied
               case_control = F,
               DA_type = "none", # none, cooling, heating, hybrid 
               seed = 2024, 
               control = list()){
  
  
  con <- list(
    max_its = 1e3, # max iteration of EM algorithm, 
    min_its = 10,  # min iteration of EM algorithm,
    priors = NULL,  # priors,
    n_interior_knots = 5, # number of interior knots used in natural splines for RS and RSR models
    termination_rule = "prob_mat", # termination rule
    tolerance = 1e-3, # tolerance for termination rules Q and prob_mat
    tolerance_ARI = 0.999, # tolerance for termination rule ARI
    tolerance_NMI = 0.999, # tolerance for termination rule NMI
    tolerance_CER = 0.01,  # tolerance for termination rule CER
    n_its_start_CA = 20, # what iteration to start computing moving averages (cumulative average of U not tracked when Q is used as termination rule)
    tolerance_diff_CA = 1e-3, # tolerance for change in cumulative average of metric and U (cumulative average of U not tracked when Q is used as termination rule)
    consecutive_diff_CA = 5, # number of consecutive instances where change in cumulative average is less than tolerance_diff_CA 
    quantile_diff = 1, # quantile for difference in prob_mat and U matrices between iterations
    beta_temp_schedule = 1, # Deterministic annealing temperature schedule
    n_control = NULL, # number of controls to be used in case-control approach
    n_start = 5, # number of starting values to try
    max_retry = 5, # number of max attempts to try if starting values causes issues with EM algo 
    IC_selection = "Total_ICL", # Information criteria to use for selection between different starting values, K, and D
    sd_random_U_GNN = 1, # standard deviation used for draws from random normal used in graphical neural network (GNN) starting values approach
    max_retry_GNN = 5, # number of max attempts for GNN approach before switching to random starting values 
    n_its_GNN = 2, # number of iterations for GNN approach
    downsampling_GNN = T # logical for whether or not to use downsampling s.t. number of links and non-links are balanced for GNN logistic regression approach
  )
  
  cl <- match.call()
  cl$A <- quote(A)
  cl$case_control <- NULL
  cl$seed <- NULL
  cl$DA_type <- NULL

  # Check for class of A
  if(!"dgCMatrix" %in% class(A)){
    A <- methods::as(A, "dgCMatrix")
  }
  
  # Check for self loops 
  if(!all(diag(A) == 0)){
    stop("Self-loop(s) detected (i.e., diagonal element(s) of A are not 0)")
  }
  
  # if A is named either col or row use that as id, else give id as row number
  if ( ( length(unique(rownames(A))) == nrow(A) ) | ( length(unique(colnames(A))) == nrow(A) )){
    ids <- if(length(unique(rownames(A))) != nrow(A)){ colnames(A) } else { rownames(A) }
  } else {
    message(paste0("Rownames/colnames for A missing or not unique, generating node IDs 1:", nrow(A)))
    ids <- 1:nrow(A)
  }
  
  # Check for isolates and remove
  isolates <- which(rowSums(A)==0 & colSums(A)==0)
  
  if (length(isolates) > 0){
    A <- A[-isolates, -isolates]
    message(paste0(length(isolates), " isolate(s) removed. Specifically node(s): ", paste(ids[isolates], collapse = ", ")))
    ids <- ids[-isolates]
  }
  
  # Check model 
  if(!model %in% c("NDH", "RS", "RSR")){
    stop("Model needs to be one of the following: NDH, RS, or RSR")
  } else{
    cl$model <- eval(model)
  }
  
  # If unsymmetric A provided for model = "NDH" or "RS" convert to symmetric A and warn
  if(!isSymmetric(A) & (model %in% c("NDH", "RS"))){
    A <- 1.0 * ( (A + t(A)) > 0.0)
    message(paste0("Unsymmetric A matrix supplied for model = ", model, ", converting to symmetric matrix"))
  }
  
  # Check initialization
  if(!is.list(initialization) && (!initialization %in% c("random", "GNN"))){
    stop("Please provide one for the following for initialization: 'random', 'GNN', or a named list with the necessary starting paramters")
  } else {
    cl$initialization <- eval(initialization)
  }
  
  # Check DA_type
  if(!DA_type %in% c("none", "cooling", "heating", "hybrid")){
    stop("Please provide one of the following for IC_slection: 'none', 'cooling', 'heating', or 'hybrid'")
  }
  
  # update con n_control if case_control is T
  if(case_control == T){
    con$n_control <- 100
  } else {
    control$n_control <- NULL
  }
  
  # update con beta_temp_schedule by DA_type argument
  if(DA_type == "cooling"){
    con$beta_temp_schedule <- seq(0.5, 1, length.out = 10)
  } else if (DA_type == "heating"){
    con$beta_temp_schedule <- seq(1.5, 1, length.out = 10)
  } else if (DA_type == "hybrid"){
    con$beta_temp_schedule <- c(seq(0.5, 1.5, length.out = 5), seq(1.4, 1, length.out = 5))
  } else {
    control$beta_temp_schedule <- 1
  }
  
  # updated con n_start if user supplied starting values
  if (is.list(initialization)){
    control$n_start <- 1
  }
  
  # check if names of elements in list match control_default
  nmsC <- names(con)
  namc <- names(control)
  
  noNms <- namc[!namc %in% nmsC]
  if (length(noNms) != 0) {
    stop("unknown names in control: ", paste(noNms, collapse = ", "))
  }
  con[namc] <- control
  
  # Check value of IC_selection
  if(!con[["IC_selection"]] %in% c("BIC_logit", "BIC_mbc", "ICL_mbc", "Total_BIC", "Total_ICL")) {
    stop("Please provide one of the following for IC_slection: 'BIC_logit', 'BIC_mbc', 'ICL_mbc', 'Total_BIC', or 'Total_ICL'")
  } 
  
  # Check value of n_control supplied and compute n_control
  if(!is.null(con[["n_control"]]) && !(con[["n_control"]]< max(nrow(A)-rowSums(A)) & con[["n_control"]]>0)){
    stop("Please supply a n_control value in (0, max(nrow(A)-rowSums(A)))")
  } 
  
  # check value of quantile_diff
  if(!(con[["quantile_diff"]]>=0 & con[["quantile_diff"]]<=1)){
    stop("Please supply a quantile_diff in [0,1]")
  }
  
  # Check termination rule supplied and create array for results
  if (!con[["termination_rule"]] %in% c("ARI", "NMI", "CER", "prob_mat","Q")){
    stop("Please provide one of the following termination rules: 'ARI', 'NMI', 'CER', 'prob_mat', or 'Q'")
  }
  
  # Check n_its_start_CA
  if (con$n_its_start_CA < 1){
    stop("Please supply a n_its_start_CA value >= 1")
  }
  
  cl$control <- eval(con)
  
  # Check initialization list if supplied
  if(is.list(initialization)){
    
    K <- length(initialization$p)
    D <- ncol(initialization$U)
    cl$K <- K
    cl$D <- D
    
    check_initial_values(list_name = initialization,
                         A = A,
                         K = K,
                         D = D,
                         n_interior_knots = con$n_interior_knots,
                         model = model)
  }
  
  # Check prior if supplied
  if(!is.null(con$priors)){
    check_priors(priors = con$priors,
                 D = D,
                 K = K,
                 n_interior_knots = con$n_interior_knots,
                 model = model)
  }
  
  combinations_2run <- as.matrix(expand.grid(K = unique(K), D = unique(D), n_start = 1:con$n_start))
  combinations_2run <- combinations_2run[order(combinations_2run[, "K"],
                                               combinations_2run[, "D"], 
                                               combinations_2run[, "n_start"]), , drop = F]
  combinations_2run <- cbind(combinations_2run, seed * combinations_2run[, "K"] * combinations_2run[, "D"] * combinations_2run[, "n_start"])
  rownames(combinations_2run) <- NULL
  colnames(combinations_2run)[ncol(combinations_2run)] <- "seed"
  cl$combinations_2run <- combinations_2run
  
  if(nrow(combinations_2run) > 1){
  
    progressr::handlers(progressr::handler_progress(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                                    complete = "=",   # Completion bar character
                                                    incomplete = "-", # Incomplete bar character
                                                    current = ">",    # Current bar character
                                                    clear = FALSE,    # If TRUE, clears the bar when finish
                                                    width = 100))
    progressr::with_progress({
      p <- progressr::progressor(steps = nrow(combinations_2run))
      parallel_res <- future.apply::future_lapply(X = 1:nrow(combinations_2run), 
                                                  FUN = function(x){
                                                    out <- inner_parallel(x = x,
                                                                          call_def = cl,
                                                                          A = A)
                                                    p()
                                                    return(out)
                                                  },
                                                  future.globals = FALSE,
                                                  future.packages = "JANE",
                                                  future.seed = seed)
    })
    
  } else {
    parallel_res <- list(inner_parallel(x = 1, call_def = cl, A = A))
  }
  
  IC_out <- do.call("rbind", lapply(parallel_res, function(x){x$EM_results$IC}))
  IC_out <- cbind(combinations_2run, IC_out)
  
  selected <- rep(0, nrow(IC_out))
  optimal_pos <- which(IC_out[, con$IC_selection] == min(IC_out[, con$IC_selection]))
  selected[optimal_pos] <- 1
  IC_out <- cbind(IC_out, selected)
  optimal_pos <- ifelse(length(optimal_pos) > 1, optimal_pos[1], optimal_pos)
  
  optimal_res <- parallel_res[[optimal_pos]]$EM_results
  if (length(optimal_res) > 1){
    optimal_res[["cluster_labels"]] <- apply(optimal_res$prob_mat, 1, which.max)
    names(optimal_res[["cluster_labels"]]) <- ids
    rownames(optimal_res$prob_matrix) <- ids
    rownames(optimal_res$U) <- ids
  } else {
    optimal_res <- NULL
  }
  
  optimal_starting <- parallel_res[[optimal_pos]]$starting_params
  if(!is.null(optimal_starting) & is.list(optimal_starting)){
    optimal_starting[["model"]] <- model
    optimal_starting[["cluster_labels"]] <- apply(optimal_starting$prob_mat, 1, which.max)
    names(optimal_starting[["cluster_labels"]]) <- ids
    rownames(optimal_starting$prob_matrix) <- ids
    rownames(optimal_starting$U) <- ids
  } else {
    optimal_starting <- NULL
  }
  
  if(is.list(initialization)){
    IC_out <- IC_out[, !colnames(IC_out) %in% c("n_start", "seed", "selected")]
  }
  return(list(optimal_res = optimal_res[sort(names(optimal_res))],
              optimal_starting = optimal_starting[sort(names(optimal_starting))],
              IC_out = IC_out))
  
}

#' @export
EM_inner <- function(A,
                     D,
                     K,
                     model,
                     starting_params,
                     control,
                     ...){
  
  extra_args <-  list(...)
  
  # Run initialize function
  current <- initialize_fun(A = A,
                            list_name = starting_params, 
                            model = model, 
                            n_interior_knots = control$n_interior_knots,
                            n_control = control$n_control, 
                            priors = control$priors, 
                            K = K,
                            D = D)
  
  current$termination_rule <- control$termination_rule
  
  current$termination_metric <- array(NA, dim = c(control$max_its, 
                                                  ncol = 6 + (!(control$termination_rule %in% c("ARI", "NMI", "CER", "Q"))) - (control$termination_rule == "Q"),
                                                  length(control$beta_temp_schedule)))
  
  colnames(current$termination_metric) <- if(ncol(current$termination_metric) == 7){
    c("ARI", "NMI", "CER", paste0(control$termination_rule, ifelse(control$termination_rule == "prob_mat", 
                                                                   paste0("_*",control$quantile_diff),
                                                                   "")), 
      "abs_diff_U", "abs_diff_MA_metric", "abs_diff_MA_U")
    
  } else if(ncol(current$termination_metric) == 5){
    
    c("ARI", "NMI", "CER", control$termination_rule, "abs_diff_MA_metric")
    
  } else { 
    
    c("ARI", "NMI", "CER", "abs_diff_U", "abs_diff_MA", "abs_diff_MA_U")
  }
  
  
  col_term_metric <- grep(control$termination_rule, 
                          colnames(current$termination_metric))
  
  if (control$termination_rule != "Q"){
    col_term_U <- grep("abs_diff_U", 
                       colnames(current$termination_metric))
  }
  
  current$convergence_ind <- matrix(0, nrow = length(control$beta_temp_schedule), ncol = 3)
  colnames(current$convergence_ind) <- c("beta_temperature", "convergence_ind", "n_iterations")
  current$convergence_ind[, "beta_temperature"] <- control$beta_temp_schedule
  current$convergence_ind[, "n_iterations"] <- rep(control$max_its, time = length(control$beta_temp_schedule))
  
  if(nrow(extra_args$combinations_2run) == 1){
    pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = control$max_its*length(control$beta_temp_schedule),
                           complete = "=",   # Completion bar character
                           incomplete = "-", # Incomplete bar character
                           current = ">",    # Current bar character
                           clear = FALSE,    # If TRUE, clears the bar when finish
                           width = 100)      # Width of the progress bar
  }
  
  # Start loop
  for (beta_temp in 1:length(control$beta_temp_schedule)){
    
    for(n_its in 1:control$max_its){
      
      if(nrow(extra_args$combinations_2run) == 1){
        pb$tick()
      }
      
      # update U
      current$fun_list$update_U(U = current$U, 
                                A = A, 
                                n_control = control$n_control,
                                mus = current$mus, 
                                omegas = current$omegas, 
                                prob_matrix = current$prob_matrix, 
                                beta = current$beta,
                                X =  current$X,
                                model = current$model)
      
      # update prob_matrix
      current$fun_list$update_prob_matrix(prob_matrix = current$prob_matrix, 
                                          mus = current$mus, omegas = current$omegas, 
                                          p = current$p, U = current$U,
                                          temp_beta = as.double(control$beta_temp_schedule[beta_temp]))
      
      # update p
      current$fun_list$update_p(prob_matrix = current$prob_matrix, p = current$p, nu = current$priors$nu)
      
      
      # update mus and omegas
      current$fun_list$update_mus_omegas(prob_matrix = current$prob_matrix,
                                         U = current$U, b = current$priors$b, a = current$priors$a,
                                         c = current$priors$c, G = current$priors$G,
                                         mus = current$mus, omegas = current$omegas)
      
      # update beta
      current$fun_list$update_beta(A = A, 
                                   n_control = control$n_control,
                                   U = current$U,
                                   beta = current$beta, 
                                   f = current$priors$f, 
                                   e = current$priors$e,
                                   X =  current$X,
                                   model = current$model)
      
      # check termination
      check_convergence <- terminate_EM(A = A,
                                        current = current,  
                                        termination_rule = control$termination_rule,
                                        tolerance = control$tolerance,
                                        tolerance_ARI = control$tolerance_ARI,
                                        tolerance_NMI = control$tolerance_NMI,
                                        tolerance_CER = control$tolerance_CER,
                                        quantile_diff = control$quantile_diff,
                                        n_control = control$n_control)
      
      
      # Additional check for convergence using moving average (MA)
      # Here we check if the MA of the termination metric is changing by iterations. If there
      # has been at least x (e.g. 10) number of consecutive iterations where 
      # the abs diff in MA is less than control$tolerance_diff_CA.
      
      if(n_its <= control$n_its_start_CA){
        
        MA <- ifelse(is.infinite(check_convergence$metric[col_term_metric]), 0,
                     check_convergence$metric[col_term_metric])
        diff_MA <- Inf
        counter <- 0
        
        if(control$termination_rule != "Q"){
          MA_U <- ifelse(is.infinite(check_convergence$metric[col_term_U]), 0,
                         check_convergence$metric[col_term_U])
          diff_MA_U <- Inf
          counter_U <- 0
        }
        
      } else {
        
        diff_MA_prev <- diff_MA
        diff_MA <- abs( ((MA - check_convergence$metric[col_term_metric]))/(n_its) )
        counter <- (counter + 1) * (diff_MA<control$tolerance_diff_CA)*(diff_MA_prev<control$tolerance_diff_CA)
        MA <- (check_convergence$metric[col_term_metric] + (n_its-1)*MA)/(n_its)
        
        if(control$termination_rule != "Q"){
          diff_MA_U_prev <- diff_MA_U
          diff_MA_U <- abs( ((MA_U - check_convergence$metric[col_term_U]))/(n_its) )
          counter_U <- (counter_U + 1) * (diff_MA_U<control$tolerance_diff_CA)*(diff_MA_U_prev<control$tolerance_diff_CA)
          MA_U <- (check_convergence$metric[col_term_U] + (n_its-1)*MA_U)/(n_its)
        }
        
      }
      
      if(control$termination_rule != "Q"){
        
        current$termination_metric[n_its,,beta_temp] <- c(check_convergence$metric, diff_MA, diff_MA_U)
        
        if((check_convergence$terminate | (counter >= control$consecutive_diff_CA & counter_U >= control$consecutive_diff_CA)) & (n_its > control$min_its)){
          current$convergence_ind[beta_temp, "convergence_ind"] <- 1L
          current$convergence_ind[beta_temp, "n_iterations"] <- n_its
          break
        }
        
      } else {
        
        current$termination_metric[n_its,,beta_temp] <- c(check_convergence$metric, diff_MA)
        
        if((check_convergence$terminate | (counter >= control$consecutive_diff_CA)) & (n_its > control$min_its)){
          current$convergence_ind[beta_temp, "convergence_ind"] <- 1L
          current$convergence_ind[beta_temp, "n_iterations"] <- n_its
          break
        }
        
      }
      
    }
    
  }
  
  current$termination_metric <- do.call("rbind", apply(current$termination_metric, 3, 
                                                       function(x){x[apply(x,1, function(x){!any(is.na(x) & !is.nan(x))}),]},
                                                       simplify = F))
  
  dont_return_names <- c("log_Q", "previous_prob_mat", "previous_U", "fun_list")
  out <- mget(x = names(current)[!(names(current) %in% dont_return_names)], envir = current)
  
  # Get IC info
  out$IC <- unlist(BICL(A = A, object = out))
  
  return(out[!(names(out) %in% "X")])
  
}

#' @export
inner_parallel <- function(x, call_def, A){
  
  combinations_2run_x <- call_def$combinations_2run[x, ]
  call_def$K <- combinations_2run_x["K"]
  call_def$D <- combinations_2run_x["D"]
  
  set.seed(combinations_2run_x["seed"]) #set.seed
  
  retry <- T
  retry_counter <- 0
  
  while(retry & retry_counter <= call_def$control$max_retry){
    
    if (!is.list(call_def$initialization)){
      
      call_def[[1]] <- as.symbol("initialize_starting_values")
      call_def$random_start <- ifelse(call_def$initialization == "GNN", F, T)
      call_def$starting_params <- eval(call_def)
      
    } else{
      call_def$starting_params <- call_def$initialization
    }
    
    run_fun <- tryCatch(
      {
        call_def[[1]] <- as.symbol("EM_inner")
        eval(call_def)
      },
      
      error = function(e) {
        message("Issues with starting values. Retrying with new starting values.\n")
        NA
      },
      warning = function(w) {
        message("Issues with starting values. Retrying with new starting values.\n")
        NA
      }
    ) 
    
    if (length(run_fun)>1){
      
      return(list(EM_results = run_fun,
                  starting_params = call_def$starting_params))
      
    } else {
      
      retry_counter <- retry_counter + 1
      
    }
    
  }
  
  if(retry){
    
    warning("Max re-try (i.e., max_retry) attempts reached. Issues with starting values. Returning Inf values. If this occurs often consider using alternative initialization.")
    
    return(list(EM_results = list(IC = c(BIC_logit = Inf,
                                         BIC_mbc = Inf,
                                         ICL_mbc = Inf,
                                         Total_BIC = Inf,
                                         Total_ICL = Inf)),
                starting_params = Inf))
    
  }
  
}

#' @useDynLib JANE  
update_U <- function(U, A, mus, omegas, prob_matrix, beta, X, n_control, model) {
  invisible(.Call('_JANE_update_U', PACKAGE = 'JANE', U, A, mus, omegas, prob_matrix, beta, X, n_control, model))
}

#' @useDynLib JANE  
update_U_CC <- function(U, n_control, A, mus, omegas, prob_matrix, beta, X, model) {
  invisible(.Call('_JANE_update_U_CC', PACKAGE = 'JANE', U, n_control, A, mus, omegas, prob_matrix, beta, X, model))
}

#' @useDynLib JANE  
update_U_RE <- function(U, A, mus, omegas, prob_matrix, beta, X, model, n_control) {
  invisible(.Call('_JANE_update_U_RE', PACKAGE = 'JANE', U, A, mus, omegas, prob_matrix, beta, X, model, n_control))
}

#' @useDynLib JANE  
update_U_RE_CC <- function(U, n_control, A, mus, omegas, prob_matrix, beta, X, model) {
  invisible(.Call('_JANE_update_U_RE_CC', PACKAGE = 'JANE', U, n_control, A, mus, omegas, prob_matrix, beta, X, model))
}

#' @useDynLib JANE  
update_beta <- function(beta, A, U, f, e, X, n_control, model) {
  invisible(.Call('_JANE_update_beta', PACKAGE = 'JANE', beta, A, U, f, e, X, n_control, model))
}

#' @useDynLib JANE  
update_beta_CC <- function(beta, A, n_control, U, f, e, X, model) {
  invisible(.Call('_JANE_update_beta_CC', PACKAGE = 'JANE', beta, A, n_control, U, f, e, X, model))
}

#' @useDynLib JANE  
update_beta_RE <- function(beta, A, U, f, e, X, model, n_control) {
  invisible(.Call('_JANE_update_beta_RE', PACKAGE = 'JANE', beta, A, U, f, e, X, model, n_control))
}

#' @useDynLib JANE  
update_beta_RE_CC <- function(beta, A, n_control, U, f, e, X, model) {
  invisible(.Call('_JANE_update_beta_RE_CC', PACKAGE = 'JANE', beta, A, n_control, U, f, e, X, model))
}

#' @useDynLib JANE  
update_mus_omegas <- function(prob_matrix, U, b, a, c, G, mus, omegas) {
  invisible(.Call('_JANE_update_mus_omegas', PACKAGE = 'JANE', prob_matrix, U, b, a, c, G, mus, omegas))
}

#' @useDynLib JANE  
update_p <- function(prob_matrix, p, nu) {
  invisible(.Call('_JANE_update_p', PACKAGE = 'JANE', prob_matrix, p, nu))
}

#' @useDynLib JANE  
update_prob_matrix_DA <- function(prob_matrix, mus, omegas, p, U, temp_beta) {
  invisible(.Call('_JANE_update_prob_matrix_DA', PACKAGE = 'JANE', prob_matrix, mus, omegas, p, U, temp_beta))
}

