
terminate_EM <- function(A, current, termination_rule = "ARI", 
                         tolerance = sqrt(.Machine$double.eps),
                         tolerance_ARI = 0.999,
                         tolerance_NMI = 0.999,
                         tolerance_CER = 0.05, 
                         quantile_diff = 1, 
                         n_control){
  
  #termination indicator
  terminate <- 0
  
  current_Z <- as.double(apply(current$prob_matrix, 1, which.max))
  previous_Z <- as.double(apply(current$previous_prob_mat, 1, which.max))
  
  current_prob_mat <- current$prob_matrix * 1.0
  previous_prob_mat <- current$previous_prob_mat * 1.0
  current$previous_prob_mat <- current_prob_mat * 1.0
  
  if(current$noise_weights){
    
    # Commented out as computing ARI, NMI and CER on density*N*(N-1) obs is slow
    # current_Z_W <- as.double(apply(current$prob_matrix_W, 1, which.max))
    # previous_Z_W <- as.double(apply(current$previous_prob_mat_W, 1, which.max))
    
    current_prob_mat_W <- current$prob_matrix_W * 1.0
    previous_prob_mat_W <- current$previous_prob_mat_W * 1.0
    current$previous_prob_mat_W <- current_prob_mat_W * 1.0
    
    # Commented out as computing ARI, NMI and CER on density*N*(N-1) obs is slow
    # ARI_metric_W <- mclust::adjustedRandIndex(current_Z_W, previous_Z_W) 
    # NMI_metric_W <- aricode::NMI(current_Z_W, previous_Z_W) 
    # CER_metric_W <- mclust::classError(current_Z_W, previous_Z_W)$errorRate 
    
  } else {
    
    ARI_metric_W <- Inf
    NMI_metric_W <- Inf
    CER_metric_W <- -Inf
    
  }
  
  ARI_metric <- mclust::adjustedRandIndex(current_Z, previous_Z) 
  NMI_metric <- aricode::NMI(current_Z, previous_Z) 
  CER_metric <- mclust::classError(current_Z, previous_Z)$errorRate 
  metric <- NA
  
  current_U <- current$U * 1.0
  previous_U <- current$previous_U * 1.0
  current$previous_U <- current_U * 1.0
  diff_U <- unname(stats::quantile(abs(current_U - previous_U), quantile_diff))
  
  if(termination_rule == "ARI"){
    
    if(all(c(ARI_metric, ARI_metric_W) > tolerance_ARI)){
      terminate <- 1
    }
    
  } else if (termination_rule == "NMI"){
    
    if(all(c(NMI_metric, NMI_metric_W) > tolerance_NMI)){
      terminate <- 1
    }
    
  } else if (termination_rule == "CER"){
    
    if(all(c(CER_metric, CER_metric_W) < tolerance_CER)){
      terminate <- 1
    }
    
  } else if (termination_rule == "prob_mat"){
    
    metric <- unname(stats::quantile(abs(current_prob_mat - previous_prob_mat), quantile_diff))
    
    if(current$noise_weights){
      metricW <- unname(stats::quantile(abs(current_prob_mat_W[, c("hat_zij1", "hat_zij2")] - previous_prob_mat_W[, c("hat_zij1", "hat_zij2")]), quantile_diff))
    } else {
      metricW <- -Inf
    }

    if(all(c(metric, metricW) < tolerance)){
      terminate <- 1
    }
    
  } else {
    
    current_log_Q <- current$fun_list$log_Q(A = A,
                                            U = current$U,
                                            mus = current$mus,
                                            omegas = current$omegas,
                                            prob_matrix = current$prob_matrix,
                                            beta = current$beta,
                                            X = current$X, 
                                            n_control = n_control,
                                            p = current$p,
                                            a = current$priors$a,
                                            b = current$priors$b,
                                            c = current$priors$c,
                                            G = current$priors$G,
                                            nu = current$priors$nu,
                                            e = current$priors$e,
                                            f = current$priors$f,
                                            model = current$model)
    
    previous_log_Q <- current$log_Q
    
    current$log_Q <- current_log_Q
    
    metric <- abs(current_log_Q - previous_log_Q)
    
    if(metric < tolerance){
      terminate <- 1
    }
    
  }
  
  if (termination_rule %in% c("ARI","NMI","CER")){
    
    if(!current$noise_weights){
      
      out_metric_table <- c(ARI_metric, NMI_metric, CER_metric, diff_U)
      
    } else {
      
      loc <- grep(termination_rule, c("ARI_metric_W", "NMI_metric_W", "CER_metric_W"))
      out_metric_table <- c(ARI_metric, NMI_metric, CER_metric,
                            c(ARI_metric_W, NMI_metric_W, CER_metric_W)[loc],
                            diff_U)
      
    }
    
  } else if(termination_rule == "prob_mat"){
    
    if(!current$noise_weights){
      
      out_metric_table <- c(ARI_metric, NMI_metric, CER_metric, metric, diff_U)
      
    } else {
      
      out_metric_table <- c(ARI_metric, NMI_metric, CER_metric, metric, metricW, diff_U)
      
    }
    
  } else {
    
    out_metric_table <- c(ARI_metric, NMI_metric, CER_metric, metric)
    
  }
  
  return(list(terminate = terminate,
              metric = out_metric_table))
  
}