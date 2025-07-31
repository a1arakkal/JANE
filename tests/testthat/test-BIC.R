
## Testing sim_A

test_that("BIC works", {
  
  library(Matrix)
  set.seed(2)
  
  sim_configs <- list(
    
    spherical_ws = list(mus = matrix(c(-5,-5, 5,-5, 5,5), nrow = 3, ncol = 2, byrow = T)/5,  
                        omegas = array(c(diag(rep(7, 2)), diag(rep(7, 2)), diag(rep(7, 2))), dim = c(2,2,3))),
    
    spherical_sws = list(mus = matrix(c(-5,-5, 5,-5, 5,5), nrow = 3, ncol = 2, byrow = T)/6,
                         omegas = array(c(diag(rep(7, 2)), diag(rep(7, 2)), diag(rep(7, 2))), dim = c(2,2,3))),
    
    ellipsoidal_ws = list(mus = rbind(c(2.0,0),
                                      c(-1,1.15),
                                      c(-1,-1.15)),
                          omegas = array(c(c(8,0,0,1), c(1,0,0,4), c(1,0,0,4)), dim = c(2,2,3))),
    
    ellipsoidal_sws = list(mus = rbind(c(1.0,0),
                                       c(-1,1.15),
                                       c(-1,-1.15)),
                           omegas = array(c(c(8,0,0,1), c(1,0,0,4), c(1,0,0,4)), dim = c(2,2,3)))
    
    
  )
  
  mus <- sim_configs$ellipsoidal_sws$mus
  omegas <- sim_configs$ellipsoidal_sws$omegas
  beta0 <- -1
  
  # set.seed(2024)
  p <- extraDistr::rdirichlet(1, alpha = rep(3, 3))[1,] # true mixture probabilities
  
  temp <- expand.grid(c("NDH", "RS", "RSR"), c(F,T), c("bernoulli", "lognormal",  "poisson"))
  temp$Var4 <- ifelse(temp$Var2, 0.1, 0.0)
  temp$Var1 <- as.character(temp$Var1)
  temp$Var3 <- as.character(temp$Var3)
  
  res <- logical(nrow(temp))
  for (m in 1:nrow(temp)) {
    
    family <- temp$Var3[m]
    model <- temp$Var1[m]
    noise_weights <- temp$Var2[m]
    q_prob <- temp$Var4[m]
    
    if(family == "lognormal"){
      params_weights <-  list(beta0 = 2,
                              precision_weights = 1)
      precision_noise_weights <- 2
      mean_noise_weights <- 2
    } else if (family == "poisson"){
      params_weights <-  list(beta0 = log(5))
      precision_noise_weights <- NULL
      mean_noise_weights <- 2
    } else {
      params_weights <- NULL
      precision_noise_weights <- NULL
      mean_noise_weights <- NULL
    }
    
    if (model == "RSR"){
      params_weights$precision_R_effects <- solve(matrix(c(0.5,0.1,0.1,0.5), nrow = 2))
    } else if (model == "RS"){
      params_weights$precision_R_effects <- 2
    } else {
      params_weights$precision_R_effects <- NULL
    }
    
    sim_data <- JANE::sim_A(N = 100,
                            model = model,
                            family = family,
                            mus = mus,
                            omegas = omegas,
                            p = p,
                            params_LR = list(beta0 = beta0),
                            params_weights = params_weights,
                            precision_noise_weights = precision_noise_weights,
                            mean_noise_weights = mean_noise_weights,
                            noise_weights_prob = q_prob,
                            remove_isolates = TRUE)
    
    control <- list(
      verbose = FALSE,
      max_its = 1e3, 
      min_its = 10,  
      priors = NULL,  
      n_interior_knots = 5, 
      termination_rule = "prob_mat",
      tolerance = 1e-3, 
      tolerance_ARI = 0.999, 
      tolerance_NMI = 0.999, 
      tolerance_CER = 0.01,  
      n_its_start_CA = 20, 
      tolerance_diff_CA = 1e-3, 
      consecutive_diff_CA = 5,
      quantile_diff = 1,
      beta_temp_schedule = 1, 
      n_control = NULL, 
      n_start = 5,
      max_retry = 5, 
      IC_selection = "Total_ICL", 
      sd_random_U_GNN = 1, 
      max_retry_GNN = 10, 
      n_its_GNN = 10,
      downsampling_GNN = TRUE 
    )
    
    starting_params <- JANE:::initialize_starting_values(A = (sim_data$W>0)*1, 
                                                         K = 3,
                                                         D = 2, 
                                                         model = model,
                                                         random_start = F,
                                                         control = control,
                                                         family = family,
                                                         noise_weights = noise_weights,
                                                         guess_noise_weights = 1,
                                                         prob_matrix_W= prob_matrix_W,
                                                         q_prob = q_prob)
    
    A = sim_data$W
    
    # Extract edge indices and weights and store as prob_matrix_W
    prob_matrix_W <- cbind(as.matrix(summary(A)), 1.0, 0.0)
    colnames(prob_matrix_W) <- c("i", "j", "weight", "hat_zij1", "hat_zij2")  
    
    # Only store upper triangular for RS and NDH as it is symmetric
    if(model != "RSR"){
      prob_matrix_W <- prob_matrix_W[prob_matrix_W[,"j"]>prob_matrix_W[,"i"], ]
      
      if(noise_weights){
        A[prob_matrix_W[, c("i", "j")]] <- rbeta(nrow(prob_matrix_W), 1,1)
        A[prob_matrix_W[, c("j", "i")]] <-  A[prob_matrix_W[, c("i", "j")]]
      }
    } else {
      if(noise_weights){
        A[prob_matrix_W[, c("i", "j")]] <- rbeta(nrow(prob_matrix_W), 1,1)
      }
    }
    
    object <- JANE:::initialize_fun(A = (sim_data$W>0)*1, 
                                    K = 3,
                                    D = 2,
                                    family = family,
                                    noise_weights = noise_weights,
                                    prob_matrix_W= prob_matrix_W,
                                    priors = NULL, 
                                    list_name = starting_params, 
                                    model = model, 
                                    n_interior_knots = control$n_interior_knots,
                                    n_control = control$n_control)
    
    object <- mget(x = names(object), envir = object)
    
    #### NDH -----------------------------------------------------------------------
    
    # noise_weights = F
    
    bic_logit_NDH_R <- function(A, object){
      
      N <- nrow(A)
      
      if(!object$noise_weights){
        
        n_edges <- 0.5*sum(A)
        p1 <- 0
        for (i in 1:N){
          for (j in 1:N){
            
            if (j != i){
              eta <- object$beta - crossprod(object$U[i,] - object$U[j,])[1,1]
              p1 <- p1 + ( ( A[i,j] * eta ) - log( 1 + exp(eta)) )
            }
            
          }
        }
        
        return(-2*(0.5*p1) + (length(object$beta)*log(n_edges)))
        
      } else {
        
        A[object$prob_matrix_W[, 1:2]] <- object$prob_matrix_W[, 3]
        A[object$prob_matrix_W[, 2:1]] <- object$prob_matrix_W[, 3]
        q_prob <- object$q_prob
        n_edges <- 0.5*sum(A>0)
        
        if (family == "bernoulli"){
          p1 <- 0
          for (i in 1:N){
            for (j in 1:N){
              
              if (j != i){
                eta <- object$beta - crossprod(object$U[i,] - object$U[j,])[1,1]
                p1 <- p1 + ( (1-(A[i,j]>0))*log(1-q_prob) - log(1+exp(eta)) + (A[i,j]>0)*log(exp(eta) + q_prob) )
              }
              
            }
          }
          n_params <- length(object$beta) + 1
        } else if (family == "poisson"){
          
          p1 <- 0
          for (i in 1:N){
            for (j in 1:N){
              
              if (j != i){
                eta <- object$beta - crossprod(object$U[i,] - object$U[j,])[1,1]
                
                eta_w <- object$beta2
                w <- A[i,j]
                density_w <- extraDistr::dtpois(x = w, lambda = exp(eta_w), a = 0.0, log = F)
                density_noise<- extraDistr::dtpois(x = w, lambda = object$guess_noise_weights, a = 0.0, log =F)
                
                if(w<=0){
                  p1 <- p1 + ( log(1-q_prob) - log(1+exp(eta)) )
                } else {
                  p1 <- p1 + (-log(1+exp(eta)) + log(exp(eta)*density_w + q_prob*density_noise))
                }
              }
              
            }
          }
          n_params <- length(object$beta) +length(object$beta2)+1
        } else {
          
          p1 <- 0
          for (i in 1:N){
            for (j in 1:N){
              
              if (j != i){
                eta <- object$beta - crossprod(object$U[i,] - object$U[j,])[1,1]
                
                eta_w <- object$beta2
                w <- A[i,j]
                density_w <- dlnorm(x = w, meanlog = eta_w, sdlog = 1/(sqrt(object$precision_weights)), log = F)
                density_noise<- dlnorm(x = w, meanlog = object$guess_noise_weights, sdlog = 1/(sqrt(object$precision_noise_weights)), log = F)
                
                if(w<=0){
                  p1 <- p1 + ( log(1-q_prob) - log(1+exp(eta)) )
                } else {
                  p1 <- p1 + (-log(1+exp(eta)) + log(exp(eta)*density_w + q_prob*density_noise))
                }
              }
              
            }
          }
          n_params <- length(object$beta) +length(object$beta2)+2+1
        }
        
        
        return(-2*(0.5*p1) + (n_params)*log(n_edges))
        
        
      }
    }
    
    # all.equal(JANE:::BICL(A = A, object = object)$BIC_model, bic_logit_NDH_R(A, object))
    
    bic_logit_RS_R <- function(A, object){
      
      N <- nrow(A)
      
      
      if(!object$noise_weights){
        
        n_edges <- 0.5*sum(A)
        p1 <- 0
        
        for (i in 1:N){
          for (j in 1:N){
            
            if (j != i){
              x_ij <- c(1, object$X[i,] + object$X[j,])
              beta_x <- (t(x_ij) %*% object$beta)[1,1]
              eta <- beta_x - crossprod(object$U[i,] - object$U[j,])[1,1]
              p1 <- p1 + ( ( A[i,j] * eta ) - log( 1 + exp(eta)) )
            }
            
          }
        }
        
        return(-2*(0.5*p1) + (length(object$beta)*log(n_edges)))
        
      } else {
        A[object$prob_matrix_W[, 1:2]] <- object$prob_matrix_W[, 3]
        A[object$prob_matrix_W[, 2:1]] <- object$prob_matrix_W[, 3]
        q_prob <- object$q_prob
        n_edges <- 0.5*sum(A>0)
        
        if (family == "bernoulli"){
          p1 <- 0
          for (i in 1:N){
            for (j in 1:N){
              
              if (j != i){
                x_ij <- c(1, object$X[i,] + object$X[j,])
                beta_x <- (t(x_ij) %*% object$beta)[1,1]
                eta <- beta_x - crossprod(object$U[i,] - object$U[j,])[1,1]
                p1 <- p1 + ( (1-(A[i,j]>0))*log(1-q_prob) - log(1+exp(eta)) + (A[i,j]>0)*log(exp(eta) + q_prob) )
              }
              
            }
          }
          n_params <- length(object$beta) + 1
        } else if (family == "poisson"){
          
          p1 <- 0
          for (i in 1:N){
            for (j in 1:N){
              
              if (j != i){
                x_ij <- c(1, object$X[i,] + object$X[j,])
                beta_x <- (t(x_ij) %*% object$beta)[1,1]
                eta <- beta_x - crossprod(object$U[i,] - object$U[j,])[1,1]
                
                x_ij <- c(1, object$X2[i,] + object$X2[j,])
                beta_x <- (t(x_ij) %*% object$beta2)[1,1]
                eta_w <- beta_x
                
                w <- A[i,j]
                density_w <- extraDistr::dtpois(x = w, lambda = exp(eta_w), a = 0.0, log = F)
                density_noise<- extraDistr::dtpois(x = w, lambda = object$guess_noise_weights, a = 0.0, log =F)
                
                if(w<=0){
                  p1 <- p1 + ( log(1-q_prob) - log(1+exp(eta)) )
                } else {
                  p1 <- p1 + (-log(1+exp(eta)) + log(exp(eta)*density_w + q_prob*density_noise))
                }
              }
              
            }
          }
          n_params <- length(object$beta) +length(object$beta2)+1
        } else {
          
          p1 <- 0
          for (i in 1:N){
            for (j in 1:N){
              
              if (j != i){
                
                x_ij <- c(1, object$X[i,] + object$X[j,])
                beta_x <- (t(x_ij) %*% object$beta)[1,1]
                eta <- beta_x - crossprod(object$U[i,] - object$U[j,])[1,1]
                
                x_ij <- c(1, object$X2[i,] + object$X2[j,])
                beta_x <- (t(x_ij) %*% object$beta2)[1,1]
                eta_w <- beta_x
                
                w <- A[i,j]
                density_w <- dlnorm(x = w, meanlog = eta_w, sdlog = 1/(sqrt(object$precision_weights)), log = F)
                density_noise<- dlnorm(x = w, meanlog = object$guess_noise_weights, sdlog = 1/(sqrt(object$precision_noise_weights)), log = F)
                
                if(w<=0){
                  p1 <- p1 + ( log(1-q_prob) - log(1+exp(eta)) )
                } else {
                  p1 <- p1 + (-log(1+exp(eta)) + log(exp(eta)*density_w + q_prob*density_noise))
                }
              }
              
            }
          }
          n_params <- length(object$beta) +length(object$beta2)+2+1
        }
        
        
      }
      return(-2*(0.5*p1) + ((n_params)*log(n_edges)))
    }
    
    # all.equal(JANE:::BICL(A = A, object = object)$BIC_model, bic_logit_RS_R(A, object))
    
    
    bic_logit_RSR_R <- function(A, object){
      
      N <- nrow(A)
      
      if(!object$noise_weights){
        
        n_edges <- sum(A)
        p1 <- 0
        
        for (i in 1:N){
          for (j in 1:N){
            
            if (j != i){
              x_ij <- c(1, object$X[i, 1:(ncol(object$X) * 0.5)], object$X[j, (ncol(object$X)*0.5 + 1):ncol(object$X) ])
              beta_x <- (t(x_ij) %*% object$beta)[1,1]
              eta <- beta_x - crossprod(object$U[i,] - object$U[j,])[1,1]
              p1 <- p1 + ( ( A[i,j] * eta ) - log( 1 + exp(eta)) )
            }
            
          }
        }
        
        return(-2*(p1) + (length(object$beta)*log(n_edges)))
        
      } else {
        A[object$prob_matrix_W[, 1:2]] <- object$prob_matrix_W[, 3]
        q_prob <- object$q_prob
        n_edges <- sum(A>0)
        
        if (family == "bernoulli"){
          p1 <- 0
          for (i in 1:N){
            for (j in 1:N){
              
              if (j != i){
                x_ij <- c(1, object$X[i, 1:(ncol(object$X) * 0.5)], object$X[j, (ncol(object$X)*0.5 + 1):ncol(object$X) ])
                beta_x <- (t(x_ij) %*% object$beta)[1,1]
                eta <- beta_x - crossprod(object$U[i,] - object$U[j,])[1,1]
                p1 <- p1 + ( (1-(A[i,j]>0))*log(1-q_prob) - log(1+exp(eta)) + (A[i,j]>0)*log(exp(eta) + q_prob) )
              }
              
            }
          }
          n_params <- length(object$beta) + 1
        } else if (family == "poisson"){
          
          p1 <- 0
          for (i in 1:N){
            for (j in 1:N){
              
              if (j != i){
                x_ij <- c(1, object$X[i, 1:(ncol(object$X) * 0.5)], object$X[j, (ncol(object$X)*0.5 + 1):ncol(object$X) ])
                beta_x <- (t(x_ij) %*% object$beta)[1,1]
                eta <- beta_x - crossprod(object$U[i,] - object$U[j,])[1,1]
                
                x_ij <- c(1, object$X2[i, 1:(ncol(object$X2) * 0.5)], object$X2[j, (ncol(object$X2)*0.5 + 1):ncol(object$X2) ])
                beta_x <- (t(x_ij) %*% object$beta2)[1,1]
                eta_w <- beta_x
                
                w <- A[i,j]
                density_w <- extraDistr::dtpois(x = w, lambda = exp(eta_w), a = 0.0, log = F)
                density_noise<- extraDistr::dtpois(x = w, lambda = object$guess_noise_weights, a = 0.0, log =F)
                
                if(w<=0){
                  p1 <- p1 + ( log(1-q_prob) - log(1+exp(eta)) )
                } else {
                  p1 <- p1 + (-log(1+exp(eta)) + log(exp(eta)*density_w + q_prob*density_noise))
                }
              }
              
            }
          }
          n_params <- length(object$beta) +length(object$beta2)+1
        } else {
          
          p1 <- 0
          for (i in 1:N){
            for (j in 1:N){
              
              if (j != i){
                
                x_ij <- c(1, object$X[i, 1:(ncol(object$X) * 0.5)], object$X[j, (ncol(object$X)*0.5 + 1):ncol(object$X) ])
                beta_x <- (t(x_ij) %*% object$beta)[1,1]
                eta <- beta_x - crossprod(object$U[i,] - object$U[j,])[1,1]
                
                x_ij <- c(1, object$X2[i, 1:(ncol(object$X2) * 0.5)], object$X2[j, (ncol(object$X2)*0.5 + 1):ncol(object$X2) ])
                beta_x <- (t(x_ij) %*% object$beta2)[1,1]
                eta_w <- beta_x
                
                w <- A[i,j]
                density_w <- dlnorm(x = w, meanlog = eta_w, sdlog = 1/(sqrt(object$precision_weights)), log = F)
                density_noise<- dlnorm(x = w, meanlog = object$guess_noise_weights, sdlog = 1/(sqrt(object$precision_noise_weights)), log = F)
                
                if(w<=0){
                  p1 <- p1 + ( log(1-q_prob) - log(1+exp(eta)) )
                } else {
                  p1 <- p1 + (-log(1+exp(eta)) + log(exp(eta)*density_w + q_prob*density_noise))
                }
              }
              
            }
          }
          n_params <- length(object$beta) +length(object$beta2)+2+1
        }
        
        
      }
      return(-2*(p1) + ((n_params)*log(n_edges)))
    }
    
    # all.equal(JANE:::BICL(A = A, object = object)$BIC_model, bic_logit_RSR_R(A, object))
    
    dmvnorm_chol_quad_log <- function(x, mean, omega) {
      d <- length(mean)
      
      diff <- x - mean
      
      # Quadratic form: diff' * omega * diff
      rss <- (t(diff) %*% omega %*% diff)[1,1]
      
      # Cholesky decomposition of omega (precision matrix)
      L <- chol(omega)
      
      # log determinant of covariance matrix = -2 * sum(log(diag(L)))
      log_det_sigma <- -2 * sum(log(diag(L)))
      
      # Log density
      log_density <- -0.5 * (d * log(2 * pi) + log_det_sigma + rss)
      
      return(log_density)
    }
    
    bic_icl_mbc_R <- function(object){
      
      mus <- object$mus
      omegas <- object$omegas
      p <- object$p
      U <- object$U
      prob_mat <- object$prob_matrix
      
      N <- nrow(U)
      D <- ncol(U)
      K <- length(p)
      
      Z_star <- matrix(0, nrow = N, ncol = K)
      Z_star[cbind( 1:N, apply(prob_mat, 1, which.max) )] <- 1
      
      # all.equal(apply(Z_star,1,which.max), apply(prob_mat, 1, which.max))
      
      n_params <- (K-1) + (K*D) + ( K * (0.5*(D*(D+1))) )
      
      p1 <- 0
      
      for (i in 1:N) {
        
        p2 <- 0

        for (k in 1:K){
          den <- dmvnorm_chol_quad_log(x = U[i,], mean = mus[k,], omega = omegas[,,k])
          p2 <- p2 + (p[k]*exp(den))
        }
        
        p1 <- p1 + log(p2)
        
      }
      
      BIC <- (-2*p1) + (n_params*log(N))
      
      exp_entropy <- -1.0*sum(Z_star*ifelse(prob_mat > 0, log(prob_mat), 0))
      
      ICL <- BIC + 2*exp_entropy
      
      return(list(BIC_mbc = BIC,
                  ICL_mbc = ICL))
    }
    
    
    # all.equal(JANE:::BICL(A = A, object = object)[c("BIC_mbc", "ICL_mbc")], bic_icl_mbc_R(object))
    
    BICL_R <- function(A, object){
      
      if(object$model == "NDH"){
        BIC_model <- bic_logit_NDH_R(A, object)
      } else if(object$model == "RS" ){
        BIC_model <- bic_logit_RS_R(A, object)
      } else {
        BIC_model <- bic_logit_RSR_R(A, object)
      }
      
      bic_icl_mbc <- bic_icl_mbc_R(object)
      
      out <- list(BIC_model = BIC_model,
                  BIC_mbc = bic_icl_mbc$BIC_mbc,
                  ICL_mbc = bic_icl_mbc$ICL_mbc,
                  Total_BIC = BIC_model + bic_icl_mbc$BIC_mbc,
                  Total_ICL = BIC_model + bic_icl_mbc$ICL_mbc)
      
      return(out)
    }
    
    if(any(is.na(unlist(JANE:::BICL(A, object))))){
      stop("NA values detected")
    }
    res[m] <- all.equal(BICL_R(A, object), JANE:::BICL(A, object))
    print(JANE:::BICL(A, object))
    
  }
  
  expect_true(all(res))
  
  
})

