## Testing JANE

test_that("JANE works", {
   
   # Simulate network
   mus <- matrix(c(-1,-1,1,-1,1,1), 
                 nrow = 3,
                 ncol = 2, 
                 byrow = TRUE)
   omegas <- array(c(diag(rep(7,2)),
                     diag(rep(7,2)), 
                     diag(rep(7,2))), 
                   dim = c(2,2,3))
   p <- rep(1/3, 3)
   beta0 <- 1.0
   sim_data <- JANE::sim_A(N = 100L, 
                           model = "NDH",
                           mus = mus, 
                           omegas = omegas, 
                           p = p, 
                           params_LR = list(beta0 = beta0),
                           remove_isolates = TRUE)
   
   # Run JANE on simulated data
   expect_no_error( JANE::JANE(A = sim_data$A,
                               D = 2L,
                               K = 3L,
                               initialization = "GNN", 
                               model = "NDH",
                               case_control = FALSE,
                               DA_type = "none") )
   
   # Run JANE on simulated data - consider multiple D and K
   expect_no_error( JANE::JANE(A = sim_data$A,
                               D = 2:3,
                               K = 2:4,
                               initialization = "GNN", 
                               model = "NDH",
                               case_control = FALSE,
                               DA_type = "none") )
   
   # Run JANE on simulated data - parallel with 2 cores
   # expect_no_error({
   #   future::plan(future::multisession, workers = 2)
   #   JANE::JANE(A = sim_data$A,
   #                   D = 2L,
   #                   K = 3L,
   #                   initialization = "GNN",
   #                   model = "NDH",
   #                   case_control = FALSE,
   #                   DA_type = "none")
   #  future::plan(future::sequential)
   # })
   
   # Run JANE on simulated data - case/control approach with 20 controls sampled for each actor
   expect_no_error( JANE::JANE(A = sim_data$A,
                               D = 2L,
                               K = 3L,
                               initialization = "GNN", 
                               model = "NDH",
                               case_control = TRUE,
                               DA_type = "none",
                               control = list(n_control = 20)) )
   
   expect_no_error( JANE::JANE(A = sim_data$A,
                               D = 2L,
                               K = 3L,
                               initialization = "GNN", 
                               model = "RS",
                               case_control = TRUE,
                               DA_type = "none",
                               control = list(n_control = 20)) )
   
   expect_no_error( JANE::JANE(A = sim_data$A,
                               D = 2L,
                               K = 3L,
                               initialization = "GNN", 
                               model = "RSR",
                               case_control = TRUE,
                               DA_type = "none",
                               control = list(n_control = 20)) )
   # Reproducibility
   expect_equal( JANE::JANE(A = sim_data$A,
                            D = 2L,
                            K = 3L,
                            initialization = "GNN", 
                            seed = 1234,
                            model = "NDH",
                            case_control = FALSE,
                            DA_type = "none"),
                 
                 JANE::JANE(A = sim_data$A,
                            D = 2L,
                            K = 3L,
                            initialization = "GNN", 
                            seed = 1234,
                            model = "NDH",
                            case_control = FALSE,
                            DA_type = "none") )
   
   # Another reproducibility example where the seed was not set. 
   # It is possible to replicate the results using the starting values due to 
   # the nature of EM algorithms
   res3 <- JANE::JANE(A = sim_data$A,
                      D = 2L,
                      K = 3L,
                      initialization = "GNN", 
                      model = "NDH",
                      case_control = FALSE,
                      DA_type = "none")
   ## Extract starting values               
   start_vals <- res3$optimal_start
   
   ## Run JANE using extracted starting values, no need to specify D and K 
   ## below as function will determine those values from start_vals
   res4 <- JANE::JANE(A = sim_data$A,
                      initialization = start_vals, 
                      model = "NDH",
                      case_control = FALSE,
                      DA_type = "none")
   
   ## Check if optimal_res are identical
   expect_equal(res3$optimal_res, res4$optimal_res)       
   
   # Check with weighted network
   sim_data <- sim_A(N = 100L,
                     model = "RS",
                     family = "poisson",
                     mus = mus,
                     omegas = omegas,
                     p = rep(1/3, 3),
                     params_LR = list(beta0 = 1),
                     params_weights = list(beta0 = 2),
                     noise_weights_prob = 0.1,
                     mean_noise_weights = 1,
                     remove_isolates = TRUE)
   
   expect_no_error(
      JANE(A = sim_data$W,
           D = 2,
           K = 5,
           model = "RS",
           noise_weights = TRUE,
           family = "poisson",
           guess_noise_weights = 1L)
   )
   
   expect_no_error(
      JANE(A = sim_data$W,
           D = 2,
           K = 5,
           model = "RS",
           noise_weights = TRUE,
           family = "lognormal",
           guess_noise_weights = -3.5) # log-scale mean
   )
   
   
   expect_no_error(
      JANE(A = sim_data$A,
           D = 2,
           K = 5,
           model = "RSR",
           noise_weights = TRUE,
           family = "bernoulli",
           guess_noise_weights = 0.2) # expected noise edge proportion
   )
   
   expect_no_error(
      JANE(A = sim_data$W,
           D = 2,
           K = 5,
           model = "RS",
           noise_weights = TRUE,
           family = "poisson",
           guess_noise_weights = 1L,
           case_control = TRUE,
           control = list(n_control = 20)) 
   )
   
   expect_no_error(
      JANE(A = sim_data$W,
           D = 2,
           K = 5,
           model = "RS",
           noise_weights = TRUE,
           family = "lognormal",
           guess_noise_weights = -3.5,
           case_control = TRUE,
           control = list(n_control = 20)) 
   )
   
   
   expect_no_error(
      JANE(A = sim_data$A,
           D = 2,
           K = 5,
           model = "RSR",
           noise_weights = TRUE,
           family = "bernoulli",
           guess_noise_weights = 0.2,
           case_control = TRUE,
           control = list(n_control = 20)) 
   )
   
})

test_that("JANE comprehensive works", {
   
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
   p <- extraDistr::rdirichlet(1, alpha = rep(3, 3))[1,]
   
   all_test <- list()
   
   for(noise in c(0, 0.1)){
      for (model in c("NDH", "RS", "RSR")){
         for (family in c("bernoulli", "poisson", "lognormal")){
            
            if(family == "lognormal"){
               params_weights <-  list(beta0 = 3,
                                       precision_weights = 1)
               precision_noise_weights <- 2
               mean_noise_weights <- 0.1
            } else if (family == "poisson"){
               params_weights <-  list(beta0 = 3)
               precision_noise_weights <- NULL
               mean_noise_weights <- 2.5
            } else {
               params_weights <- NULL
               precision_noise_weights <- NULL
               mean_noise_weights <- NULL
            }
            
            if (model == "RSR"){
               # params_weights$precision_R_effects <- solve(matrix(c(0.5,0.1,0.1,0.5), nrow = 2))
               params_weights$precision_R_effects <- solve(matrix(c(0.1,0.05,0.05,0.1), nrow = 2))
            } else if (model == "RS"){
               params_weights$precision_R_effects <- 1/0.1
            } else {
               params_weights$precision_R_effects <- NULL
            }
            
            if (model == "RSR"){
               precision_R_effects <- solve(matrix(c(0.5,0.1,0.1,0.5), nrow = 2))
            } else if (model == "RS"){
               precision_R_effects <- 2
            } else {
               precision_R_effects <- NULL
            }
            
            params_LR <- list(beta0 = beta0,
                              precision_R_effects=precision_R_effects)
            
            sim_data <- JANE::sim_A(N = 100,
                                    model = model,
                                    family = family,
                                    mus = mus,
                                    omegas = omegas,
                                    p = p,
                                    params_LR = params_LR,
                                    params_weights = params_weights,
                                    precision_noise_weights = precision_noise_weights,
                                    mean_noise_weights = mean_noise_weights,
                                    noise_weights_prob = noise,
                                    remove_isolates = TRUE)
            
            # Extract edge indices and weights and store as prob_matrix_W
            prob_matrix_W <- cbind(as.matrix(summary(sim_data$W)), 1.0, 0.0)
            colnames(prob_matrix_W) <- c("i", "j", "weight", "hat_zij1", "hat_zij2")  
            
            # Only store upper triangular for RS and NDH as it is symmetric
            if(model != "RSR"){
               
               prob_matrix_W <- prob_matrix_W[prob_matrix_W[,"j"]>prob_matrix_W[,"i"], ]
               
            }
            
            control <- list(
               verbose = F,
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
                                                                 noise_weights = ifelse(noise == 0, F, T),
                                                                 guess_noise_weights = 1,
                                                                 prob_matrix_W= prob_matrix_W,
                                                                 q_prob = noise)
            
            
            
            current <- JANE:::initialize_fun(A = (sim_data$W>0)*1, 
                                             K = 3,
                                             D = 2,
                                             family = family,
                                             noise_weights = T,
                                             prob_matrix_W= prob_matrix_W,
                                             priors = NULL, 
                                             list_name = starting_params, 
                                             model = model, 
                                             n_interior_knots = control$n_interior_knots,
                                             n_control = control$n_control)
            
            A =(sim_data$W>0)*1
            
            # prob_mat -----------------------------------------------------------------
            
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
            
            prob_mat_R <- function(mus, omegas, p, U){
               K <- nrow(mus)
               N <- nrow(U)
               test_p <- matrix(NA, nrow = N, ncol = K)
               
               for (i in 1:N){
                  out <- numeric(K)
                  for(k in 1:K){
                     out[k] <- log(p[k])+ dmvnorm_chol_quad_log(x = U[i, ], mean = mus[k, ], omega = omegas[,,k])
                  }
                  
                  max_val <- max(out)
                  p2 <- log(sum(exp(out - max_val)))
                  
                  test_p[i,] <- exp(out - max_val - p2)
                  
               }
               
               return(test_p)
            }
            
            JANE:::update_prob_matrix_DA(prob_matrix = current$prob_matrix, 
                                         mus = current$mus, omegas = current$omegas, 
                                         p = current$p, U = current$U, temp_beta = 1)
            
            
            
            all_test[[as.character(noise)]][[model]][[family]][["updated_prob_matrix_U"]] <- all.equal(prob_mat_R(mus = current$mus,
                                                                                                    omegas = current$omegas,
                                                                                                    p = current$p, 
                                                                                                    U = current$U),
                                                                                         current$prob_matrix)
            
            if(noise != 0){
            # prob_mat_W -----------------------------------------------------------------
            
            prob_mat_W_R <- function(edge_indicies, model, beta, U, X, q, X2, family, beta2, tau, tau_noise, guess_noise_weights){
               
               test_p <- edge_indicies
               test_p[, 4] <- 0
               test_p[, 5] <- 0
               
               for (m in 1:nrow(edge_indicies)){
                  
                  i <- edge_indicies[m, 1]
                  j <- edge_indicies[m, 2]
                  diff <- U[i,] - U[j,]
                  cross_prod <- crossprod(diff)[1,1]
                  
                  if (model == "NDH"){
                     exp_eta <- exp(beta-cross_prod)
                     p <- (exp_eta/(1 + exp_eta))
                     
                     if (family !="bernoulli"){
                        eta_w <- beta2
                     }
                  } else if (model == "RS"){
                     x_ij <- c(1, X[i,] + X[j,])
                     beta_x <- (t(x_ij) %*% beta)[1,1]
                     exp_eta <- exp(beta_x-cross_prod) 
                     p <- (exp_eta/(1 + exp_eta))
                     
                     if (family !="bernoulli"){
                        x_ij <- c(1, X2[i,] + X2[j,])
                        beta_x <- (t(x_ij) %*% beta2)[1,1]
                        eta_w <- beta_x
                     }
                     
                  }else {
                     x_ij <- c(1, X[i, 1:(ncol(X) * 0.5)], X[j, (ncol(X)*0.5 + 1):ncol(X) ])
                     beta_x <- (t(x_ij) %*% beta)[1,1]
                     exp_eta <- exp(beta_x-cross_prod) 
                     p <- (exp_eta/(1 + exp_eta))
                     
                     if (family !="bernoulli"){
                        x_ij <- c(1, X2[i, 1:(ncol(X2) * 0.5)], X2[j, (ncol(X2)*0.5 + 1):ncol(X2) ])
                        beta_x <- (t(x_ij) %*% beta2)[1,1]
                        eta_w <- beta_x
                     }
                     
                  }
                  
                  if (family == "bernoulli"){
                     log_density_noise <- 0
                     log_density_w <- 0
                     z_hat1 <- exp( -1.0*log( 1.0 + exp( log(q*(1.0-p)) + log_density_noise - log(p) -  log_density_w) )  )
                  } else if (family == "poisson"){
                     w <- edge_indicies[m, 3]
                     log_density_w <- extraDistr::dtpois(x = w, lambda = exp(eta_w), a = 0.0, log = T)
                     log_density_noise<- extraDistr::dtpois(x = w, lambda = guess_noise_weights, a = 0.0, log = T)
                     z_hat1 <- exp( -1.0*log( 1.0 + exp( log(q*(1.0-p)) + log_density_noise - log(p) -  log_density_w) )  )
                  } else{
                     w <- edge_indicies[m, 3]
                     log_density_w <- dlnorm(x = w, meanlog = eta_w, sdlog = 1/sqrt(tau), log = T)
                     log_density_noise<-  dlnorm(x = w, meanlog = guess_noise_weights, sdlog = 1/sqrt(tau_noise), log = T)
                     z_hat1 <- exp( -1.0*log( 1.0 + exp( log(q*(1.0-p)) + log_density_noise - log(p) -  log_density_w) )  )
                  } 
                  
                  if(is.nan(z_hat1)){
                     z_hat1 <- 0.5
                  }
                  
                  test_p[m, 4] <- z_hat1
                  test_p[m, 5] <- 1-z_hat1
                  
                  
                  
               }
               
               return(test_p)
            }
            
            temp <- prob_mat_W_R(edge_indicies = current$prob_matrix_W,
                                 model = current$model,
                                 beta = current$beta,
                                 U = current$U,
                                 X = current$X,
                                 q = current$q_prob,
                                 X2 = current$X2,
                                 family = current$family,
                                 beta2 = current$beta2,
                                 tau_noise = current$precision_noise_weights,
                                 tau = current$precision_weights,
                                 guess_noise_weights = current$guess_noise_weights)
            
            
            JANE:::update_prob_matrix_W_DA(prob_matrix_W = current$prob_matrix_W,
                                           model = current$model,
                                           family = current$family,
                                           beta = current$beta,
                                           beta2 = if(family == "bernoulli"){matrix(0.0, 1, 1)}else{current$beta2},
                                           precision_weights = ifelse(family != "lognormal", 0.0, current$precision_weights),
                                           precision_noise_weights = ifelse(family != "lognormal", 0.0, current$precision_noise_weights),
                                           guess_noise_weights = current$guess_noise_weights,
                                           U = current$U,
                                           X = current$X,
                                           X2 = if(family == "bernoulli"){matrix(0.0, 1, 1)}else{current$X2},
                                           q = current$q_prob,
                                           temp_beta =1)
            
            
            all_test[[as.character(noise)]][[model]][[family]][["updated_prob_matrix_W"]] <- all.equal(temp,
                                                                                         current$prob_matrix_W)
            
            A2 <- A
            A2[current$prob_matrix_W[, 1:2]] <- current$prob_matrix_W[,4]
            
            if(model !="RSR"){
               A2[current$prob_matrix_W[, 2:1]] <- current$prob_matrix_W[,4]
            }
            
            # q_prob -----------------------------------------------------------------
            
            q_prob_R <- function(A, prob_matrix_W, model, h,l){
               
               
               if (model == "RSR"){
                  
                  p1 = 0
                  p2 = 0
                  
                  for (i in 1:nrow(prob_matrix_W)){
                     for (j in 1:nrow(prob_matrix_W)){
                        if (i != j){
                           p1 <- p1+ (((A[i,j]>0)*1)*(1-prob_matrix_W[i,j]))
                           p2 <- p2 + (1-((A[i,j]>0)*1))
                        }
                     }
                  }
                  
                  q <- (p1+h-1)/((p1+h-1)+(p2+l-1))
                  
               } else {
                  
                  p1 = 0
                  p2 = 0
                  
                  for (i in 1:nrow(prob_matrix_W)){
                     for (j in 1:nrow(prob_matrix_W)){
                        if (j>i){
                           p1 <- p1+ (((A[i,j]>0)*1)*(1-prob_matrix_W[i,j]))
                           p2 <- p2 + (1-((A[i,j]>0)*1))
                        }
                     }
                  }
                  
                  q <- (p1+h-1)/((p1+h-1)+(p2+l-1))
               }
               
               return(q)
               
            }
            
            
            
            
            JANE:::update_q_prob(q_prob = current$q_prob,
                                 prob_matrix_W = current$prob_matrix_W, 
                                 N = nrow(A),
                                 model = current$model, 
                                 h = current$priors$h,
                                 l= current$priors$h)
            
            
            all_test[[as.character(noise)]][[model]][[family]][["updated_q"]] <- all.equal(q_prob_R(A,
                                                                                      prob_matrix_W = A2, 
                                                                                      model = current$model, 
                                                                                      h = current$priors$h,
                                                                                      l= current$priors$h),
                                                                             current$q_prob)
            } else {
               A2 <- A
               A2[current$prob_matrix_W[, 1:2]] <- current$prob_matrix_W[,4]
               
               if(model !="RSR"){
                  A2[current$prob_matrix_W[, 2:1]] <- current$prob_matrix_W[,4]
               }
            }
            
            # mu_omega ----------------------------------------------------------
            
            mus_omegas_R <- function(prob_matrix, U, b, a, c, G, mus, omegas){
               
               D <- ncol(U)
               K <- ncol(prob_matrix)
               N <- nrow(U)
               mus_out <- matrix(NA, nrow = K, ncol = D)
               omegas_out <- array(NA, dim = c(D,D,K))
               for (k in 1:K){
                  sum_hold <- numeric(D)
                  sum_hold2 <- 0
                  sum_hold3 <- matrix(0, D, D)
                  for (i in 1:N){
                     sum_hold <- sum_hold + prob_matrix[i,k]*U[i,]
                     sum_hold2 <- sum_hold2 + prob_matrix[i,k]
                     sum_hold3 <- sum_hold3 + prob_matrix[i,k]* (U[i,] %*% t(U[i,]))
                  }
                  mus_out[k,] <- (sum_hold + b*a[1,])/(sum_hold2 + b)
                  
                  temp1 <- sum_hold2 + c - D
                  temp2 <- G + sum_hold3 + (b* t(a) %*% a) - ((sum_hold2 + b) * mus_out[k,] %*% t(mus_out[k,]))
                  omegas_out[,,k] <- temp1 * chol2inv(chol(temp2))
               }
               return(list(mus = mus_out,
                           omegas = omegas_out))
            }  
            
            temp_updated_mus_omegas <- mus_omegas_R(prob_matrix = current$prob_matrix, U = current$U,
                                                    b = current$priors$b, a = current$priors$a, c = current$priors$c, G = current$priors$G, 
                                                    mus = current$mus, omegas = ocurrent$megas)
            
            JANE:::update_mus_omegas(prob_matrix = current$prob_matrix,
                                     U = current$U, b = current$priors$b, a = current$priors$a,
                                     c = current$priors$c, G = current$priors$G,
                                     mus = current$mus, omegas = current$omegas)
            
            all_test[[as.character(noise)]][[model]][[family]][["updated_mus"]] <- all.equal(temp_updated_mus_omegas$mus,
                                                                               current$mus, check.attributes = F)
            
            all_test[[as.character(noise)]][[model]][[family]][["updated_omegas"]] <- all.equal(temp_updated_mus_omegas$omegas,
                                                                                  current$omegas)
            
            # update p ---------------------------------------------------------------------
            
            update_p_R <- function(prob_matrix, nu){
               
               N <- nrow(prob_matrix)
               K <- ncol(prob_matrix)
               p_new <- numeric(K)
               for (k in 1:K){
                  
                  p1 <- 0
                  for (i in 1:N){
                     p1 <- p1 + prob_matrix[i,k]
                  }
                  p_new[k] <- p1 + nu[k] - 1
               }
               
               p_new <- p_new/sum(p_new)
               
               return(p_new)
            }
            
            temp_update_p <- update_p_R(prob_matrix = current$prob_matrix, nu = current$priors$nu)
            JANE:::update_p(prob_matrix = current$prob_matrix, p = current$p, nu = current$priors$nu)
            all_test[[as.character(noise)]][[model]][[family]][["updated_p"]] <- all.equal(temp_update_p,
                                                                                           current$p)
            # Update U -----------------------------------------------------------------
            
            if(model == "NDH"){
               
               update_U_R <- function(U, A, mus, omegas, prob_matrix, beta0){
                  
                  N <- nrow(U)
                  D <- ncol(U)
                  K <- nrow(mus)
                  
                  U_prime <- U
                  
                  for (i in 1:N){
                     
                     p1_1 <- 0
                     for (k in 1:K){
                        p1_1 <- p1_1 + prob_matrix[i,k] * omegas[,,k]
                     }
                     
                     p1_2 <- 0
                     for (j in 1:N){
                        if (j != i){
                           p1_2 <- p1_2 + 2*A[i,j]
                        }
                     }
                     
                     p1 <- chol2inv(chol(p1_1 + p1_2 *diag(D) ))
                     
                     p2_1 <- numeric(D)
                     for (k in 1:K){
                        p2_1 <- p2_1 + prob_matrix[i,k] * (omegas[,,k] %*% mus[k,])[,1]
                     }
                     
                     p2_2 <- numeric(D)
                     for (j in 1:N){
                        if (j != i){
                           p2_2 <- p2_2 + 2*(A[i,j]) * U_prime[j,]
                        }
                     }
                     
                     p2_3 <-  numeric(D)
                     for (j in 1:N){
                        if(j != i){
                           diff <- U_prime[i,] - U_prime[j,]
                           cross_prod <- crossprod(diff)[1,1]
                           exp_eta <- exp(beta0-cross_prod)
                           p <- (exp_eta/(1 + exp_eta))
                           p2_3 <- p2_3 + 2*p*diff
                        }
                     }
                     
                     U_prime[i, ] <- (p1 %*% (p2_1 + p2_2 + p2_3))[,1]
                     
                  }
                  
                  return(U_prime)
               }
               
               temp_update_U <- update_U_R(U = current$U, A = A, 
                                           mus = current$mus, omegas = current$omegas,
                                           prob_matrix = current$prob_matrix, 
                                           beta0 = current$beta)
               
               JANE:::update_U(U = current$U, A = A,  
                               mus = current$mus, omegas = current$omegas, prob_matrix = current$prob_matrix, 
                               beta = current$beta,
                               X = current$X,
                               n_control = NULL,
                               model = NULL)
               
               all_test[[as.character(noise)]][[model]][[family]][["updated_U"]] <- all.equal(temp_update_U,
                                                                                              current$U)
               
            } else if( model == "RS"){
               
               update_U_RS_R <- function(U, A, mus, omegas, prob_matrix, beta, X){
                  
                  N <- nrow(U)
                  D <- ncol(U)
                  K <- nrow(mus)
                  
                  U_prime <- U
                  
                  for (i in 1:N){
                     
                     p1_1 <- 0
                     for (k in 1:K){
                        p1_1 <- p1_1 + prob_matrix[i,k] * omegas[,,k]
                     }
                     
                     p1_2 <- 0
                     for (j in 1:N){
                        if (j != i){
                           p1_2 <- p1_2 + 2*A[i,j]
                        }
                     }
                     
                     p1 <- chol2inv(chol(p1_1 + p1_2 *diag(D) ))
                     
                     p2_1 <- numeric(D)
                     for (k in 1:K){
                        p2_1 <- p2_1 + prob_matrix[i,k] * (omegas[,,k] %*% mus[k,])[,1]
                     }
                     
                     p2_2 <- numeric(D)
                     for (j in 1:N){
                        if (j != i){
                           p2_2 <- p2_2 + 2*(A[i,j]) * U_prime[j,]
                        }
                     }
                     
                     p2_3 <-  numeric(D)
                     for (j in 1:N){
                        if(j != i){
                           diff <- U_prime[i,] - U_prime[j,]
                           cross_prod <- crossprod(diff)[1,1]
                           x_ij <- c(1, X[i,] + X[j,])
                           beta_x <- (t(x_ij) %*% beta)[1,1]
                           exp_eta <- exp(beta_x-cross_prod )
                           p <- (exp_eta/(1 + exp_eta))
                           p2_3 <- p2_3 + 2*p*diff
                        }
                     }
                     
                     U_prime[i, ] <- (p1 %*% (p2_1 + p2_2 + p2_3))[,1]
                     
                  }
                  
                  return(U_prime)
               }
               
               temp_update_U <- update_U_RS_R(U = current$U, A = A, 
                                              mus = current$mus, omegas = current$omegas, prob_matrix = current$prob_matrix, 
                                              beta = current$beta,
                                              X = current$X)
               
               JANE:::update_U_RE(U = current$U, A = A, 
                                  mus = current$mus, omegas = current$omegas, prob_matrix = current$prob_matrix, 
                                  beta = current$beta, 
                                  X = current$X,
                                  n_control = NULL,
                                  model = current$model)
               
               all_test[[as.character(noise)]][[model]][[family]][["updated_U"]] <- all.equal(temp_update_U,
                                                                                              current$U)
               
               
            } else {
               
               update_U_RSR_R <- function(U, A, mus, omegas, prob_matrix, beta, X){
                  
                  N <- nrow(U)
                  D <- ncol(U)
                  K <- nrow(mus)
                  
                  U_prime <- U
                  
                  for (i in 1:N){
                     
                     p1_1 <- 0
                     for (k in 1:K){
                        p1_1 <- p1_1 + prob_matrix[i,k] * omegas[,,k]
                     }
                     
                     p1_2 <- 0
                     for (j in 1:N){
                        if (j != i){
                           p1_2 <- p1_2 + 2*(A[i,j] + A[j,i])
                        }
                     }
                     
                     p1 <- chol2inv(chol(p1_1 + p1_2 *diag(D) ))
                     
                     p2_1 <- numeric(D)
                     for (k in 1:K){
                        p2_1 <- p2_1 + prob_matrix[i,k] * (omegas[,,k] %*% mus[k,])[,1]
                     }
                     
                     p2_2 <- numeric(D)
                     for (j in 1:N){
                        if (j != i){
                           p2_2 <- p2_2 + 2*(A[i,j] + A[j,i]) * U_prime[j,]
                        }
                     }
                     
                     p2_3 <-  numeric(D)
                     for (j in 1:N){
                        if(j != i){
                           
                           diff <- U_prime[i,] - U_prime[j,]
                           cross_prod <- crossprod(diff)[1,1]
                           
                           x_ij <- c(1, X[i, 1:(ncol(X)*0.5)], X[j, (ncol(X)*0.5 +1):ncol(X)])
                           beta_x <- (t(x_ij) %*% beta)[1,1]
                           exp_eta <- exp(beta_x-cross_prod )
                           p_ij <- (exp_eta/(1 + exp_eta))
                           
                           
                           x_ji <- c(1, X[j, 1:(ncol(X)*0.5)], X[i, (ncol(X)*0.5 +1):ncol(X)])
                           beta_x <- (t(x_ji) %*% beta)[1,1]
                           exp_eta <- exp(beta_x-cross_prod )
                           p_ji <- (exp_eta/(1 + exp_eta))
                           
                           p2_3 <- p2_3 + (2*(p_ij + p_ji)*diff)
                           
                        }
                     }
                     
                     U_prime[i, ] <- (p1 %*% (p2_1 + p2_2 + p2_3))[,1]
                     
                  }
                  
                  return(U_prime)
               }
               
               temp_update_U <- update_U_RSR_R(U = current$U, A = A, 
                                               mus = current$mus, omegas = current$omegas, prob_matrix = current$prob_matrix, 
                                               beta = current$beta,
                                               X = current$X)
               
               JANE:::update_U_RE(U = current$U, A = A, 
                                  mus = current$mus, omegas = current$omegas, prob_matrix = current$prob_matrix, 
                                  beta = current$beta, 
                                  X = current$X,
                                  n_control = NULL,
                                  model = current$model)
               all_test[[as.character(noise)]][[model]][[family]][["updated_U"]] <- all.equal(temp_update_U,
                                                                                              current$U)
               
            }
            
            # beta -----------------------------------------------------------------
            
            if (model == "NDH"){
               update_beta_R <- function(A, U, beta0, f, e){
                  
                  N <- nrow(U)
                  p1 <- f*e
                  p2 <- 0
                  for (i in 1:N){
                     for (j in 1:N){
                        if(j!=i){
                           p2 <- p2 + 0.5*A[i,j]
                        }
                     }
                  }
                  
                  p3 <- 0
                  for (i in 1:N){
                     for (j in 1:N){
                        if(j!=i){
                           diff <- U[i,] - U[j,]
                           cross_prod <- crossprod(diff)[1,1]
                           exp_eta <- exp(beta0-cross_prod)
                           p <- (exp_eta/(1 + exp_eta))
                           
                           p3 <- p3 + ( 0.5* ( (p*(1-p))*beta0 - p ) )
                        }
                     }
                  }
                  
                  p4 <- 0
                  for (i in 1:N){
                     for (j in 1:N){
                        if(j!=i){
                           diff <- U[i,] - U[j,]
                           cross_prod <- crossprod(diff)[1,1]
                           exp_eta <- exp(beta0-cross_prod)
                           p <- (exp_eta/(1 + exp_eta))
                           
                           p4 <- p4 + ( 0.5* ( (p*(1-p))) )
                        }
                     }
                  }
                  
                  beta_0_prime <- (p1 + p2 + p3)/(f + p4)
                  
                  return(beta_0_prime)
               }
               
               temp_update_beta <- update_beta_R(A = A, U = current$U, beta0 = current$beta, 
                                                 f = current$priors$f, e = current$priors$e)
               
               JANE:::update_beta(A = A, n_control = NULL, 
                                  U = current$U, beta = current$beta, 
                                  f = current$priors$f, e = current$priors$e,
                                  X = NULL,
                                  model = NULL)
               all_test[[as.character(noise)]][[model]][[family]][["updated_beta_LR"]] <- all.equal(temp_update_beta,
                                                                                                    current$beta)
               
               
            } else if(model == "RS"){
               update_beta_RS_R <- function(A, U, beta, f, e, X){
                  
                  N <- nrow(U)
                  
                  p1 <- matrix(0, length(beta), length(beta))
                  p2 <- matrix(0, nrow = length(beta), ncol = 1 )
                  p3 <-  matrix(0, nrow = length(beta), ncol = 1 )
                  
                  for (i in 1:N){
                     for(j in 1:N){
                        if(j !=i){
                           diff <- U[i,] - U[j,]
                           cross_prod <- crossprod(diff)[1,1]
                           x_ij <- c(1, X[i,] + X[j,])
                           beta_x <- (t(x_ij) %*% beta)[1,1]
                           exp_eta <- exp(beta_x-cross_prod) 
                           p <- (exp_eta/(1 + exp_eta))
                           
                           p1 <- p1 + ( (0.5*(p*(1-p))) * (x_ij %*% t(x_ij)) ) 
                           p2 <- p2 + ( 0.5*A[i,j]*x_ij)
                           p3 <- p3 + ( 0.5* ( ( (p*(1-p))* (x_ij %*% t(x_ij) %*% beta) ) - (p*x_ij) ) )
                        }
                     }
                  }
                  
                  beta_0_prime = chol2inv(chol(f + p1)) %*% (f%*%e + p2 + p3)
                  
                  return(beta_0_prime)
               }
               
               temp_update_beta <- update_beta_RS_R(A = A2, U = current$U, beta = current$beta, 
                                                    X = current$X,
                                                    f = current$priors$f, e = current$priors$e)[,1]
            
               JANE:::update_beta_RE(A = A2, U = current$U, beta = current$beta, 
                                     X = current$X,
                                     f = current$priors$f, e = current$priors$e,
                                     n_control = NULL,
                                     model = current$model)
               all_test[[as.character(noise)]][[model]][[family]][["updated_beta_LR"]] <- all.equal(temp_update_beta,
                                                                                                    current$beta)
               
            }else {
               
               update_beta_RSR_R <- function(A, U, beta, f, e, X){
                  
                  N <- nrow(U)
                  
                  p1 <- matrix(0, length(beta), length(beta))
                  p2 <- matrix(0, nrow = length(beta), ncol = 1 )
                  p3 <-  matrix(0, nrow = length(beta), ncol = 1 )
                  
                  for (i in 1:N){
                     for(j in 1:N){
                        if(j !=i){
                           diff <- U[i,] - U[j,]
                           cross_prod <- crossprod(diff)[1,1]
                           x_ij <- c(1, X[i, 1:(ncol(X) * 0.5)], X[j, (ncol(X)*0.5 + 1):ncol(X) ])
                           beta_x <- (t(x_ij) %*% beta)[1,1]
                           exp_eta <- exp(beta_x-cross_prod) 
                           p <- (exp_eta/(1 + exp_eta))
                           
                           p1 <- p1 + ( ((p*(1-p))) * (x_ij %*% t(x_ij)) ) 
                           p2 <- p2 + ( A[i,j]*x_ij)
                           p3 <- p3 + (  ( ( (p*(1-p))* (x_ij %*% t(x_ij) %*% beta) ) - (p*x_ij) ) )
                        }
                     }
                  }
                  
                  beta_0_prime = chol2inv(chol(f + p1)) %*% (f%*%e + p2 + p3)
                  
                  return(beta_0_prime)
               }
               
               temp_update_beta <- update_beta_RSR_R(A = A, U = current$U, beta = current$beta, 
                                                     X = current$X,
                                                     f = current$priors$f, e = current$priors$e)[,1]
               #modifies in place
               JANE:::update_beta_RE(A = A, U = current$U, beta = current$beta, 
                                     X = current$X,
                                     f = current$priors$f, e = current$priors$e,
                                     n_control = NULL,
                                     model = current$model)
               all_test[[as.character(noise)]][[model]][[family]][["updated_beta_LR"]] <- all.equal(temp_update_beta,
                                                                                                    current$beta)
            }
            
            if(noise != 0){
               # beta2 -----------------------------------------------------------------
               
               update_beta2_R <- function(A2, W, model, family, beta2, X2, f, e){
                  
                  if(family == "lognormal"){
                     
                     if (model == "NDH"){
                        
                        p1 <- 0
                        p2 <- 0
                        
                        for (i in 1:nrow(W)){
                           for (j in 1:nrow(W)){
                              
                              if (j != i){
                                 
                                 ind <- (W[i,j] >0.0)*1
                                 if(ind){
                                    p1 <- p1 + ind*A2[i,j]*log(W[i,j])
                                    p2 <- p2 + ind*A2[i,j]
                                 }
                              }
                              
                           }
                        }
                        
                        beta2_out <- (0.5*p1 + e*f)/(0.5*p2 + f)
                        
                     } else if (model == "RS"){
                        
                        p1 <- rep(0, ncol(X2)+1)
                        p2 <- matrix(0,ncol(X2)+1,ncol(X2)+1)
                        
                        for (i in 1:nrow(W)){
                           for (j in 1:nrow(W)){
                              
                              if (j != i){
                                 
                                 x2_ij <- c(1, X2[i,] + X2[j,])
                                 
                                 ind <- (W[i,j] >0.0)*1
                                 if(ind){
                                    p1 <- p1 + ind*A2[i,j]*log(W[i,j])*x2_ij
                                    p2 <- p2 + ind*A2[i,j]*(tcrossprod(x2_ij))
                                 }
                              }
                              
                           }
                        }
                        
                        beta2_out <- chol2inv(chol(0.5*p2 + f)) %*% (0.5*p1 + f%*%e )
                        
                     } else {
                        
                        p1 <- rep(0, ncol(X2)+1)
                        p2 <- matrix(0,ncol(X2)+1,ncol(X2)+1)
                        
                        for (i in 1:nrow(W)){
                           for (j in 1:nrow(W)){
                              
                              if (j != i){
                                 
                                 x2_ij <- c(1, X2[i, 1:(ncol(X2) * 0.5)], X2[j, (ncol(X2)*0.5 + 1):ncol(X2) ])
                                 ind <- (W[i,j] >0.0)*1
                                 if(ind){
                                    p1 <- p1 + ind*A2[i,j]*log(W[i,j])*x2_ij
                                    p2 <- p2 + ind*A2[i,j]*(tcrossprod(x2_ij))
                                 }
                              }
                              
                           }
                        }
                        
                        beta2_out <- chol2inv(chol(p2 + f)) %*% (p1 + f%*%e )
                        
                        
                     }
                     
                  } else {
                     
                     
                     if (model == "NDH"){
                        
                        p1 <- 0
                        p2 <- 0
                        
                        for (i in 1:nrow(W)){
                           for (j in 1:nrow(W)){
                              
                              if (j != i){
                                 
                                 ind <- (W[i,j] >0.0)*1
                                 if(ind){
                                    lambda <- exp(beta2)
                                    g1 <- lambda/(1-exp(-1*lambda))
                                    g2 <- (lambda * (1-exp(-1*lambda)*(1+lambda)))/(1-exp(-1*lambda))^2 
                                    
                                    
                                    p1 <- p1 + ind*A2[i,j]*(W[i,j] + g2*beta2 - g1)
                                    p2 <- p2 + ind*A2[i,j]*g2
                                 }
                              }
                              
                           }
                        }
                        
                        beta2_out <- (0.5*p1 + e*f)/(0.5*p2 + f)
                        
                     } else if (model == "RS"){
                        
                        p1 <- rep(0, ncol(X2)+1)
                        p2 <- matrix(0,ncol(X2)+1,ncol(X2)+1)
                        
                        for (i in 1:nrow(W)){
                           for (j in 1:nrow(W)){
                              
                              if (j != i){
                                 
                                 x2_ij <- c(1, X2[i,] + X2[j,])
                                 
                                 ind <- (W[i,j] >0.0)*1
                                 if(ind){
                                    lambda <- exp(t(x2_ij)%*%beta2)[1,1]
                                    g1 <- (lambda/(1-exp(-1*lambda)))*x2_ij
                                    g2 <- ((lambda * (1-exp(-1*lambda)*(1+lambda)))/(1-exp(-1*lambda))^2 ) * tcrossprod(x2_ij)
                                    
                                    p1 <- p1 + ind*A2[i,j]*(W[i,j]*x2_ij + g2 %*% beta2 - g1)
                                    p2 <- p2 + ind*A2[i,j]*g2
                                 }
                              }
                              
                           }
                        }
                        
                        beta2_out <- chol2inv(chol(0.5*p2 + f)) %*% (0.5*p1 + f%*%e )
                        
                     } else {
                        
                        p1 <- rep(0, ncol(X2)+1)
                        p2 <- matrix(0,ncol(X2)+1,ncol(X2)+1)
                        
                        for (i in 1:nrow(W)){
                           for (j in 1:nrow(W)){
                              
                              if (j != i){
                                 
                                 x2_ij <- c(1, X2[i, 1:(ncol(X2) * 0.5)], X2[j, (ncol(X2)*0.5 + 1):ncol(X2) ])
                                 ind <- (W[i,j] >0.0)*1
                                 if(ind){
                                    lambda <- exp(t(x2_ij)%*%beta2)[1,1]
                                    g1 <- (lambda/(1-exp(-1*lambda)))*x2_ij
                                    g2 <- ((lambda * (1-exp(-1*lambda)*(1+lambda)))/(1-exp(-1*lambda))^2 ) * tcrossprod(x2_ij)
                                    
                                    p1 <- p1 + ind*A2[i,j]*(W[i,j]*x2_ij + g2 %*% beta2 - g1)
                                    p2 <- p2 + ind*A2[i,j]*g2
                                 }
                              }
                              
                           }
                        }
                        
                        beta2_out <- chol2inv(chol(p2 + f)) %*% (p1 + f%*%e )
                        
                        
                     }
                     
                     
                     
                  }
                  
                  return(beta2_out)
                  
               }
               
               if(family != "bernoulli"){
                  
                  temp <- update_beta2_R(A2 = A2,
                                         W = sim_data$W,
                                         model = model,
                                         family = family,
                                         beta2 = current$beta2, 
                                         X2 = current$X2, 
                                         f = current$priors$f_2, 
                                         e = current$priors$e_2)
                  
                  if (model != "NDH"){
                     temp <- temp[,1]
                  }
                  JANE:::update_beta2(beta2 = current$beta2,
                                      prob_matrix_W = current$prob_matrix_W,
                                      X2 = current$X2, 
                                      model = current$model, 
                                      family = current$family,
                                      f_2 = if(model == "NDH"){matrix(current$priors$f_2, 1, 1)}else{current$priors$f_2}, 
                                      e_2 = if(model == "NDH"){matrix(current$priors$e_2, 1, 1)}else{current$priors$e_2})
                  
                  all_test[[as.character(noise)]][[model]][[family]][["updated_beta2"]] <- all.equal(temp,
                                                                                       current$beta2)
                  
               }
               
               # tau -----------------------------------------------------------------
               
               update_tau_R <- function(A2, 
                                        W, 
                                        model,
                                        beta2,
                                        X2,
                                        m_1,
                                        o_1,
                                        f, 
                                        e){
                  
                  if(model == "NDH"){
                     p1 <- 0
                     p2 <- 0
                     
                     for (i in 1:nrow(W)){
                        for (j in 1:nrow(W)){
                           
                           if(j>i){
                              
                              ind <- (W[i,j] >0.0)*1
                              if(ind){
                                 
                                 p1 <- p1 + ind*A2[i,j]
                                 p2 <- p2 + ind*A2[i,j]*log(W[i,j])^2 - 2*beta2*ind*A2[i,j]*log(W[i,j])+ (beta2^2)*ind*A2[i,j]
                                 
                              }
                              
                           }
                           
                        }
                     }
                     
                     tau_out <- (p1 +m_1-1)/(p2 + f*((beta2 - e)^2) + o_1)
                     
                     
                  } else if (model == "RS"){
                     
                     p1 <- 0
                     p2 <- 0
                     
                     for (i in 1:nrow(W)){
                        for (j in 1:nrow(W)){
                           
                           if(j>i){
                              x2_ij <- c(1, X2[i,] + X2[j,])
                              ind <- (W[i,j] >0.0)*1
                              if(ind){
                                 
                                 p1 <- p1 + ind*A2[i,j]
                                 p2 <- p2 + ind*A2[i,j]*log(W[i,j])^2 + ind*A2[i,j]*t(beta2)%*%x2_ij%*%t(x2_ij)%*% beta2 - 2*ind*A2[i,j]*log(W[i,j])*t(beta2)%*%x2_ij
                                 
                              }
                              
                           }
                           
                        }
                     }
                     
                     tau_out <- (p1+m_1+length(beta2)-2)/(p2 + t(beta2-e)%*% f %*% (beta2-e)+ o_1)
                     
                  } else {
                     
                     p1 <- 0
                     p2 <- 0
                     
                     for (i in 1:nrow(W)){
                        for (j in 1:nrow(W)){
                           
                           if(j!=i){
                              x2_ij <- c(1, X2[i, 1:(ncol(X2) * 0.5)], X2[j, (ncol(X2)*0.5 + 1):ncol(X2) ])
                              ind <- (W[i,j] >0.0)*1
                              if(ind){
                                 
                                 p1 <- p1 + ind*A2[i,j]
                                 p2 <- p2 + ind*A2[i,j]*log(W[i,j])^2 +  ind*A2[i,j]*t(beta2)%*%x2_ij%*%t(x2_ij)%*% beta2 - 2*ind*A2[i,j]*log(W[i,j])*t(beta2)%*%x2_ij
                                 
                              }
                              
                           }
                           
                        }
                     }
                     
                     tau_out <- (p1+m_1+length(beta2)-2)/(p2 + (t(beta2-e)%*% f %*% (beta2-e)) + o_1)
                     
                  }
                  
                  return(tau_out)
                  
               }
               
               
               if(family == "lognormal"){
                  
                  temp <- update_tau_R(A2 =A2,
                                       W = sim_data$W,
                                       model = current$model,
                                       beta2 = current$beta2,
                                       X2 = current$X2,
                                       f = current$priors$f_2, 
                                       e = current$priors$e_2,
                                       m_1 = current$priors$m_1,
                                       o_1 = current$priors$o_1)
                  
                  
                  JANE:::update_precision_weights(precision_weights = current$precision_weights,
                                                  beta2 = current$beta2,
                                                  prob_matrix_W = current$prob_matrix_W,
                                                  X2 = current$X2, 
                                                  model = current$model, 
                                                  f_2 = if(current$model == "NDH"){matrix(current$priors$f_2, 1, 1)}else{current$priors$f_2}, 
                                                  e_2 = if(current$model == "NDH"){matrix(current$priors$e_2, 1, 1)}else{current$priors$e_2},
                                                  m_1 = current$priors$m_1,
                                                  o_1 = current$priors$o_1)
                  
                  all_test[[as.character(noise)]][[model]][[family]][["update_precision_weights"]] <- all.equal(c(temp),
                                                                                                  current$precision_weights)
                  
               }
               
               
               # tau_noise -----------------------------------------------------------------
               
               update_tau_noise_R <- function(A2, 
                                              W, 
                                              model,
                                              guess_noise_weights,
                                              m_2,
                                              o_2){
                  
                  if(model != "RSR"){
                     p1 <- 0
                     p2 <- 0
                     for (i in 1:nrow(W)){
                        for (j in 1:nrow(W)){
                           
                           if(j>i){
                              
                              ind <- (W[i,j] >0.0)*1
                              if(ind){
                                 
                                 p1 <- p1 + ind*(1-A2[i,j])
                                 p2 <- p2 + ind*(1-A2[i,j])*(log(W[i,j]) - guess_noise_weights)^2
                                 
                              }
                              
                           }
                           
                        }
                     }
                     
                     tau_out <- (p1+m_2-2)/(p2+o_2)
                     
                  } else {
                     p1 <- 0
                     p2 <- 0
                     for (i in 1:nrow(W)){
                        for (j in 1:nrow(W)){
                           
                           if(j!=i){
                              
                              ind <- (W[i,j] >0.0)*1
                              if(ind){
                                 
                                 p1 <- p1 + ind*(1-A2[i,j])
                                 p2 <- p2 + ind*(1-A2[i,j])*(log(W[i,j]) - guess_noise_weights)^2
                                 
                              }
                              
                           }
                           
                        }
                     }
                     
                     tau_out <- (p1+m_2-2)/(p2+o_2)
                     
                     
                  }
                  
                  return(tau_out)
               }
               
               
               if(family == "lognormal"){
                  
                  temp <- update_tau_noise_R(A2 = A2, 
                                             W = sim_data$W, 
                                             model = current$model,
                                             guess_noise_weights = current$guess_noise_weights,
                                             m_2 = current$priors$m_2,
                                             o_2 = current$priors$o_2)
                  
                  
                  JANE:::update_precision_noise_weights(precision_noise_weights = current$precision_noise_weights,
                                                        prob_matrix_W = current$prob_matrix_W,
                                                        guess_noise_weights = current$guess_noise_weights,
                                                        m_2 = current$priors$m_2,
                                                        o_2 = current$priors$o_2)
                  
                  all_test[[as.character(noise)]][[model]][[family]][["update_precision_noise_weights"]] <- all.equal(c(temp),
                                                                                                        current$precision_noise_weights)
                  
               }
               
            }
            
            if(any(!unlist(all_test))){
               break
            }
            
         }
      }
      
   }
   
   expect_true(all(unlist(all_test)))
   
})

