
## Testing sim_A

test_that("sim_A works", {
  
  expect_no_error({
  
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
    
     # Simulate an undirected, unweighted network, with no noise and no degree heterogeneity
     JANE::sim_A(N = 100L, 
                 model = "NDH",
                 mus = mus, 
                 omegas = omegas, 
                 p = p, 
                 params_LR = list(beta0 = beta0),
                 remove_isolates = TRUE)
    
     # Simulate a directed, weighted network, with no noise and degree heterogeneity
     JANE::sim_A(N = 100L, 
                 model = "RSR",
                 family = "lognormal",
                 mus = mus, 
                 omegas = omegas, 
                 p = p, 
                 params_LR = list(beta0 = beta0),
                 params_weights = list(beta0 = 2,
                                       precision_weights = 1),
                 remove_isolates = TRUE)
    
     # Simulate an undirected, weighted network, with noise and degree heterogeneity
     JANE::sim_A(N = 100L, 
                 model = "RS",
                 family = "poisson",
                 mus = mus, 
                 omegas = omegas, 
                 p = p, 
                 params_LR = list(beta0 = beta0),
                 params_weights = list(beta0 = 2),
                 noise_weights_prob = 0.1,
                 mean_noise_weights = 1,
                 remove_isolates = TRUE)
                 
  })
  
})


test_that("sim_A comprehensive works", {
  
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
  p
  
  main <- c()
  for (i in c("NDH", "RS", "RSR")){
    for (j in c("bernoulli", "poisson", "lognormal")){
      for (k in c(0, 0.1)){
        model <- i
        family <- j
        noise_weights_prob <- k
        
        if(family == "lognormal"){
          params_weights <-  list(beta0 = 2,
                                  precision_weights = 1)
          precision_noise_weights <- 2
          mean_noise_weights <- 2
        } else if (family == "poisson"){
          params_weights <-  list(beta0 = 2)
          precision_noise_weights <- NULL
          mean_noise_weights <- 2
        } else {
          params_weights <- NULL
          precision_noise_weights <- NULL
          mean_noise_weights <- NULL
        }
        
        test <- c()
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
                                noise_weights_prob = noise_weights_prob,
                                remove_isolates = TRUE)
        
        if (noise_weights_prob > 0){
          if(model != "RSR"){
            test <- c(test, all.equal(dim(summary(sim_data$A)[summary(sim_data$A)$j>summary(sim_data$A)$i, ])[1], 
                                      unname(table(sim_data$Z_W[, "Z_W"])[1])))
            test <- c(test, all.equal(dim(summary(sim_data$W)[summary(sim_data$W)$j>summary(sim_data$W)$i, ])[1], unname(sum(table(sim_data$Z_W[, "Z_W"])))))
            
            q_prob <- table(sim_data$Z_W[, "Z_W"])[2]/( 0.5*nrow(sim_data$W)*(nrow(sim_data$W)-1) - table(sim_data$Z_W[, "Z_W"])[1])
            
            den_A <- sum(sim_data$A)/(nrow(sim_data$A)*(nrow(sim_data$A)-1))
            test <- c(test, all.equal((1-den_A)/((1-den_A) + (den_A/q_prob) ), prop.table(table(sim_data$Z_W[, "Z_W"]))[2]))
            
          } else {
            test <- c(test, all.equal(dim(summary(sim_data$A))[1], unname(table(sim_data$Z_W[, "Z_W"])[1])))
            test <- c(test, all.equal(dim(summary(sim_data$W))[1], unname(sum(table(sim_data$Z_W[, "Z_W"])))))
            
            q_prob <- table(sim_data$Z_W[, "Z_W"])[2]/( nrow(sim_data$W)*(nrow(sim_data$W)-1) - table(sim_data$Z_W[, "Z_W"])[1])
            
            den_A <- sum(sim_data$A)/(nrow(sim_data$A)*(nrow(sim_data$A)-1))
            test <- c(test, all.equal((1-den_A)/((1-den_A) + (den_A/q_prob) ), prop.table(table(sim_data$Z_W[, "Z_W"]))[2]))
            
          }
        }
        
        A <- sim_data$A
        W <- sim_data$W
        
        if (noise_weights_prob > 0){
          if (model !="RSR"){
            test <- c(test, isSymmetric(A) == T & isSymmetric(W) == T)
            
            Z_W <- sim_data$Z_W
            Z_W <- Z_W[Z_W[, "Z_W"] ==2,]
            
            W[Z_W[,c(1,2)]] <- 0
            W[Z_W[,c(2,1)]] <- 0
            
            Z_W <- sim_data$Z_W
            Z_W <- Z_W[Z_W[, "Z_W"] ==1,]
            A[Z_W[,c(1,2)]] <- Z_W[,3]
            A[Z_W[,c(2,1)]] <- Z_W[,3]
            test <- c(test, all.equal(A, W))
            
          } else {
            test <- c(test, isSymmetric(A) == F & isSymmetric(W) == F)
            
            Z_W <- sim_data$Z_W
            Z_W <- Z_W[Z_W[, "Z_W"] ==2,]
            
            W[Z_W[,c(1,2)]] <- 0
            
            Z_W <- sim_data$Z_W
            Z_W <- Z_W[Z_W[, "Z_W"] ==1,]
            A[Z_W[,c(1,2)]] <- Z_W[,3]
            test <- c(test, all.equal(A, W))
          }
        } else {
          
          test <- c(test, all.equal(sim_data$A, (sim_data$W>0.0)*1.0))
        }
        
        
        if (family == "poisson"){
          
          test <- c(test,all(sim_data$W@x %% 1.0 == 0 & sim_data$W@x >0))
          
        } else {
          
          test <- c(test,all(sim_data$W@x >0))
          
        } 
        
        if (noise_weights_prob > 0){
          if (model !="RSR"){
            
            test <- c(test,
                      all.equal(sum(sim_data$A)/(nrow(sim_data$A)*(nrow(sim_data$A)-1)) + (2*table( sim_data$Z_W[, "Z_W"])[2])/(nrow(sim_data$A)*(nrow(sim_data$A)-1)),
                                sum(sim_data$W>0.0)/(nrow(sim_data$W)*(nrow(sim_data$W)-1)), check.attributes=F))
            
          } else {
            
            test <- c(test,
                      all.equal(sum(sim_data$A)/(nrow(sim_data$A)*(nrow(sim_data$A)-1)) + (table( sim_data$Z_W[, "Z_W"])[2])/(nrow(sim_data$A)*(nrow(sim_data$A)-1)),
                                sum(sim_data$W>0.0)/(nrow(sim_data$W)*(nrow(sim_data$W)-1)), check.attributes=F))
          }
        } else {
          test <- c(test,
                    all.equal(sum(sim_data$A)/(nrow(sim_data$A)*(nrow(sim_data$A)-1)),
                              sum(sim_data$W>0.0)/(nrow(sim_data$W)*(nrow(sim_data$W)-1)), check.attributes=F))
        }
        main <- c(main, all(test))
        
      }
    }
  }
  
  expect_true(all(main))
  
  
})

