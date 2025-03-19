
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
