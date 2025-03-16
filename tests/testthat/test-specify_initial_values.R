
## Testing specify_initial_values

test_that("specify_initial_values works", {
  
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
   beta0 <- -1
   sim_data <- JANE::sim_A(N = 100L, 
                           model = "RSR",
                           mus = mus, 
                           omegas = omegas, 
                           p = p, 
                           params_LR = list(beta0 = beta0), 
                           remove_isolates = TRUE)
   
   # Specify starting values
   D <- 3L
   K <- 5L
   N <- nrow(sim_data$A)
   n_interior_knots <- 5L
   
   U <- matrix(stats::rnorm(N*D), nrow = N, ncol = D)
   omegas <- stats::rWishart(n = K, df = D+1, Sigma = diag(D))
   mus <- matrix(stats::rnorm(K*D), nrow = K, ncol = D)
   p <- extraDistr::rdirichlet(n = 1, rep(3,K))[1,]
   Z <-  extraDistr::rdirichlet(n = N, alpha = rep(1, K))
   beta <- stats::rnorm(n = 1 + 2*(1 + n_interior_knots))
   
   expect_no_error({
     my_starting_values <- JANE::specify_initial_values(A = sim_data$A,
                                                        D = D,
                                                        K = K,
                                                        model = "RSR",
                                                        n_interior_knots = n_interior_knots,
                                                        U = U,
                                                        omegas = omegas, 
                                                        mus = mus, 
                                                        p = p, 
                                                        Z = Z,
                                                        beta = beta)         
     
     # Run JANE using my_starting_values (no need to specify D and K as function will 
     # determine those values from my_starting_values)
     res <- JANE::JANE(A = sim_data$A,
                       initialization = my_starting_values,
                       model = "RSR")
   })
   
})
