
## Testing specify_priors

test_that("specify_priors works", {
  
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
                           model = "RS",
                           mus = mus, 
                           omegas = omegas, 
                           p = p, 
                           params_LR = list(beta0 = beta0),
                           remove_isolates = TRUE)
                           
                           
   # Specify prior hyperparameters
   D <- 3L
   K <- 5L
   n_interior_knots <- 5L
   
   a <- rep(1, D)
   b <- 3
   c <- 4
   G <- 10*diag(D)
   nu <- rep(2, K)
   e <- rep(0.5, 1 + (n_interior_knots + 1))
   f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
   
   expect_no_error({
     my_prior_hyperparameters <- specify_priors(D = D,
                                                K = K,
                                                model = "RS",
                                                n_interior_knots = n_interior_knots,
                                                a = a,
                                                b = b,
                                                c = c,
                                                G = G,
                                                nu = nu,
                                                e = e,
                                                f = f)
     
     # Run JANE on simulated data using supplied prior hyperparameters
     res <- JANE::JANE(A = sim_data$A,
                       D = D,
                       K = K,
                       initialization = "GNN",
                       model = "RS",
                       case_control = FALSE,
                       DA_type = "none",
                       control = list(priors = my_prior_hyperparameters))
     
     # Specify first set of prior hyperparameters
     D <- 3
     K <- 5
     n_interior_knots <- 5L
     a <- rep(1, D)
     b <- 3
     c <- 4
     G <- 10*diag(D)
     nu <- rep(2, K)
     e <- rep(0.5, 1 + (n_interior_knots + 1))
     f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
     
     my_prior_hyperparameters_1 <- specify_priors(D = D,
                                                  K = K,
                                                  model = "RS",
                                                  n_interior_knots = n_interior_knots,
                                                  a = a,
                                                  b = b,
                                                  c = c,
                                                  G = G,
                                                  nu = nu,
                                                  e = e,
                                                  f = f)
     
     # Specify second set of prior hyperparameters
     D <- 2
     K <- 3
     n_interior_knots <- 5L
     a <- rep(1, D)
     b <- 3
     c <- 4
     G <- 10*diag(D)
     nu <- rep(2, K)
     e <- rep(0.5, 1 + (n_interior_knots + 1))
     f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
     
     my_prior_hyperparameters_2 <- specify_priors(D = D,
                                                  K = K,
                                                  model = "RS",
                                                  n_interior_knots = n_interior_knots,
                                                  a = a,
                                                  b = b,
                                                  c = c,
                                                  G = G,
                                                  nu = nu,
                                                  e = e,
                                                  f = f)
     
     # Create nested list
     my_prior_hyperparameters <- list(my_prior_hyperparameters_1,
                                      my_prior_hyperparameters_2)
     
     # Run JANE on simulated data using supplied prior hyperparameters
     res <- JANE::JANE(A = sim_data$A,
                       D = D,
                       K = K,
                       initialization = "GNN",
                       model = "RS",
                       case_control = FALSE,
                       DA_type = "none",
                       control = list(priors = my_prior_hyperparameters))
     
     # Specify prior hyperparameters as unevaluated calls
     n_interior_knots <- 5L
     e <- rep(0.5, 1 + (n_interior_knots + 1))
     f <- diag(c(0.1, rep(0.5, n_interior_knots + 1)))
     my_prior_hyperparameters <- specify_priors(model = "RS",
                                                n_interior_knots = n_interior_knots,
                                                a = quote(rep(1, D)),
                                                b = b,
                                                c = quote(D + 1),
                                                G = quote(10*diag(D)),
                                                nu = quote(rep(2, K)),
                                                e = e,
                                                f = f)

     # Run JANE using supplied prior hyperparameters
     JANE(A = sim_data$A,
          D = 2:3,
          K = 2:3,
          initialization = "GNN",
          model = "RS",
          case_control = FALSE,
          noise_weights = FALSE,
          control = list(priors = my_prior_hyperparameters))
   })

   
})
