
## Testing plot.JANE

test_that("plot.JANE works", {

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
  res <- JANE::JANE(A = sim_data$A,
                    D = 2L,
                    K = 3L,
                    initialization = "GNN", 
                    model = "NDH",
                    case_control = FALSE,
                    DA_type = "none")
 
  # plot trace plot
  expect_no_error( plot(res, type = "trace_plot") )
                     
  # plot network
  expect_no_error( plot(res) )
  
  # plot network - misclassified
  expect_no_error( plot(res, type = "misclassified", true_labels = apply(sim_data$Z_U, 1, which.max)) )
  
  # plot network - uncertainty and swap axes
  expect_no_error( plot(res, type = "uncertainty", swap_axes = TRUE) )
  
  # plot network - but only show contours of MVNs
  expect_no_error( plot(res, swap_axes = TRUE, alpha_edge = 0, alpha_node = 0) )
  
  # plot using starting values of EM algorithm
  expect_no_error( plot(res, initial_values = TRUE) )
  
})

## Testing summary.JANE

test_that("summary.JANE works", {
  
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
  res <- JANE::JANE(A = sim_data$A,
                    D = 2L,
                    K = 3L,
                    initialization = "GNN", 
                    model = "NDH",
                    case_control = FALSE,
                    DA_type = "none")
                    
  # Summarize fit 
  expect_no_error( summary(res) )
  
  # Summarize fit and compare to true cluster labels
  expect_no_error( summary(res, true_labels = apply(sim_data$Z_U, 1, which.max)) )
  
  # Summarize fit using starting values of EM algorithm
  expect_no_error( summary(res, initial_values = TRUE) )
  
})
