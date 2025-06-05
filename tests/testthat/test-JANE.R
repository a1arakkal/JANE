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
  
})
