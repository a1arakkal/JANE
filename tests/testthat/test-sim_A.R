
## Testing sim_A

test_that("sim_A works", {
  
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
  expect_no_error(
    JANE::sim_A(N = 100L, 
                model = "NDH",
                mus = mus, 
                omegas = omegas, 
                p = p, 
                beta0 = beta0, 
                remove_isolates = TRUE) 
  )
  
})
