
#' @export
specify_initial_values <- function(A,
                                   D,
                                   K,
                                   model,
                                   n_interior_knots,
                                   U, # latent positions (NxD matrix)
                                   omegas, # cluster precisions (DxDxK array)
                                   mus, # cluster centers (KXD matrix)
                                   p, # mixture probabilities (vector of length K) 
                                   prob_matrix, # probability of node belonging to cluster k (NxD matrix)
                                   beta){ # beta parameters (vector of length 1 + (model =="RS")*(n_interior_knots + 1) +  (model =="RSR")*2*(n_interior_knots + 1))
  
  list_name <- list(U = U, 
                    omegas = omegas, 
                    mus = mus, 
                    p = p, 
                    prob_matrix = prob_matrix,
                    beta = beta)
  
  check_initial_values(list_name = list_name,
                       A = A,
                       K = K,
                       D = D,
                       n_interior_knots = n_interior_knots,
                       model = model)
  return(list_name)
  
}