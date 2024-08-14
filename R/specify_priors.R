
#' @export
specify_priors <- function(D, 
                           K,
                           model,
                           n_interior_knots,
                           a, # prior on mean of mus (vector of length D)
                           b, # prior on precision of mus (scalar)
                           c, # prior on df for omegas (scalar)
                           G, # prior on sigma for omegas (DxD matrix)
                           nu, # prior for dirichlet vector of length k
                           e, # prior mean on beta (vector of length 1 + (model =="RS")*(n_interior_knots + 1) +  (model =="RSR")*2*(n_interior_knots + 1))
                           f){ # prior precision on beta (square matrix of dim 1 + (model =="RS")*(n_interior_knots + 1) +  (model =="RSR")*2*(n_interior_knots + 1))
  
  priors <- list(
    a = t(a),
    b = b,
    c = c,
    G = G,
    nu = nu,
    e = e,
    f = f
  )
  
  check_priors(priors = priors,
               D = D,
               K = K,
               n_interior_knots = n_interior_knots,
               model = model)
  
  return(priors)
}