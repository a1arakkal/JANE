#' @useDynLib JANE  
log_Q <- function(A, U, mus, omegas, prob_matrix, beta, p, a, b, c, G, nu, e, f, X, n_control, model) {
  .Call('_JANE_log_Q', PACKAGE = 'JANE', A, U, mus, omegas, prob_matrix, beta, p, a, b, c, G, nu, e, f, X, n_control, model)
}

#' @useDynLib JANE  
log_Q_RE <- function(A, U, mus, omegas, prob_matrix, beta, p, a, b, c, G, nu, e, f, X, model, n_control) {
  .Call('_JANE_log_Q_RE', PACKAGE = 'JANE', A, U, mus, omegas, prob_matrix, beta, p, a, b, c, G, nu, e, f, X, model, n_control)
}