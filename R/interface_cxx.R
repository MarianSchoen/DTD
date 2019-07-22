solve_fista_goertler <- function(model, lambda, maxiter) {
  # check input params...
  result <- .Call("_dtd_solve_fista_goertler", model, lambda, maxiter, PACKAGE="DTD")
  return(result)
}