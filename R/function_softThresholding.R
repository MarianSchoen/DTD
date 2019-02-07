#' Soft thresholding
#'
#' implementation of the proximal operator for l1 penalty
#' (see e.g. Bubeck 2015)
#'
#' @param x list of numerics, x = g - step.size * gradient(g)
#' @param lambda float or list of floats, regularization parameter
#'
#' @export
#' @return vector of numerics, same length as x
#' @examples
#' set.seed(1)
#' soft_thresholding(x = runif(10) - rnorm(10), lambda = 0.2)
soft_thresholding <- function(x, lambda) {
  if (length(lambda) == length(x)) {
    ret <- sign(x) * pmax(abs(x) - lambda, 0)
  }
  if (length(lambda) == 1) {
    ret <- sign(x) * pmax(abs(x) - rep(lambda, length(x)), 0)
  }
  if (!exists("ret")) {
    stop("Error: soft_thresholding, length of lambda must either be 1 or same as x.")
  }
  return(ret)
}
