#' Nesterov Faktor 2/(k+1)
#'
#' The FISTA algorithm extrapolates/correctes after the gradient step. The "nesterov_faktor" function
#' returns how far the FISTA will exptrapolate in a given step "k"
#'
#' @param k integer, iteration step of the FISTA algorithm
#'
#' @return integer
#'
#' @export
#'
#' @examples
#' nesterov_faktor(2)
nesterov_faktor <- function(k) {
  return(2 / (k + 1))
}
#' Wrapper for pmax
#'
#' positive_subspace_pmax expects a vector x, and applies pmax(x, 0) on it.
#' It is needed as a wrapper within the \code{\link{descent_generalized_fista}} implementation.
#'
#' @param x numeric vector
#'
#' @return pmax(x,0), numeric vector with same length as x, but without any negative entries
#' @export
#'
#' @examples
#' set.seed(2202)
#' vec.with.neg <- rnorm(10)
#' range(vec.with.neg)
#' vec.without.neg <- positive_subspace_pmax(vec.with.neg)
#' range(vec.without.neg)
positive_subspace_pmax <- function(x) {
  return(pmax(x, 0))
}
