#' Nesterov Factor 2/(k+1)
#'
#' The FISTA algorithm extrapolates/correctes after the gradient step. The "nesterov_factor" function
#' returns how far the FISTA will exptrapolate in a given step "k"
#'
#' @param k integer, iteration step of the FISTA algorithm
#'
#' @return integer
#'
#' @examples
#' nesterov_factor(2)
nesterov_factor <- function(k) {
  test <- test_integer(test.value = k,
                       output.info = c("nesterov_factor", "k"),
                       min = 1,
                       max = Inf)
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
#'
#' @examples
#' set.seed(2202)
#' vec.with.neg <- rnorm(10)
#' range(vec.with.neg)
#' vec.without.neg <- positive_subspace_pmax(vec.with.neg)
#' range(vec.without.neg)
positive_subspace_pmax <- function(x) {
  test <- test_tweak_vec(tweak.vec = x,
                         output.info = c("positive_subspace_pmax", "x"))
  return(pmax(x, 0))
}
