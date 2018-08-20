#' Soft thresholding
#'
#' implementation of the proximal operator for l1 penalty
#' (see Bubeck 2015)
#'
#' @param x list of numerics, x = g - step.size * gradient(g)
#' @param lambda float or list of floats, regularization parameter
#'
#' @export
#' @return
#'
soft_thresholding <- function(x, lambda){
  if(length(lambda) == length(x)){
    x <- sign(x) * pmax(abs(x) - lambda, 0)
  }
  if(length(lambda) == 1){
    x <- sign(x) * pmax(abs(x) - rep(lambda, length(x)), 0)
  }
  return(x)
}
