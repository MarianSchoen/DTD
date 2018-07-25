#' DOCUMENTATION
#'
#' @param x
#' @param lambda
#'
#' @return
#'
#' @examples
soft_thresholding <- function(x, lambda){
  if(length(lambda) == length(x)){
    x <- sign(x) * pmax(abs(x) - lambda, 0)
  }
  if(length(lambda) == 1){
    x <- sign(x) * pmax(abs(x) - rep(lambda, length(x)), 0)
  }
  return(x)
}
