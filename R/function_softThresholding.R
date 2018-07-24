#' DOCUMENTATION
#'
#' @param x
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
soft_thresholding <- function(x, lambda){
  if(length(lambda) == length(x)){
    x <- sign(x) * pmax(abs(x) - lambda, 0)
  }
  if(length(lambda) == 1){
    x <- sign(x) * pmax(abs(x) - rep(lambda, length(x)), 0)
  }
  ### adjusted, such that g_i <= 0
#  x[x < 0] <- 0
  return(x)
}
