#' Identity
#'
#' @param x, any R object
#'
#' @return x
#' @export
#'
#' @examples
#' identity(1:10)
identity <- function(x){
  return(x)
}

#' Norm a vector to ||x||_2 = 1
#'
#' @param x numeric vector
#'
#' @return numeric vector with same length as x
#' @export
#'
#' @examples
#' set.seed(1)
#' same.values <- rnorm(n=10)
#' norm(same.values, type = "2")
#' normed.values <- n2normed(same.values)
#' norm(normed.values, type = "2")
#'
n2normed <- function(x){
  n2 <- norm(x, type = "2")
  tmp <- x/n2
  return(tmp)
}
