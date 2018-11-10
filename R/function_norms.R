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

#' Norm a vector such that ||x||_2 = 1 is equal to the length of
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
  n2 <- norm(matrix(x, ncol=1), type = "2")
  tmp <- (length(x) * x)/n2
  return(tmp)
}

#' Norm a vector such that ||x||_1 is equal to the length of x
#'
#' @param x numeric vector
#'
#' @return numeric vector with same length as x
#' @export
#'
#' @examples
#' set.seed(1)
#' same.values <- rnorm(n=10)
#' norm(as.matrix(same.values), type = "O")
#' normed.values <- n1normed(same.values)
#' norm(as.matrix(normed.values), type = "O")
#'
n1normed <- function(x){
  n2 <- norm(matrix(x, ncol=1), type="O")
  tmp <- (length(x)*x)/n2
  return(tmp)
}

