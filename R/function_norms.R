#' Identity
#'
#' @param x, any R object
#'
#' @return x
#'
#' @examples
#' DTD::identity(1:10)
identity <- function(x) {
  return(x)
}

#' Norm a vector such that ||x||_2 = 1 is equal to the length of
#'
#' @param x numeric vector
#' @param to numeric, to this value, the norm gets scaled
#'
#' @return numeric vector with same length as x
#'
#' @examples
#' set.seed(1)
#' same.values <- rnorm(n = 10)
#' print(norm(same.values, type = "2"))
#' normed.values <- n2normed(same.values)
#' print(norm(normed.values, type = "2"))
n2normed <- function(x, to=NA) {
  if( is.na(to) ){
    to <- length(x)
  }
  n2 <- norm(matrix(x, ncol = 1), type = "2")
  if(n2 == 0){
    ret <- x
  }else{
    ret <- (to * x) / n2
  }
  return(ret)
}

#' Norm a vector such that ||x||_1 is equal to the length of x
#'
#' @param x numeric vector
#'
#' @return numeric vector with same length as x
#' @param to numeric, to this value, the norm gets scaled
#'
#' @examples
#' set.seed(1)
#' same.values <- rnorm(n = 10)
#' print(norm(as.matrix(same.values), type = "O"))
#' normed.values <- n1normed(same.values)
#' print(norm(as.matrix(normed.values), type = "O"))
n1normed <- function(x, to=NA) {
  if( is.na(to) ){
    to <- length(x)
  }  
  n1 <- norm(matrix(x, ncol = 1), type = "O")
  if(n1 == 0){
    ret <- x
  }else{
    ret <- (to * x) / n1
  }
  return(ret)
}
