#' Estimating C
#'
#' Given a reference matrix X, a matrix of bulks and a g-vector, est.cs finds the solution of
#' \deqn{arg min || diag(g) (Y - XC) ||_2}
#' over all C using direct analytical solution
#'
#' @param X numeric matrix with cells as columns, and features as rows. Reference matrix X of the DTD problem.
#' @param Y numeric matrix with samples as columns, and features as rows.
#' @param gamma.vec numeric vector with length of nrow(X). In the equation above the gamma.vec is denotated as g.
#'
#' @return numeric matrix with ncol(X) rows, and ncol(Y) columns
#'
#' @export
#'
#' @examples
#' library(DTD)
#'
#' # simulate random data:
#' random.data <- generate.random.data(nTypes = 5,
#'                                    nSamples.perType = 1,
#'                                    nFeatures = 100)
#'
#' # simulate a true c
#' # (this is not used by the est.cs function, it is only used to show the result!)
#' true.c <- rnorm(n = ncol(random.data), mean = 1, sd = 0.5)
#'
#' # calculate bulk y = Xc * some_error
#' bulk <- random.data %*% true.c * rnorm(n = nrow(random.data), mean = 1, sd = 0.01)
#'
#' #estimate c
#' estimated.c <- est.cs(X = random.data,
#'                       Y = bulk,
#'                       gamma.vec = rep(1, nrow(random.data)))
#'
#' # visualize that the estimated c are close to the true c
#' plot(true.c, estimated.c)
est.cs <- function(X, Y, gamma.vec){
  Gamma <- diag(gamma.vec)
  temp.sol <- solve(t(X)%*%Gamma%*%X)%*%t(X)%*%Gamma
  sol <- temp.sol%*%Y
  return(sol)
}
