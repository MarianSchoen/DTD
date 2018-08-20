#' Gradient of \eqn{L(g) = \sum cor(true_C, estimated_C(g) )}
#'
#' This function returns the value of the gradient of our Loss-function L.\cr
#' For the mathematical theory see Goertler et al 2018.
#' For examples see wrapper in \code{\link{descent_generalized_fista}}
#'
#' @param X numeric matrix with features as rows, and reference samples as columns
#' @param Y numeric matrix with features as rows, and samples as columns
#' @param C numeric matrix with reference samples as rows, and samples as columns
#' @param tweak numeric vector with length of nrow(X).
#' In the Loss function above tweak is named "g"
#'
#' @export
#'
#' @import matrixStats
#' @import Matrix
#'
#' @return numeric list, same length as "tweak"

Trace.H.gradient <- function(X, Y, C, tweak){

  Gamma <- Matrix::Matrix(diag(tweak))
  estimates.cs <- est.cs(X, Y, tweak)
  hat.sigmas <- matrixStats::rowSds(estimates.cs, na.rm = T)
  sigmas <- matrixStats::rowSds(C, na.rm = T)
  mean.hat.cs <- Matrix::rowMeans(estimates.cs)
  mean.cs <- Matrix::rowMeans(C)
  N <- ncol(C)
  cov.cs.hat.cs <- NULL

  for(j in 1:nrow(C)){
    cov.cs.hat.cs <- c(cov.cs.hat.cs, stats::cov(estimates.cs[j,],C[j,]) )
  }

  A <- matrix(NA, nrow=ncol(X), ncol=ncol(C))
  for(j in 1:ncol(X)){
    faktor <- 1/(hat.sigmas[j] * sigmas[j])
    cov.part <- cov.cs.hat.cs[j] / (N * hat.sigmas[j]^2)
    for(k in 1:ncol(C)){
      A[j, k] <- faktor * (cov.part * (estimates.cs[j,k] - mean.hat.cs[j]) - ((C[j,k] - mean.cs[j])/N))
    }
  }

  alpha <- solve(t(X)%*%Gamma%*%X)%*%t(X)
  beta <- (diag(1, nrow=nrow(Y), ncol=nrow(Y)) - X %*% alpha %*% Gamma) %*% Y

  B <- as.matrix(beta %*% t(A) %*% alpha)

  d.elems <- diag(B)

  return(d.elems)
}

