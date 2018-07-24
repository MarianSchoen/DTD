library(matrixStats)
library(Matrix)

#' Title
#'
#' @param X 
#' @param Y 
#' @param C 
#' @param gamma.vec 
#'
#' @return numeric list, same length as "gamma.vec" 
#' @export
#'
#' @examples Trace.H.gradient(tweak)
Trace.H.gradient <- function(X = x_mat, Y = y_mat, C = c_mat, gamma.vec = v_vec){

  Gamma <- Matrix(diag(gamma.vec))
  estimates.cs <- est.cs(X, Y, gamma.vec)
  hat.sigmas <- rowSds(estimates.cs, na.rm = T)
  sigmas <- rowSds(C, na.rm = T)
  mean.hat.cs <- rowMeans(estimates.cs)
  mean.cs <- rowMeans(C)
  N <- ncol(C)
  cov.cs.hat.cs <- NULL
  for(j in 1:nrow(C)){
    cov.cs.hat.cs <- c(cov.cs.hat.cs, cov(estimates.cs[j,],C[j,]) )
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
  
  B <- beta %*% t(A) %*% alpha
  d.elems <- diag(B)
  
  return(d.elems)
}
