#' Gradient of \eqn{L(g) = \sum cor(true_C, estimated_C(g) )}
#'
#' This function returns the value of the gradient of our Loss-function L.\cr
#' For the mathematical theory see Goertler et al, 2018.
#' For examples see code of \code{\link{train_deconvolution_model}}
#'
#' @param X numeric matrix with features as rows, and reference samples as columns
#' @param Y numeric matrix with features as rows, and samples as columns
#' @param C numeric matrix with reference samples as rows, and samples as columns
#' @param tweak numeric vector with length of nrow(X).
#' In the Loss function above tweak is named "g"
#' @param estimate.c.type string, either "non_negative", or "direct".
#' Indicates how the algorithm finds the solution of
#'
#' @export
#'
#' @import matrixStats
#' @import Matrix
#'
#' @return numeric list, same length as "tweak"

gradient_cor_trace <- function(X, Y, C, tweak, estimate.c.type) {
  # safety check: tweak
  test <- test_tweak_vec(tweak.vec = tweak,
                         output.info = c("gradient_cor_trace", "tweak"))
  # end -> tweak
  # safety check X
  if(!any(is.numeric(X))){
    stop("In gradient_cor_trace: X is not numeric")
  }
  # end --> X
  # safety check C
  if(!any(is.numeric(C))){
    stop("In gradient_cor_trace: C is not numeric")
  }
  # end --> C
  # safety check Y
  if(!any(is.numeric(Y))){
    stop("In gradient_cor_trace: Y is not numeric")
  }
  # end --> Y

  if(estimate.c.type %in% c("non_negative", "direct")){
    if(estimate.c.type == "non_negative"){
      ESTIMATE.C.FUN <- estimate_nn_c
    }else{
      ESTIMATE.C.FUN <- estimate_c
    }
  }else{
    stop("In train_correlation_model: estimate.c.type does not match 'non_negative' or 'direct.")
  }

  Gamma <- Matrix::Matrix(diag(tweak))
  estimates.cs <- ESTIMATE.C.FUN(
    X.matrix = X,
    new.data = Y,
    DTD.model = tweak
  )
  hat.sigmas <- matrixStats::rowSds(estimates.cs, na.rm = T)
  sigmas <- matrixStats::rowSds(C, na.rm = T)
  mean.hat.cs <- Matrix::rowMeans(estimates.cs)
  mean.cs <- Matrix::rowMeans(C)
  N <- ncol(C)
  cov.cs.hat.cs <- NULL

  for (j in 1:nrow(C)) {
    cov.cs.hat.cs <- c(cov.cs.hat.cs, stats::cov(estimates.cs[j, ], C[j, ]))
  }

  A <- matrix(NA, nrow = ncol(X), ncol = ncol(C))
  for (j in 1:ncol(X)) {
    faktor <- 1 / (hat.sigmas[j] * sigmas[j])
    cov.part <- cov.cs.hat.cs[j] / (N * hat.sigmas[j]^2)
    for (k in 1:ncol(C)) {
      A[j, k] <- faktor * (cov.part * (estimates.cs[j, k] - mean.hat.cs[j]) - ((C[j, k] - mean.cs[j]) / N))
    }
  }

  alpha <- solve(t(X) %*% Gamma %*% X) %*% t(X)
  beta <- (diag(1, nrow = nrow(Y), ncol = nrow(Y)) - X %*% alpha %*% Gamma) %*% Y

  B <- as.matrix(beta %*% t(A) %*% alpha)

  d.elems <- diag(B)

  # due to the constraint that all g have to be positive:
  d.elems[d.elems > 0] <- 0
  return(d.elems)
}
