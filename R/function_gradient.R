#' Gradient of \eqn{L(g) = \sum cor(true_C, estimated_C(g) )}
#'
#' This function returns the value of the gradient of our Loss-function L.\cr
#' For the mathematical theory see Goertler et al, 2018.
#' For examples see code of \code{\link{train_deconvolution_model}}
#'
#' @param X.matrix numeric matrix, with features/genes as rows,
#' and cell types as column. Each column of X.matrix is a reference
#' expression profile
#' @param Y numeric matrix with samples as columns,
#' and features/genes as rows. Holding bulk gene expressions
#' @param C numeric matrix with reference samples as rows,
#' and samples as columns. Holding cellular composionts
#' @param tweak numeric vector with length of nrow(X).
#' In the Loss function above tweak is named "g"
#' @param estimate.c.type string, either "non_negative", or "direct".
#' Indicates how the algorithm finds the solution of
#' \eqn{arg min_C ||diag(g)(Y - XC)||_2}.
#' \itemize{
#'    \item If 'estimate.c.type' is set to "direct",
#'  there is no regularization (see \code{\link{estimate_c}}),
#'    \item if 'estimate.c.type' is set to "non_negative",
#'  the estimates "C" must not be negative (non-negative least squares)
#' (see (see \code{\link{estimate_nn_c}}))
#' }
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

  Gamma <- Matrix::Matrix(diag(tweak))
  estimates.cs <- estimate_c(
    X.matrix = X
    , new.data =Y
    , DTD.model = tweak
    , estimate.c.type = estimate.c.type
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
  
  alpha <- chol2inv(x = chol(x = as.matrix(t(X) %*% Gamma %*% X))) %*% t(X)
  beta <- (diag(1, nrow = nrow(Y), ncol = nrow(Y)) - X %*% alpha %*% Gamma) %*% Y

  B <- as.matrix(beta %*% t(A) %*% alpha)

  d.elems <- diag(B)

  # due to the constraint that all g have to be positive:
  d.elems[d.elems > 0] <- 0
  return(d.elems)
}
