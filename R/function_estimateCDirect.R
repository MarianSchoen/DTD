#' Estimating C
#'
#' Given a reference matrix X, a matrix of bulks Y and a g-vector, "estimate_c" finds
#' the solution of \deqn{arg min || diag(g) (Y - XC) ||_2} over all C using
#' direct analytical solution: (see Goertler et al. 2018)
#' \deqn{ C(g) = (X^T \Gamma X )^(-1) X^T \Gamma Y}
#' with \eqn{\Gamma} = diag(g)
#'
#' @param X.matrix numeric matrix with cells as columns, and features as rows.
#'  Reference matrix X of the DTD problem. X.matrix can be set to NA (default), if the DTD.model
#'  includes the reference matrix X (default for \code{\link{train_deconvolution_model}})
#' @param new.data numeric matrix with samples as columns, and features as rows.
#' In the formula above denoated as Y.
#' @param Gamma either a numeric vector with length of nrow(X),
#' or a list returned by \code{\link{train_deconvolution_model}}, \code{\link{DTD_cv_lambda}},
#' or\code{\link{descent_generalized_fista}}. In the equation above
#'   the DTD.model provides the vector g.
#'
#' @return numeric matrix with ncol(X) rows, and ncol(Y) columns
estimate_direct_c <- function(
  X.matrix, Gamma, new.data
){
  sol <- solve(t(X.matrix) %*% Gamma %*% X.matrix) %*% t(X.matrix) %*% Gamma %*% new.data
  return(sol)
}
