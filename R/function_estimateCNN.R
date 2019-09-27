#' Estimating non negative C
#'
#' Given a reference matrix X, a matrix of bulks Y and a g-vector, "estimate_nn_c" finds
#' the solution of \deqn{arg min || diag(g) (Y - XC) ||_2} with non negative constraint over C,
#' using the package 'nnls' (non-negative least squares)
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
#'
#' @import nnls
estimate_nn_c <- function(X.matrix = NA, new.data, Gamma){
  estimates <- matrix(NA, nrow = ncol(X.matrix), ncol = ncol(new.data))
  rownames(estimates) <- colnames(X.matrix)
  colnames(estimates) <- colnames(new.data)
  g.X.matrix <- Gamma %*% X.matrix
  for(l.mix in colnames(new.data)){
    g.new.data <- Gamma %*% new.data[, l.mix]
    estimates.l.mix <- nnls::nnls(g.X.matrix, g.new.data)
    estimates[, l.mix] <- estimates.l.mix$x
  }
  return(estimates)
}
