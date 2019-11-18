#' Estimating non negative C
#'
#' Given a reference matrix X, a matrix of bulks Y and a g-vector,
#' "estimate_nn_c" finds the solution of
#' \deqn{arg min || diag(g) (Y - XC) ||_2} with non negative constraint over C,
#' using the package 'nnls' (non-negative least squares)
#'
#' @param X.matrix numeric matrix, with features/genes as rows,
#' and cell types as column. Each column of X.matrix is a reference
#' expression profile
#' @param new.data numeric matrix with samples as columns,
#' and features/genes as rows. In the formula above denoated as Y.
#' @param Gamma numeric matrix, with nrow(X) x nrow(X)
#'
#' @return numeric matrix with ncol(X) rows, and ncol(Y) columns
#'
#' @import nnls
estimate_nn_c <- function(X.matrix = NA, new.data, Gamma){
  if(!is.matrix(Gamma)){
    stop('in estimate_c_direct: Gamm must be a matrix')
  }
  estimates <- matrix(NA, nrow = ncol(X.matrix), ncol = ncol(new.data))
  rownames(estimates) <- colnames(X.matrix)
  colnames(estimates) <- colnames(new.data)
  g.X.matrix <- Gamma %*% X.matrix
  for(l.mix in 1:ncol(new.data)){
    g.new.data <- Gamma %*% new.data[, l.mix]
    estimates.l.mix <- nnls::nnls(g.X.matrix, g.new.data)
    estimates[, l.mix] <- estimates.l.mix$x
  }
  return(estimates)
}
