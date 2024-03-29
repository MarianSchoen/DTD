#' Estimating C
#'
#' Given a reference matrix X, a matrix of bulks Y and a g-vector, "estimate_c" finds
#' the solution of \deqn{arg min || diag(g) (Y - XC) ||_2} over all C using
#' direct analytical solution: (see Goertler et al. 2018)
#' \deqn{ C(g) = (X^T \Gamma X )^(-1) X^T \Gamma Y}
#' with \eqn{\Gamma} = diag(g)
#'
#' @param X.matrix numeric matrix, with features/genes as rows,
#' and cell types as column. Each column of X.matrix is a reference
#' expression profile.
#' @param new.data numeric matrix with samples as columns,
#' and features/genes as rows. In the formula above denoated as Y.
#' @param Gamma numeric matrix, with nrow(X) x nrow(X)
#'
#' @return numeric matrix with ncol(X) rows, and ncol(Y) columns
estimate_direct_c <- function(
  X.matrix, Gamma, new.data
){
  library(Matrix)
  
  if(!is.matrix(Gamma)){
    stop("in 'estimate_direct_c': Gamma must be a matrix")
  }
  
  xt.g.x <- t(X.matrix) %*% Gamma %*% X.matrix
  
  if(rankMatrix(x = xt.g.x) == 0){
    warning(
      paste0(
        "in 'estimate_direct_c': g/tweak/Gamma is to small, leading to a ", 
        "matrix of rank 0, for X^T G X. Returning C = 0."
        )
    )
    sol <- matrix(0, nrow = ncol(X.matrix), ncol = ncol(new.data))
    rownames(sol) <- colnames(X.matrix)
    return(sol)
  }
  
  sol <- chol2inv(chol(xt.g.x)) %*% t(X.matrix) %*% Gamma %*% new.data
  rownames(sol) <- colnames(X.matrix)
  return(sol)
}
