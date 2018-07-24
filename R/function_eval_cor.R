#' Evaluate correlation
#' 
#' The loss-function learning digital tissue deconvolution finds a vector g which minimizes the Loss function L\cr
#' \deqn{L(g) = \sum cor(true_C, estimatd_C(g))}
#' For a given vector g the evaluate_cor function returns the value of the loss function. 
#'
#' @param X numeric matrix with cells as columns, and features as rows. Reference matrix X of the DTD problem. 
#' @param Y numeric matrix with samples as columns, and features as rows. Each sample in Y is a bulk measurement, 
#' for which the quantity of the cells in X are known (and saved in C)
#' @param C numeric matrix with cells as rows, and mixtures as columns. 
#' Each row of C holds the distribution of the cell over all mixtures. 
#' @param tweak numeric vector with length of nrow(X). 
#'
#' @return
#' @export
#'
#' @examples
evaluate_cor <- function(X, Y, C, tweak){
  esti.cs <- est.cs(X, Y, tweak)
  tmp <- rep(NA, nrow(C))
  for(l1 in 1:nrow(C)){
    tmp[l1] <- cor(C[l1, ], esti.cs[l1, ])
  }
  return(1 - mean(tmp))
}