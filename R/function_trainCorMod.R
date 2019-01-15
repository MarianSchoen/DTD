#' Title
#'
#' @param tweak
#' @param X.matrix
#' @param train.list
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
train_correlation_model <- function(tweak, X.matrix, train.data.list, ...){
  # safety checks:
  if(!all(c("mixtures", "quantities") %in% names(train.data.list))){
    stop("In train_correlation_model: train.list must include the entries 'mixtures' and 'quantites'")
  }
  if(nrow(X.matrix) != nrow(train.data.list$mixtures) || ncol(X.matrix) != nrow(train.data.list$quantities)){
    stop("In train_correlation_model: Dimension of 'X.matrix' does not fit 'train.list'")
  }
  # define wrapper functions for gradient and correlation evaluation
  DTD.grad.wrapper <- function(tweak,
                               X = X.matrix,
                               train.list = train.data.list){
    Y <- train.list$mixtures
    C <- train.list$quantities
    grad <- gradient_cor_trace(X = X, Y = Y, C = C, tweak = tweak)
    return(grad)
  }
  DTD.evCor.wrapper <- function(tweak,
                                X = X.matrix,
                                train.list = train.data.list){
    Y <- train.list$mixtures
    C <- train.list$quantities
    loss <- evaluate_cor(X = X, Y = Y, C = C, tweak = tweak)/ncol(X)
    return(loss)
  }
  catch <- DTD_cv_lambda(tweak.start = tweak,
                         cv.verbose = TRUE,
                         train.list = train.data.list,
                         F.GRAD.FUN = DTD.grad.wrapper,
                         EVAL.FUN = DTD.evCor.wrapper,
                         ...
  )
  return(catch)
}
