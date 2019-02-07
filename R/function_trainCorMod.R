#' Train a DTD model based on correlation loss function
#'
#' This function works as a wrapper for the correlation loss function and its gradient.
#' It provides a workaround for digital tissue deconvolution, such that the user only
#' has to provide a starting vector (in literature g, in code tweak), a reference matrix X,
#' and a training list, including the training mixtures, and their cell type quantities.
#'
#' For an example see `browseVignettes("DTD")`
#'
#' @param tweak numeric vector, starting point for the FISTA optimization.
#' Notice, the names of the vector will be kept, and are of use later on.
#' @param X.matrix numeric matrix, with features/genes as rows, and cell types as column.
#' Each column of X.matrix is a reference expression profile
#' @param train.data.list list, with two entries: 'quantities', numeric matrix with as many rows as columns in X.
#' Each column of the 'quantities' matrix holds the quantities of the cell types of X in the training mixtures.
#' 'mixtures', numeric matrix with as many rows as rows in X.matrix. Each column of 'mixtures' is the
#' expression profile of a training mixture.
#' @param ... parameters passed to DTD_cv_lambda
#'
#' @return
#' @export
train_correlation_model <- function(tweak, X.matrix, train.data.list, ...) {
  # safety checks:
  if(any(is.na(tweak))){
    stop("train_correlation_model: Tweak vector includes NAs")
  }
  if (!all(c("mixtures", "quantities") %in% names(train.data.list))) {
    stop("In train_correlation_model: train.list must include the entries 'mixtures' and 'quantites'")
  }
  if (nrow(X.matrix) != nrow(train.data.list$mixtures) || ncol(X.matrix) != nrow(train.data.list$quantities)) {
    stop("In train_correlation_model: Dimension of 'X.matrix' does not fit 'train.list'")
  }
  # define wrapper functions for gradient and correlation evaluation
  DTD.grad.wrapper <- function(tweak,
                                 X = X.matrix,
                                 train.list = train.data.list) {
    Y <- train.list$mixtures
    C <- train.list$quantities
    grad <- gradient_cor_trace(X = X, Y = Y, C = C, tweak = tweak)
    return(grad)
  }
  DTD.evCor.wrapper <- function(tweak,
                                  X = X.matrix,
                                  train.list = train.data.list) {
    Y <- train.list$mixtures
    C <- train.list$quantities
    loss <- evaluate_cor(
      X.matrix = X,
      new.data = Y,
      true.compositions = C,
      DTD.model = tweak
    ) / ncol(X)
    return(loss)
  }
  catch <- DTD_cv_lambda(
    tweak.start = tweak,
    train.list = train.data.list,
    F.GRAD.FUN = DTD.grad.wrapper,
    EVAL.FUN = DTD.evCor.wrapper,
    ...
  )
  catch$reference.X <- X.matrix
  return(catch)
}
