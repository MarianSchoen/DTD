#' Train a DTD model based on correlation loss function
#'
#' Loss-function learning Digital Tissue Deconvolution (DTD) adjustes a deconvolution model to
#' a biological context. 'train_deconvolution_model' is the main function of the DTD package.
#' As input it takes the reference matrix X, a list of training data and a start vector 'tweak'.
#' Then, it iteratively finds that vector 'g' that best deconvolutes based on the loss fucntion:
#' \deqn{L(g) = - \sum cor(C_{j,.} \widehat(C_{j,.}(g))) + \lambda ||g||_1}
#' The 'train_deconvolution_model' function calls the cross validation function \code{\link{DTD_cv_lambda}}
#' to find the optimal lambda. Afterwards it optimizes a model on the complete dataset.
#'
#' This function works as a wrapper for the correlation loss function and its gradient.
#' It provides a workaround for digital tissue deconvolution, such that the user only
#' has to provide a starting vector (in literature "g", in code 'tweak'), a reference matrix X,
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
#' expression profile of a training mixture. This data is used to train the deconvolution model
#' @param ... parameters passed to DTD_cv_lambda
#' @param test.data.list list, with two entries, same structure as train.data.list:
#' 'quantities', numeric matrix with as many rows as columns in X. Each column of the 'quantities' matrix
#'  holds the quantities of the cell types of X in the test mixtures.
#' 'mixtures', numeric matrix with as many rows as rows in X.matrix. Each column of 'mixtures' is the
#'  expression profile of a test mixture. After learning, this data is used to validate the model
#' @param estimate.c.type string, either "non_negative", or "direct". Indicates how the algorithm finds the solution of
#' \eqn{arg min_C ||diag(g)(Y - XC)||_2}. If estimate.c.type is set to "direct" there is no regularization
#' (see \code{\link{estimate_c}}),
#' if estimate.c.type is set to "non_negative" the estimates "C" must not be negative (non-negative least squares) (see (see \code{\link{estimate_nn_c}}))
#' @param useImplementation string, either "R" or "cxx". chooses between the R reference implementation and the faster c++ implementation.
#'
#' @return list, including 4 entries:
#' \itemize{
#'     \item cv.obj' (see \code{\link{DTD_cv_lambda}})
#'     \item 'best.model' (see \code{\link{DTD_cv_lambda}})
#'     \item 'reference.X'
#'     \item 'pics' (see `browseVignettes("DTD")`)
#' }
#' @export
train_deconvolution_model <- function(tweak,
                                    X.matrix,
                                    train.data.list,
                                    test.data.list = NULL,
                                    estimate.c.type,
                                    useImplementation = "cxx",
                                    ...){

  if(length(tweak) != nrow(X.matrix)){
    stop("In train_deconvolution_model: 'tweak' does not fit 'X.matrix'. 'length(tweak)' must be 'nrow(X.matrix)'")
  }

  if(is.null(names(tweak))){
    warning("Name the initial tweak.start vector.
            The naming will be kept during learning.
            Therefore, it is easier to keep track which 'g'/'tweak' belongs to which feature")
  }else{
    if(all(rownames(X.matrix) %in% names(tweak))){
      tweak <- tweak[rownames(X.matrix)]
    }else{
      stop("In train_deconvolution_model: There are features within 'X.matrix' where no entry in 'tweak' can be found")
    }
  }
  test <- test_tweak_vec(tweak.vec = tweak,
                        output.info = c("train_deconvolution_model", "tweak"))
  # end -> tweak

  # safety check for X.matrix:
  if(!(is.matrix(X.matrix) && !any(is.na(X.matrix)))){
    stop("In train_deconvolution_model: X.matrix must be a matrix without any NA.")
  }
  # end -> X.matrix

  # safety check for train.data.list:
  if(is.list(train.data.list) && length(train.data.list) == 2){
    if(!all(c("quantities", "mixtures") %in%  names(train.data.list))){
      stop("In train_deconvolution_model: entries of train.data.list must be named 'quantities' and 'mixtures'")
    }else{
      if(!is.matrix(train.data.list$mixtures)){
        stop("In train_deconvolution_model: 'train.data.list$mixtures' is not a matrix")
      }
      if(!is.matrix(train.data.list$quantities)){
        stop("In train_deconvolution_model: 'train.data.list$quantities' is not a matrix")
      }
    }
  }else{
    stop("In train_deconvolution_model: train.data.list must be provided as a list with two entries: 'quantities' and 'mixtures'")
  }
  # end -> train.data.list

  # safety check for test.data.list
  if(!is.null(test.data.list)){
    if(is.list(test.data.list) && length(test.data.list) == 2){
      if(!all(c("quantities", "mixtures") %in%  names(test.data.list))){
        stop("In train_deconvolution_model: entries of test.data.list must be named 'quantities' and 'mixtures'")
      }else{
        if(!is.matrix(test.data.list$mixtures)){
          stop("In train_deconvolution_model: 'test.data.list$mixtures' is not a matrix")
        }
        if(!is.matrix(test.data.list$quantities)){
          stop("In train_deconvolution_model: 'test.data.list$quantities' is not a matrix")
        }
      }
    }else{
      stop("In train_deconvolution_model: test.data.list must be provided as a list with two entries: 'quantities' and 'mixtures'")
    }
  }
  # end -> test.data.list

  # saftey check if X.matrix and test/train are compatible:
  if(nrow(train.data.list$mixtures) != nrow(X.matrix)){
    stop("In train_deconvolution_model: 'nrow(train.data.list$mixtures)' does not match 'nrow(X.matrix)'")
  }
  if(nrow(train.data.list$quantities) != ncol(X.matrix)){
    stop("In train_deconvolution_model: 'ncol(train.data.list$quantities)' does not match 'ncol(X.matrix)'")
  }

  if(!is.null(test.data.list)){
    if(nrow(test.data.list$mixtures) != nrow(X.matrix)){
      stop("In train_deconvolution_model: 'nrow(test.data.list$mixtures)' does not match 'nrow(X.matrix)'")
    }
    if(nrow(test.data.list$quantities) != ncol(X.matrix)){
      stop("In train_deconvolution_model: 'ncol(test.data.list$quantities)' does not match 'ncol(X.matrix)'")
    }
  }
  # safety checks if X.matrix, test and train are sorted in the same way:
  feature.names <- rownames(X.matrix)
  if(all(feature.names %in% rownames(train.data.list$mixtures))){
    train.data.list$mixtures <- train.data.list$mixtures[feature.names, ]
  }else{
    stop("In train_deconvolution_model: there are features in 'X.matrix' that are not in 'train.data.list$mixtures'")
  }
  type.names <- colnames(X.matrix)
  if(all(type.names %in% rownames(train.data.list$quantities))){
    train.data.list$quantities <- train.data.list$quantities[type.names, ]
  }else{
    stop("In train_deconvolution_model: there are cell types (=> columns) in 'X.matrix' that are not in 'train.data.list$quantities")
  }

  if(!is.null(test.data.list)){
    if(all(feature.names %in% rownames(test.data.list$mixtures))){
      test.data.list$mixtures <- test.data.list$mixtures[feature.names, ]
    }else{
      stop("In train_deconvolution_model: there are features (=> rows) in 'X.matrix' that are not in 'test.data.list$mixtures'")
    }
    if(all(type.names %in% rownames(test.data.list$quantities))){
      test.data.list$quantities <- test.data.list$quantities[type.names, ]
    }else{
      stop("In train_deconvolution_model: there are cell types (=> columns) in 'X.matrix' that are not in 'test.data.list$quantities")
    }
  }
  # end => compatible test

  if( useImplementation == "cxx" && estimate.c.type != "direct" ) {
    stop("only \"direct\" C estimation is implemented within cxx. use R for other ways to estimate C")
  }

  ESTIMATE.C.FUN <- test_c_type(test.value = estimate.c.type,
                                output.info = c("train_deconvolution_model", "estimate.c.type"))

  catch <- list()

  if( useImplementation == "R" ) {

    # define wrapper functions for gradient and correlation evaluation
    DTD.grad.wrapper <- function(tweak,
                                X = X.matrix,
                                train.list = train.data.list,
                                esti.c.type = estimate.c.type) {
      Y <- train.list$mixtures
      C <- train.list$quantities
      grad <- gradient_cor_trace(
        X = X,
        Y = Y,
        C = C,
        tweak = tweak,
        estimate.c.type = esti.c.type)
      return(grad)
    }

    DTD.evCor.wrapper <- function(tweak,
                                  X = X.matrix,
                                  train.list = train.data.list,
                                  esti.c.type = estimate.c.type) {
      Y <- train.list$mixtures
      C <- train.list$quantities
      loss <- evaluate_cor(
        X.matrix = X,
        new.data = Y,
        true.compositions = C,
        DTD.model = tweak,
        estimate.c.type = esti.c.type
      ) / ncol(X)
      return(loss)
    }

    catch <- DTD_cv_lambda_R(
      tweak.start = tweak,
      train.data.list = train.data.list,
      F.GRAD.FUN = DTD.grad.wrapper,
      EVAL.FUN = DTD.evCor.wrapper,
      ...
    )
  } else if( useImplementation == "cxx" || useImplementation == "cpp" ) {
    cat("using cxx implementation!")
    catch <- DTD_cv_lambda_cxx(
      tweak.start = tweak,
      train.data.list = train.data.list,
      X.matrix = X.matrix,
      ...
    )
  } else {
    stop(paste("cannot use implementation \"", useImplementation, "\": not implemented."))
  }
  catch$reference.X <- X.matrix

  pics <- vector(mode = "list")
  pics$cv <- DTD::ggplot_cv(DTD.model = catch)
  if(is.null(test.data.list)){
    pics$convergence <- DTD::ggplot_convergence(estimate.c.type = estimate.c.type,
                                                DTD.model = catch,
                                                X.matrix = X.matrix,
                                                test.data = NA
                                                )
  }else{
    pics$convergence <- DTD::ggplot_convergence(estimate.c.type = estimate.c.type,
                                                DTD.model = catch,
                                                X.matrix = X.matrix,
                                                test.data = test.data.list
    )
  }

  pics$path <- DTD::ggplot_gpath(catch)$gPath

  pics$histogram <- DTD::ggplot_ghistogram(DTD.model = catch)

  pics$Xheatmap <- ggplot_heatmap(DTD.model = catch,
                                  X.matrix = X.matrix)

  if(!is.null(test.data.list)){
    estimates <- ESTIMATE.C.FUN(new.data = test.data.list$mixtures,
                                DTD.model = catch)
    pics$true_vs_esti <- DTD::ggplot_true_vs_esti(DTD.model = catch,
                                                  test.data = test.data.list,
                                                  estimate.c.type = estimate.c.type,
                                                  X.matrix = X.matrix)
  }
  catch$pics <- pics
  return(catch)
}
