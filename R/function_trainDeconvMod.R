#' Train a DTD model based on correlation loss function
#'
#' Loss-function learning Digital Tissue Deconvolution (DTD) adapts a
#' deconvolution model to its biological context. 'train_deconvolution_model'
#' is the main function of the DTD package.\cr
#' As input it takes the reference matrix X, a list of training data and a
#' start vector 'tweak'. Then, it iteratively finds that vector 'g' that
#' deconvolutes best based on the loss fucntion:
#' \deqn{L(g) = - \sum cor(C_{j,.} \widehat C_{j,.}(g) ) + \lambda ||g||_1}
#' The 'train_deconvolution_model' function calls the cross validation function
#' \code{\link{DTD_cv_lambda_cxx}} (or \code{\link{DTD_cv_lambda_R}},
#' depending on 'use.implementation') to find the optimal lambda.
#' After the cross validation, it optimizes a model on the complete dataset
#' with the optimal \eqn{\lambda}.
#'
#' For an example see `browseVignettes("DTD")`
#'
#' @param tweak numeric vector with length of nrow(X).
#' In the Loss function above tweak is named "g"
#' Notice, the names of the vector will be kept, and are of use later on.
#' @param X.matrix numeric matrix, with features/genes as rows,
#' and cell types as column. Each column of X.matrix is a reference
#' expression profile
#' @param train.data.list list, with two entries, a numeric matrix each,
#' named 'mixtures' and 'quantities'
#' Within this list the train/test cross validation will be done.
#' (see Vignette `browseVignettes("DTD")` for details).
#' Generate 'train.data.list' using \code{\link{mix_samples}}
#' or \code{\link{mix_samples_with_jitter}}.
#' @param ... parameters passed to \code{\link{DTD_cv_lambda_cxx}}, or
#' \code{\link{DTD_cv_lambda_R}}
#' @param test.data.list list, with two entries, a numeric matrix each,
#' named 'mixtures' and 'quantities'
#' On this data, the trained model will be tested. Notice, this data is not
#'  shown to the optimization.
#' (see Vignette `browseVignettes("DTD")` for details).
#' Generate 'test.data.list' using \code{\link{mix_samples}}
#' or \code{\link{mix_samples_with_jitter}}.
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
#' @param use.implementation string, either "R" or "cxx".
#' Chooses between the R reference implementation and the faster c++
#' implementation.
#' Notice, if 'use.implementation' is set to "R" the cross validation
#' function \code{\link{DTD_cv_lambda_R}} is used.
#'
#' @return list, including 5 entries:
#' \itemize{
#'     \item cv.obj' (see \code{\link{DTD_cv_lambda_cxx}})
#'     \item 'best.model' (see \code{\link{DTD_cv_lambda_cxx}})
#'     \item 'reference.X'
#'     \item 'estimate.c.type'
#'     \item 'pics' (see `browseVignettes("DTD")`)
#' }
#' @export
train_deconvolution_model <- function(
  tweak,
  X.matrix,
  train.data.list,
  test.data.list = NULL,
  estimate.c.type,
  use.implementation = "cxx",
  ...
  ){
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

  # safety check: estimate.c.type
  test <- test_c_type(
    test.value = estimate.c.type,
    output.info = c("train_deconvolution_model", "estimate.c.type")
    )
  # end estimte.c.type

  if( use.implementation == "R" ) {

    # define wrapper functions for gradient and correlation evaluation
    DTD.grad.wrapper <- function(
      tweak,
      X = X.matrix,
      train.list = train.data.list,
      esti.c.type = estimate.c.type
      ) {
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

    DTD.evCor.wrapper <- function(
      tweak,
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
  } else if( use.implementation == "cxx" || use.implementation == "cpp" ) {
    catch <- DTD_cv_lambda_cxx(
      tweak.start = tweak,
      train.data.list = train.data.list,
      X.matrix = X.matrix,
      estimate.c.type = estimate.c.type,
      ...
    )
  } else {
    stop(paste("cannot use implementation \"", use.implementation, "\": not implemented."))
  }
  catch$reference.X <- X.matrix
  catch$estimate.c.type <- estimate.c.type

  pics <- vector(mode = "list")
  if("cv.obj" %in% names(catch)){
    pics$cv <- DTD::ggplot_cv(DTD.model = catch)
  }
  if(is.null(test.data.list)){
    pics$convergence <- DTD::ggplot_convergence(
      estimate.c.type = estimate.c.type
      , DTD.model = catch
      , X.matrix = X.matrix
      , test.data = NA
      )
  }else{
    pics$convergence <- DTD::ggplot_convergence(
      estimate.c.type = estimate.c.type
      , DTD.model = catch
      , X.matrix = X.matrix
      , test.data = test.data.list
      )
  }

  if("History" %in% names(catch)){
    pics$path <- DTD::ggplot_gpath(catch)$gPath
  }
  pics$histogram <- DTD::ggplot_ghistogram(DTD.model = catch)

  pics$Xheatmap <- DTD::ggplot_heatmap(DTD.model = catch,
                                  X.matrix = X.matrix)

  if(!is.null(test.data.list)){
    pics$true_vs_esti <- DTD::ggplot_true_vs_esti(
      DTD.model = catch
      , test.data = test.data.list
      , estimate.c.type = estimate.c.type
      , X.matrix = X.matrix)
  }
  catch$pics <- pics
  return(catch)
}
