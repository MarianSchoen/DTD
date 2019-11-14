#' Cross-validation for digital tissue deconvolution
#'
#' Our descent generalized FISTA implementation includes a l1 regularization term
#' (see \code{\link{train_deconvolution_model}}.
#' This function performs a k-fold cross validation to find the
#' best fitting regularization parameter.
#' For an example see `browseVignettes("DTD")`
#'
#' @param lambda.seq numeric vector or NULL: Over this series of lambdas the FISTA will be optimized.
#' If 'lambda.seq' is set to NULL, a generic series of lambdas - depending on the dimensions
#' of the training set -  will be generated. Default: NULL
#' @param tweak.start numeric vector, starting vector for the DTD algorithm.
#' @param n.folds integer, number of buckets in the cross validation. Defaults to 10
#' @param lambda.length integer, how many lambdas will be generated
#' (only used if lambda.seq is NULL). Defaults to 10
#' @param train.data.list list, that can be passed to the F.GRAD.FUN and EVAL.FUN.
#' Within this list the train/test cross validation will be done.
#' Notice, that the train.data.list must have an entry named "mixtues".
#' In this entry, the matrix containing the training samples (as columns)
#' and all features (as rows) must be present. (see Vignette for details)
#' @param F.GRAD.FUN gradient function, see \code{\link{descent_generalized_fista}}
#' @param EVAL.FUN evaluation function, see \code{\link{descent_generalized_fista}}
#' @param cv.verbose logical, should information about the cv process be printed to the screen? (Defaults to TRUE)
#' @param ... all parameters that are passed to the \code{\link{descent_generalized_fista}} function.
#' E.g. maxiter, tweak_vec etc ...
#' @param warm.start logical, should the solution of a previous model of the cross validation be used as a start
#'  in the next model. Defaults to TRUE. Notice, that the warm.start starts with the most unpenalized model.
#'
#' @return list of length 2.
#' \itemize{
#'    \item 'cv.obj', list of lists. DTD model for each lambda, and every folds.
#'    \item 'best.model', list. DTD model optimized on the complete data set with the best lambda from the cross validation.
#' }
#' A cross validation matrix as entry "cv.obj", and the model with minimal loss function
#' retrained on the complete dataset as "best.model
#' @export
#'
DTD_cv_lambda_R <- function(
  lambda.seq = NULL,
  tweak.start,
  n.folds = 5,
  lambda.length = 10,
  train.data.list,
  F.GRAD.FUN,
  EVAL.FUN,
  cv.verbose = TRUE,
  warm.start = FALSE,
  ...
  ) {

  DTD_cv_lambda_test_input_generic(
    lambda.seq = lambda.seq
    , tweak.start = tweak.start
    , n.folds = n.folds
    , lambda.length = lambda.length
    , train.data.list = train.data.list
    , cv.verbose = cv.verbose
    , warm.start = warm.start
    )

  # First, all possible training samples get assigned to a bucket:
  # extract Y:
  train.Y <- train.data.list$mixtures
  bucket.indicator <- make_buckets(
    train.Y = train.Y,
    folds = n.folds)

  if (!is.numeric(lambda.length)) {
    lambda.length <- 20
  }
  lambda.seq <- lambda_sequence(lambda.seq, lambda.length, train.Y)

  # after cross validation, a model is trained on the complete training
  # data using the best lambda of the cross validation.
  # The starting tweak for the end model must not be the result of a warm start,
  # therefore store it:
  tweak.start.end.model <- tweak.start

  # in older versions, I only kept the resulting loss for every fold and each lambda.
  # I think, for visualizing the cross validation, and comparing, it is better to return the complete models
  cv.object <- list()
  # Start of cross validation:
  for (lambda in lambda.seq) {
    if (cv.verbose) {
      pos <- which(lambda.seq == lambda) - 1
      cat("\ndoing lambda: ", lambda, ", completed ", pos, " of ", length(lambda.seq), ", ", 100 * pos / length(lambda.seq), "% \n")
    }
    if (cv.verbose) {
      cat("doing l.fold: ")
    }
    lambda.fold <- list()
    for (l.fold in 1:n.folds) {
      if (cv.verbose) {
        cat(l.fold, "\t")
      }

      # Split the complete training data into test and train:
      test.samples <- names(which(bucket.indicator == l.fold))
      train.samples <- names(which(bucket.indicator != l.fold))

      # reduce the train to only include the cv train samples ...
      tmp.train.list <- lapply(train.data.list, select.fun, samples = train.samples)

      # ... and reset the default values of the gradient and evaluation functions:
      tmp.grad.fun <- function(tmp.tweak, tmp.list = tmp.train.list) {
        return(F.GRAD.FUN(tmp.tweak, train.list = tmp.list))
      }
      tmp.eval.fun <- function(tmp.tweak, tmp.list = tmp.train.list) {
        return(EVAL.FUN(tmp.tweak, train.list = tmp.list))
      }
      # Now, try to train a model on the reduced training set:
      catch <- try(descent_generalized_fista(
        lambda = lambda,
        tweak.vec = tweak.start,
        F.GRAD.FUN = tmp.grad.fun,
        EVAL.FUN = tmp.eval.fun,
        ...
      ),
      silent = TRUE
      )
      # If the regularization parameter lambda is too big,
      # the fista algorithm does not find a model, and throws an error
      if (any(grepl(pattern = "Error", catch))) {
        lambda.fold[[as.character(l.fold)]] <- "could not build a model"
        next
      }
      # warm start, after learning a model, keep last tweak vec as start for next model:
      if (warm.start) {
        tweak.start <- catch$Tweak
      }
      # Evaluate the reached minimum on the test set:
      tmp.test.list <- lapply(train.data.list, select.fun, samples = test.samples)
      tmp.eval.fun.test <- function(tmp.tweak, tmp.list = tmp.test.list) {
        return(EVAL.FUN(tmp.tweak, train.list = tmp.list))
      }
      catch$cor.test <- tmp.eval.fun.test(catch$Tweak)
      lambda.fold[[as.character(l.fold)]] <- catch
    }
    cv.object[[as.character(lambda)]] <- lambda.fold
  }
  if (cv.verbose) {
    cat("\ncross validation completed, \n starting to build model on complete data, \n with  best lambda \n")
  }

  # after the cross validation, find the lambda with best evaluation score
  # pick the average mean per lambda:
  test.result.per.lambda <- unlist(
    lapply(
      cv.object
      , pick.mean.test.results.function
    )
  )

  if(all(is.na(test.result.per.lambda))){
    stop(
      "in 'DTD_cv_lambda_cxx': all lambdas crashed. \n
      Did you provide a custom lambda sequence? \n
      If yes: provide smaller lambdas \n
      If no: provide more features to the model,
      or provide a custom lambda sequence
      (via 'lambda.seq' in train_deconvolution_model(...). \n
      in 'train_deconvolution_model' set 'cv.verbose=TRUE' for logging information"
    )
  }

  mean.test.results <- test.result.per.lambda

  # and rebuild a model on the complete dataset:
  lmin.pos <- which.min(mean.test.results)
  lmin <- as.numeric(names(mean.test.results)[lmin.pos])

  bestModel <- descent_generalized_fista(
    lambda = lmin,
    tweak.vec = tweak.start.end.model,
    F.GRAD.FUN = F.GRAD.FUN,
    EVAL.FUN = EVAL.FUN,
    save.all.tweaks = TRUE,
    ...
  )

  # return the cv.object for plotting, and the model with best lambda
  ret <- list(cv.obj = cv.object, best.model = bestModel)
  return(ret)
}





#' Cross-validation for digital tissue deconvolution
#'
#' Our descent generalized FISTA implementation includes a l1 regularization term (see \code{\link{train_deconvolution_model}}.
#' This function performs a k-fold cross validation to find the best fitting regularization parameter.
#' For an example see `browseVignettes("DTD")`
#' This routine calls the cxx variant of DTD_cv_lambda. See that function for dynamic dispatch, based on the useImplementation parameter
#'
#' @param lambda.seq numeric vector or NULL: Over this series of lambdas the FISTA will be optimized.
#' If lambda.seq is set to NULL, a generic series of lambdas - depending on the dimensions
#' of the training set -  will be generated. Default: NULL
#' @param tweak.start numeric vector, starting vector for the DTD algorithm.
#' @param X.matrix reference matrix, X, each column i holds the gene expression of celltype i.
#' @param n.folds integer, number of buckets in the cross validation. Defaults to 10
#' @param lambda.length integer, how many lambdas will be generated (only used if lambda.seq is NULL). Defaults to 10
#' @param train.data.list list, that can be passed to the F.GRAD.FUN and EVAL.FUN.
#' Within this list the train/test cross validation will be done.
#' Notice, that the train.data.list must have an entry named "mixtues". In this entry, the matrix containing the
#' training samples (as columns) and all features (as rows) must be present. (see Vignette for details)
#' @param cv.verbose logical, should information about the cv process be printed to the screen? (Defaults to TRUE)
#' @param ... all parameters that are passed to the \code{\link{descent_generalized_fista}} function.
#' E.g. maxiter, tweak_vec etc ...
#' @param warm.start logical, should the solution of a previous model of the cross validation be used as a start
#'  in the next model. Defaults to TRUE. Notice, that the warm.start starts with the most unpenalized model.
#' @param estimate.c.type string, either "non_negative", or "direct".
#' Indicates how the algorithm finds the solution of \eqn{arg min_C ||diag(g)(Y - XC)||_2}.
#' If 'estimate.c.type' is set to "direct" there is no regularization (see \code{\link{estimate_c}}),
#' if 'estimate.c.type' is set to "non_negative" the estimates "C" must not be negative (non-negative least squares)
#' (see (see \code{\link{estimate_nn_c}}))
#'
#' @return list of length 2.
#' \itemize{
#'    \item 'cv.obj', list of lists. DTD model for each lambda, and every folds.
#'    \item 'best.model', list. DTD model optimized on the complete data set with the best lambda from the cross validation.
#' }
#' A cross validation matrix as entry "cv.obj", and the model with minimal loss function
#' retrained on the complete dataset as "best.model
#' @export
#'
DTD_cv_lambda_cxx <- function(
  lambda.seq = NULL,
  tweak.start,
  X.matrix,
  n.folds = 5,
  lambda.length = 10,
  train.data.list,
  cv.verbose = TRUE,
  warm.start = FALSE,
  estimate.c.type = "direct",
  NORM.FUN = "identity",
  NESTEROV.FUN = "positive",
  ST.FUN = "softmax",
  ...
  ) {
  DTD_cv_lambda_test_input_generic(lambda.seq, tweak.start, n.folds, lambda.length, train.data.list, cv.verbose, warm.start)
  if( ! is.numeric(X.matrix) ||
      ! is.matrix(X.matrix) ||
      nrow(X.matrix) != nrow(train.data.list$mixtures) ||
      ncol(X.matrix) != nrow(train.data.list$quantities) ) {
    cat("X is a ", nrow(X.matrix), " times ", ncol(X.matrix), " matrix.")
    cat("train.data.list$mixtures: ", nrow(train.data.list$mixtures), " x ", ncol(train.data.list$mixtures))
    cat("train.data.list$quantities: ", nrow(train.data.list$quantities), " x ", ncol(train.data.list$quantities))
    stop("X has wrong type or size (should be a numeric ngenes x ncells matrix.")
  }
  if( ! is.character(NORM.FUN) ) {
    stop("NORM.FUN must be a character string when using the cxx implementation (use the R impl to use user defined functions).")
  }
  if( ! is.character(NESTEROV.FUN) ) {
    stop("NESTEROV.FUN must be a character string when using the cxx implementation (use the R impl to use user defined functions).")
  }
  if( ! is.character(ST.FUN) ) {
    stop("ST.FUN must be a character string when using the cxx implementation (use the R impl to use user defined functions).")
  }
  # short hand:
  train.Y <- train.data.list$mixtures
  bucket.indicator <- make_buckets(
    train.Y = train.Y,
    folds = n.folds)

  if (!is.numeric(lambda.length)) {
    lambda.length <- 20
  }
  lambda.seq <- lambda_sequence(lambda.seq, lambda.length, train.Y)

  # after cross validation, a model is trained on the complete training
  # data using the best lambda of the cross validation.
  # The starting tweak for the end model must not be the result of a warm start,
  # therefore store it:
  tweak.start.end.model <- tweak.start

  # in older versions, I only kept the resulting loss for every fold and each lambda.
  # I think, for visualizing the cross validation, and comparing, it is better to return the complete models
  cv.object <- list()

  # prepare the model:
  model <- empty_model()
  model$X <- X.matrix
  model <- set_model_coeff_estimation(model, estimate.c.type)
  model <- set_model_normfunction(model, NORM.FUN)
  model <- set_model_subspacefunction(model, NESTEROV.FUN)
  model <- set_model_threshfunction(model, ST.FUN)
  # Start of cross validation:
  for (lambda in lambda.seq) {
    if (cv.verbose) {
      pos <- which(lambda.seq == lambda) - 1
      cat("\ndoing lambda: ", lambda, ", completed ", pos, " of ", length(lambda.seq), ", ", 100 * pos / length(lambda.seq), "% \n")
    }
    if (cv.verbose) {
      cat("doing l.fold: ")
    }
    lambda.fold <- list()
    for (l.fold in 1:n.folds) {
      if (cv.verbose) {
        cat(l.fold, "\t")
      }

      # Split the complete training data into test and train:
      test.samples <- names(which(bucket.indicator == l.fold))
      train.samples <- names(which(bucket.indicator != l.fold))

      # reduce the train to only include the cv train samples ...
      tmp.train.list <- lapply(train.data.list, select.fun, samples = train.samples)

      # Now, try to train a model on the reduced training set:
      model$Y <- tmp.train.list$mixtures
      model$C <- tmp.train.list$quantities
      model$tweak <- tweak.start

      ## This is a useful debug option
      ## options(try.outFile=stdout())
      catch <- try(descent_generalized_fista_cxx(
                              model = model,
                              lambda = lambda,
                              ...
                            ),
      silent = TRUE
      )


      # If the regularization parameter lambda is too big,
      # the fista algorithm does not find a model, and throws an error
      if (any(grepl(pattern = "Error", catch))) {
        lambda.fold[[as.character(l.fold)]] <- paste("could not build a model", catch)
        next
      }
      # warm start, after learning a model, keep last tweak vec as start for next model:
      if (warm.start) {
        tweak.start <- catch$Tweak
      }
      # Evaluate the reached minimum on the test set:
      tmp.test.list <- lapply(train.data.list, select.fun, samples = test.samples)
      tmp.eval.fun.test <- function(tmp.tweak, tmp.list = tmp.test.list) {
        thismodel <- empty_model()
        thismodel$Y <- tmp.list$mixtures
        thismodel$C <- tmp.list$quantities
        thismodel$X <- X.matrix
        thismodel$tweak <- tmp.tweak
        thisloss <- try(cxx_evaluate_model(thismodel), silent=TRUE)
        if (any(grepl(pattern = "Error", thisloss))) {
          print(thisloss)
          return(NA)
        }
        return(thisloss)
      }
      catch$cor.test <- tmp.eval.fun.test(catch$Tweak)
      lambda.fold[[as.character(l.fold)]] <- catch
    }
    cv.object[[as.character(lambda)]] <- lambda.fold
  }


  # pick the average mean per lambda:
  test.result.per.lambda <- unlist(
    lapply(
      cv.object
      , pick.mean.test.results.function
      )
  )

  if(all(is.na(test.result.per.lambda))){
    stop(
      "in 'DTD_cv_lambda_cxx': all lambdas crashed. \n
      Did you provide a custom lambda sequence? \n
      If yes: provide smaller lambdas \n
      If no: provide more features to the model,
      or provide a custom lambda sequence
      (via 'lambda.seq' in train_deconvolution_model(...). \n
      in 'train_deconvolution_model' set 'cv.verbose=TRUE' for logging information"
    )
  }

  mean.test.results <- test.result.per.lambda

  # and rebuild a model on the complete dataset:
  lmin.pos <- which.min(mean.test.results)
  lmin <- as.numeric(names(mean.test.results)[lmin.pos])

  if (cv.verbose) {
    cat("\ncross validation completed, \n starting to build model on complete data, \n with  best lambda = ", lmin, "\n")
  }

  model$tweak <- tweak.start.end.model
  model$Y <- train.data.list$mixtures
  model$C <- train.data.list$quantities

  bestModel <- descent_generalized_fista_cxx(
    model,
    lambda = lmin,
    save.all.tweaks = TRUE,
    ...
  )

  # return the cv.object for plotting, and the model with best lambda
  ret <- list(cv.obj = cv.object, best.model = bestModel)
  return(ret)
}
