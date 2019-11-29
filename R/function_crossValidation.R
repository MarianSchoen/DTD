#' Cross-validation for digital tissue deconvolution
#'
#' Our descent generalized FISTA implementation includes a l1 regularization
#' term (see \code{\link{train_deconvolution_model}}).
#' This function performs a 'n.folds'-fold cross validation to find the
#' best fitting regularization parameter.
#'
#' For an example see `browseVignettes("DTD")`.
#'
#' Notice, there is an R and a C++ implementation of our optimizer.
#' Hence, there are two cross validation implementations,
#' calling either the R or C++ implementation:
#' \code{\link{DTD_cv_lambda_R}} and \code{\link{DTD_cv_lambda_cxx}}.
#'
#' @param lambda.seq numeric vector or NULL or "none": Over this series of lambdas the
#' FISTA will be optimized. If 'lambda.seq' is set to NULL, a generic series of
#' lambdas - depending on the dimensions of the training set -
#' will be generated. If 'lambda.seq' is "none", no cross validation is done.
#' Only one model with lambda = 0 is trained on the complete data set.
#' @param tweak.start numeric vector, starting vector for the DTD algorithm.
#' @param n.folds integer, number of buckets in the cross validation.
#' @param lambda.length integer, how many lambdas will be generated
#' (only used if lambda.seq is NULL)
#' @param train.data.list list, with two entries, a numeric matrix each,
#' named 'mixtures' and 'quantities'
#' Within this list the train/test cross validation will be done.
#' (see Vignette `browseVignettes("DTD")` for details)
#' @param F.GRAD.FUN gradient function, see
#' \code{\link{descent_generalized_fista}}
#' The 'F.GRAD.FUN' and 'EVAL.FUN' parameters are only present in the
#' R cross validation. With these parameters an alternativ gradient,
#' and evaluation function can be provided. Both functions are called
#' using only the tweak vector as first argument.
#' @param EVAL.FUN evaluation function,
#' see \code{\link{descent_generalized_fista}}
#' @param cv.verbose logical, should information about the cv process be
#' printed to the screen?
#' @param ... all parameters that are passed to the
#' \code{\link{descent_generalized_fista}} function.
#' E.g. 'maxiter', 'NORM.FUN',  'cycles' etc ...
#' @param warm.start logical, should the solution of a previous model of
#' the cross validation be used as start in the next model.
#' Notice, that the warm.start starts with the most unpenalized model.
#'
#' @return list of length 2:
#' \itemize{
#'    \item 'cv.obj', list of lists. DTD model for each lambda, and every folds.
#'    \item 'best.model', list. DTD model optimized on the complete data set
#'     with the best lambda from the cross validation.
#' }
#'
#' @export
DTD_cv_lambda_R <- function(
  lambda.seq = "none",
  tweak.start,
  n.folds = 5,
  lambda.length = 10,
  train.data.list,
  cv.verbose = TRUE,
  warm.start = FALSE,
  F.GRAD.FUN,
  EVAL.FUN,
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

  if( lambda.seq[1] == "none" ){
    the.args <- list(...)
    if(!'save.all.tweaks' %in% names(the.args)){ # save.all.tweaks is not set
      bestModel <- descent_generalized_fista(
        lambda = 0,
        tweak.vec = tweak.start,
        F.GRAD.FUN = F.GRAD.FUN,
        EVAL.FUN = EVAL.FUN,
        save.all.tweaks = TRUE,
        ...
      )
    } else { # save.all.tweaks is not set
      bestModel <- descent_generalized_fista(
        lambda = 0,
        tweak.vec = tweak.start,
        F.GRAD.FUN = F.GRAD.FUN,
        EVAL.FUN = EVAL.FUN,
        ...
      )
    }

    return(list(best.model = bestModel))
  }

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
      "in 'DTD_cv_lambda_R': all lambdas crashed. \n
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
#' Our descent generalized FISTA implementation includes a l1 regularization
#' term (see \code{\link{train_deconvolution_model}}).
#' This function performs a 'n.folds'-fold cross validation to find the
#' best fitting regularization parameter.#'
#'
#' For an example see `browseVignettes("DTD")`
#'
#' Notice, there is an R and a C++ implementation of our optimizer.
#' Hence, there are two cross validation implementations,
#' calling either the R or C++ implementation:
#'
#' \code{\link{DTD_cv_lambda_R}} and \code{\link{DTD_cv_lambda_cxx}}.
#'
#' @param lambda.seq numeric vector or NULL or "none": Over this series of lambdas the
#' FISTA will be optimized. If 'lambda.seq' is set to NULL, a generic series of
#' lambdas - depending on the dimensions of the training set -
#' will be generated. If 'lambda.seq' is "none", no cross validation is done.
#' Only one model with lambda = 0 is trained on the complete data set.
#' @param tweak.start numeric vector, starting vector for the DTD algorithm.
#' @param X.matrix numeric matrix, with features/genes as rows,
#' and cell types as column. Each column of X.matrix is a reference
#' expression profile
#' @param n.folds integer, number of buckets in the cross validation.
#' @param lambda.length integer, how many lambdas will be generated
#' (only used if lambda.seq is NULL)
#' @param train.data.list list, with two entries, a numeric matrix each,
#' named 'mixtures' and 'quantities'
#' Within this list the train/test cross validation will be done.
#' (see Vignette `browseVignettes("DTD")` for details)
#' @param cv.verbose logical, should information about the cv process be
#' printed to the screen?
#' @param ... all parameters that are passed to the c++ optimization function.
#' @param warm.start logical, should the solution of a previous model of
#' the cross validation be used as start in the next model.
#' Notice, that the warm.start starts with the most unpenalized model.
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
#' @param NORM.FUN string, after each gradient descent and nesterov step, the
#' the resulting tweak/g-vector can be normed. There are three implemenations:
#' \itemize{
#'    \item 'identity': No normalization, every entry stays as it is.
#'    \item 'n2normed': the vector is scaled to \eqn{||g||_2 = 1}
#'    \item 'n1normed': the vector is scaled to \eqn{||g||_1 = 1}
#' }
#' @param NESTEROV.FUN string, sets the nesterov function. The current
#' implementation restricts the result to be positive
#'  (due to the optimization constraint \eqn{g_i \ge 0})
#' @param ST.FUN string, sets the soft thresholding function.
#' @param inv.precision numeric, for the least squares solution (X^T G X)^-1
#' must be inverted.
#'
#' @return list of length 2:
#' \itemize{
#'    \item 'cv.obj', list of lists. DTD model for each lambda, and every folds.
#'    \item 'best.model', list. DTD model optimized on the complete data set
#'     with the best lambda from the cross validation.
#' }
#'
#' @export
DTD_cv_lambda_cxx <- function(
  lambda.seq = "none",
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
  inv.precision = 1e-12,
  ...
  ) {
  DTD_cv_lambda_test_input_generic(lambda.seq, tweak.start, n.folds, lambda.length, train.data.list, cv.verbose, warm.start)

  if(!is.null(lambda.seq)){
    if( lambda.seq[1] == "none" ){
      # prepare the model:
      model <- empty_model()
      model$X <- X.matrix
      model <- set_model_coeff_estimation(model, estimate.c.type)
      model <- set_model_normfunction(model, NORM.FUN)
      model <- set_model_subspacefunction(model, NESTEROV.FUN)
      model <- set_model_threshfunction(model, ST.FUN)
      model <- set_model_inversion_precision(model, inv.precision)

      model$Y <- train.data.list$mixtures
      model$C <- train.data.list$quantities
      model$tweak <- tweak.start

      the.args <- list(...)
      if(!'save.all.tweaks' %in% names(the.args)){
        bestModel <- descent_generalized_fista_cxx(
          model = model,
          lambda = 0,
          save.all.tweaks = TRUE,
          ...
        )
      }else{
        bestModel <- descent_generalized_fista_cxx(
          model = model,
          lambda = 0,
          ...
        )
      }
      return(list(best.model = bestModel))
    }
  }


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

  if(any(is.na(inv.precision))){
    inv.prec.not.set <- TRUE
    inv.precision <- 1e-12
  }else{
    inv.prec.not.set <- FALSE
    test <- test_numeric(
      test.value = inv.precision
      , output.info = c("DTD_cv_lambda_cxx", "inv.precision")
      , min = 0
      , max = Inf
    )
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
  model <- set_model_inversion_precision(model, inv.precision)

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

  # to ensure that after the cross validation,
  # a slightly to small precision does not crashes to model
  if(inv.prec.not.set){
    model <- set_model_inversion_precision(1e-10)
  }

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
