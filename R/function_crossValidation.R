#' Cross-validation for digital tissue deconvolution
#'
#' Our descent generalized FISTA implementation includes a l1 regularization term (see \code{\link{train_deconvolution_model}}.
#' This function performs a k-fold cross validation to find the best fitting regularization parameter.
#' For an example see `browseVignettes("DTD")`
#'
#' @param lambda.seq numeric vector or NULL: Over this series of lambdas the FISTA will be optimized.
#' If lambda.seq is set to NULL, a generic series of lambdas - depending on the dimensions
#' of the training set -  will be generated. Default: NULL
#' @param tweak.start numeric vector, starting vector for the DTD algorithm.
#' @param n.folds integer, number of buckets in the cross validation. Defaults to 10
#' @param lambda.length integer, how many lambdas will be generated (only used if lambda.seq is NULL). Defaults to 10
#' @param train.data.list list, that can be passed to the F.GRAD.FUN and EVAL.FUN.
#' Within this list the train/test cross validation will be done.
#' Notice, that the train.data.list must have an entry named "mixtues". In this entry, the matrix containing the
#' training samples (as columns) and all features (as rows) must be present. (see Vignette for details)
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
DTD_cv_lambda_R <- function(lambda.seq = NULL,
                            tweak.start,
                            n.folds = 5,
                            lambda.length = 10,
                            train.data.list,
                            F.GRAD.FUN,
                            EVAL.FUN,
                            cv.verbose = TRUE,
                            warm.start = FALSE,
                            ...) {

  DTD_cv_lambda_test_input_generic(lambda.seq, tweak.start, n.folds, lambda.length, train.data.list, cv.verbose, warm.start)

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
    cat("\ncross validation completed, starting to build model on complete data, with  best lambda\n")
  }

  # after the cross validation, find the lambda with best evaluation score
  # pick the average mean per lambda:
  test.result.per.lambda <- lapply(
    cv.object,
    pick.mean.test.results.function
  )

  # unlist it => keep names
  mean.test.results <- unlist(test.result.per.lambda)

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

#' input tests for cross validation
#'
#' tests common input parameters to the cxx and R implementation of DTD_cv_lambda.
#'
#'
#' @param lambda.seq numeric vector or NULL: Over this series of lambdas the FISTA will be optimized.
#' If lambda.seq is set to NULL, a generic series of lambdas - depending on the dimensions
#' of the training set -  will be generated. Default: NULL
#' @param tweak.start numeric vector, starting vector for the DTD algorithm.
#' @param n.folds integer, number of buckets in the cross validation. Defaults to 10
#' @param lambda.length integer, how many lambdas will be generated (only used if lambda.seq is NULL). Defaults to 10
#' @param train.data.list list, that can be passed to the F.GRAD.FUN and EVAL.FUN.
#' Within this list the train/test cross validation will be done.
#' Notice, that the train.data.list must have an entry named "mixtues". In this entry, the matrix containing the
#' training samples (as columns) and all features (as rows) must be present. (see Vignette for details)
#' @param cv.verbose logical, should information about the cv process be printed to the screen? (Defaults to TRUE)
#' @param warm.start logical, should the solution of a previous model of the cross validation be used as a start
#'  in the next model. Defaults to TRUE. Notice, that the warm.start starts with the most unpenalized model.
#'
#' @export
DTD_cv_lambda_test_input_generic <- function(lambda.seq,
                                 tweak.start,
                                 n.folds,
                                 lambda.length,
                                 train.data.list,
                                 cv.verbose,
                                 warm.start) {

  # safety check: tweak.start
  test <- test_tweak_vec(tweak.vec = tweak.start,
                         output.info = c("DTD_cv_lambda", "tweak.start"))
  # end -> tweak

  # safety check: n.folds
  test <- test_integer(test.value = n.folds,
                       output.info = c("DTD_cv_lambda", "n.folds"),
                       min = 1,
                       max = Inf
                       )
  # end -> n.folds

  # safety check: lambda.length
  test <- test_integer(test.value = lambda.length,
                       output.info = c("DTD_cv_lambda", "lambda.length"),
                       min = 1,
                       max = Inf
  )
  # end -> lambda.length

  # safety check for train.data.list:
  if(is.list(train.data.list) && length(train.data.list) == 2){
    if(!all(c("quantities", "mixtures") %in%  names(train.data.list))){
      stop("In DTD_cv_lambda: entries of train.data.list must be named 'quantities' and 'mixtures'")
    }else{
      if(!is.matrix(train.data.list$mixtures)){
        stop("In DTD_cv_lambda: 'train.data.list$mixtures' is not a matrix")
      }
      if(!is.matrix(train.data.list$quantities)){
        stop("In DTD_cv_lambda: 'train.data.list$quantities' is not a matrix")
      }
    }
  }else{
    stop("In DTD_cv_lambda: train.data.list must be provided as a list with two entries: 'quantities' and 'mixtures'")
  }
  # end -> train.data.list

  # safety check => tweak.start and train.data.list
  if(nrow(train.data.list$mixtures) != length(tweak.start)){
    stop("In DTD_cv_lambda: 'nrow(train.data.list$mixtures)' does not match 'length(tweak.start)'")
  }
  if(!is.null(names(tweak.start))){
    if(all(rownames(train.data.list$mixtures) %in% names(tweak.start))){
      tweak.start <- tweak.start[rownames(train.data.list$mixtures)]
    }else{
      stop("In DTD_cv_lambda: there are features in 'train.data.list$mixtures' that don't match with 'names(tweak.start)'")
    }
  }
  # end -> tweak and train.data.list

  # safety check: cv.verbose
  test <- test_logical(test.value = cv.verbose,
                       output.info = c("DTD_cv_lambda", "cv.verbose"))
  # end -> cv.verbose

  # safety check: warm.start
  test <- test_logical(test.value = warm.start,
                       output.info = c("DTD_cv_lambda", "warm.start"))
  # end -> warm.start
  ####################### end safety check
}

#' map every sample to a bucket
#'
#' takes training data and sorts it into buckets, one for each fold.
#'
#' @param train.Y training data
#' @param folds number of buckets to build.
#'
#' @return list of lists of same structure as train.Y, but with names indicating the bucket.
#'
make_buckets <- function(train.Y, folds) {
  bucket.indicator <- sample(rep(1:folds,
                             each = ceiling(ncol(train.Y) / folds)
                             ))[1:ncol(train.Y)]
  names(bucket.indicator) <- colnames(train.Y)
  return(bucket.indicator)
}

#' generates a sequence of lambdas
#'
#' depending on the input, generates a sequence of lambdas or makes it usable (sorts it)
#'
#' @param lambda.seq sequence of lambdas. may be null.
#' @param lambda.length number of lambdas to generate
#' @param train.Y training data, used for heuristics involving the dimensionality of the problem
#'
#' @return sequence of lambdas
#'
lambda_sequence <- function(lambda.seq, lambda.length, train.Y) {
  # if necessary, generate a generic sequence,
  # based on "p" and "n" of the mixture matrix:
  if (is.null(lambda.seq) || is.na(lambda.seq) || !is.numeric(lambda.seq)) {
    p <- nrow(train.Y)
    n <- ncol(train.Y)
    lambda.0 <- sqrt(log(p) / n)
    lambda.seq <- lambda.0 * 2^seq(2, -20, length.out = lambda.length)
  }
  # cross validation can be called with warm.start.
  # these warm starts should start with the most unregularized scenario:
  if(min(lambda.seq) != lambda.seq[1]){
    lambda.seq <- sort(lambda.seq, decreasing = TRUE)
  }
  return(lambda.seq)
}
#' selects a specific sample from a list of entries
#'
#' @param list.entry list of entries (matrix or list)
#' @param samples samples, serving as an index to list.entry.
#'
#' @return specific entry of list.entry
#'
select.fun <- function(list.entry, samples) {
  # internal training samples selection function:
  if (is.matrix(list.entry)) {
    return(list.entry[, samples])
  }
  if (!is.null(names(list.entry))) {
    return(list.entry[samples])
  }
  return(list.entry)
}

#' returns mean of the individual test results, setting unconverged or unavailable results to Inf.
#'
#' @param lambda.list
#'
#' @return list of means.
#'
pick.mean.test.results.function <- function(lambda.list){
  tmp <- lapply(lambda.list, function(each.fold){
    if("cor.test" %in% names(each.fold)){
      return(each.fold$cor.test)
    }else{
      return(Inf)
    }
  })
  test.vec <- mean(unlist(tmp, use.names = FALSE), na.rm = TRUE)
  return(test.vec)
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
DTD_cv_lambda_cxx <- function(lambda.seq = NULL,
                              tweak.start,
                              X.matrix,
                              n.folds = 5,
                              lambda.length = 10,
                              train.data.list,
                              cv.verbose = TRUE,
                              warm.start = FALSE,
                              ...) {
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
  model <- list()
  model$X <- X.matrix
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
      model$tweak = tweak.start

      catch <- try(descent_generalized_fista_cxx(
                              model = model,
                              lambda = lambda,
                              ...
                            ),
      silent = TRUE
      )


      # TODO: won't happen. think about error handling in cpp and propagating them to R!
      # If the regularization parameter lambda is to big,
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
        thismodel <- list()
        thismodel$Y <- tmp.list$mixtures
        thismodel$C <- tmp.list$quantities
        thismodel$X <- X.matrix
        thismodel$tweak <- tmp.tweak
        return(evaluate_model(thismodel))
      }
      catch$cor.test <- tmp.eval.fun.test(catch$Tweak)
      lambda.fold[[as.character(l.fold)]] <- catch
    }
    cv.object[[as.character(lambda)]] <- lambda.fold
  }
  if (cv.verbose) {
    cat("\ncross validation completed, starting to build model on complete data, with  best lambda\n")
  }


  # pick the average mean per lambda:
  test.result.per.lambda <- lapply(cv.object,
                                   pick.mean.test.results.function)

  # unlist it => keep names
  mean.test.results <- unlist(test.result.per.lambda)

  # and rebuild a model on the complete dataset:
  lmin.pos <- which.min(mean.test.results)
  lmin <- as.numeric(names(mean.test.results)[lmin.pos])

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
#' Cross-validation for digital tissue deconvolution
#'
#' Our descent generalized FISTA implementation includes a l1 regularization term (see \code{\link{train_deconvolution_model}}.
#' This function performs a k-fold cross validation to find the best fitting regularization parameter.
#' For an example see `browseVignettes("DTD")`
#' This function chooses at runtime either the cxx or the R implementation of cross validation, based on the useImplementation parameter
#'
#' @param useImplementation if set to "R" calls DTD_cv_lambda_R, if set to "cxx" calls DTD_cv_lambda_cxx, passing on all arguments. See the respective implementations of these functions.
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
DTD_cv_lambda <- function(useImplementation = "R",
                          ...) {
  if( useImplementation == "R" ) {
    return(DTD_cv_lambda_R(...))
  } else if ( useImplementation == "cxx" || useImplementation == "cpp" ) {
    return(DTD_cv_lambda_cxx(...))
  } else {
    stop(paste("unknown implementation: ", useImplementation))
  }
}

