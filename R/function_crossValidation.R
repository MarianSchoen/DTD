#' Cross-validation for digital tissue deconvolution
#'
#' Our descent generalized FISTA implementation includes a l1 regularization term.
#' This function performs a k-fold cross validation to find the best fitting regularization parameter.
#' For an example see `browseVignettes("DTD")`
#'
#' @param lambda.seq numeric vector or NULL: Over this series of lambdas the FISTA will be optimized.
#' If lambda.seq is set to NULL, a generic series of lambdas - depending on the dimensions
#' of the training set -  will be generated. Default: NULL
#' @param tweak.start numeric vector, starting vector for the DTD algorithm. Default: NULL.
#' @param n.folds integer, number of buckets in the cross validation. Defaults to 10
#' @param lambda.length integer, how many lambdas will be generated (only used if lambda.seq is NULL). Defaults to 20
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
#'  in the next model. Defaults to FALSE
#'
#' @return list of length 2. A cross validation matrix as entry "cv.obj", and the model with minimal loss function
#' retrained on the complete dataset as "best.model
#' @export
#'
DTD_cv_lambda <- function(lambda.seq = NULL,
                          tweak.start = NULL,
                          n.folds = 5,
                          lambda.length = 10,
                          train.data.list,
                          F.GRAD.FUN,
                          EVAL.FUN,
                          cv.verbose = TRUE,
                          warm.start = FALSE,
                          ...) {
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

  # First, all possible training samples get assigned to a bucket:
  # extract Y:
  train.Y <- train.data.list$mixtures
  # map every sample to a bucket:
  bucket.indicator <- sample(rep(1:n.folds,
    each = ceiling(ncol(train.Y) / n.folds)
  ))[1:ncol(train.Y)]
  names(bucket.indicator) <- colnames(train.Y)

  # Next, if no lambda.seq is provided, a generic sequence is created,
  # based on "p" and "n" of the mixture matrix:
  if (is.null(lambda.seq) || is.na(lambda.seq) || !is.numeric(lambda.seq)) {
    if (!is.numeric(lambda.length)) {
      lambda.length <- 20
    }
    p <- nrow(train.Y)
    n <- ncol(train.Y)
    lambda.0 <- sqrt(log(p) / n)
    lambda.seq <- lambda.0 * 2^seq(2, -20, length.out = lambda.length)
  }

  # internal training samples selection function:
  select.fun <- function(list.entry, samples) {
    if (is.matrix(list.entry)) {
      return(list.entry[, samples])
    }
    if (!is.null(names(list.entry))) {
      return(list.entry[samples])
    }
    return(list.entry)
  }
  # Initialise the cv.object:
  cv.object <- data.frame(matrix(
    nrow = (n.folds + 2),
    ncol = length(lambda.seq)
  ))
  colnames(cv.object) <- lambda.seq
  rownames(cv.object) <- c(1:n.folds, "nonZero", "nFoundModels")
  # Start of cross validation:
  for (lambda in lambda.seq) {
    cor.test.vec <- c()
    non_zeros <- c()
    foundMods <- 0
    if (cv.verbose) {
      pos <- which(lambda.seq == lambda) - 1
      cat("\ndoing lambda: ", lambda, ", completed ", pos, " of ", length(lambda.seq), ", ", 100 * pos / length(lambda.seq), "% \n")
    }
    if (cv.verbose) {
      cat("doing l.fold: ")
    }
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

      # If the regularization parameter lambda is to big,
      # the fista algorithm does not find a model, and throws an error
      if (any(grepl(pattern = "Error", catch))) {
        cor.test.vec <- c(cor.test.vec, NA)
        non_zeros <- c(non_zeros, NA)
        next
      } else {
        # if the fista algorithm did not throw an error => increment the number of found Models ...
        foundMods <- foundMods + 1
        # ... and check how many coefficients are unequal zero
        tmp_non_zero <- sum(catch$Tweak != 0)
        non_zeros <- c(non_zeros, tmp_non_zero)
        # warm start, after learning a model, keep last tweak vec as start for next model:
        if (warm.start) {
          tweak.start <- catch$Tweak
        }
      }
      # Evaluate the reached minimum on the test set:
      tmp.test.list <- lapply(train.data.list, select.fun, samples = test.samples)
      tmp.eval.fun.test <- function(tmp.tweak, tmp.list = tmp.test.list) {
        return(EVAL.FUN(tmp.tweak, train.list = tmp.list))
      }
      cor.test.vec <- c(cor.test.vec, tmp.eval.fun.test(catch$Tweak))
    }
    # fill the cv.object data.frame:
    cv.object[[as.character(lambda)]] <- c(cor.test.vec, mean(non_zeros), foundMods / n.folds)
  }
  if (cv.verbose) {
    cat("\ncross validation completed, starting to build model on complete data, with  best lambda\n")
  }

  # after the cross validation, find the lambda with best evaluation score,
  # and rebuild a model on the complete dataset:
  lmin.pos <- which.min(apply(cv.object[1:n.folds, , drop = FALSE],
                              2,
                              mean, na.rm = TRUE)) # it may occur (due to hapless folds) that one fold leads to NA in all models
  lmin <- as.numeric(colnames(cv.object)[lmin.pos])

  bestModel <- descent_generalized_fista(
    lambda = lmin,
    tweak.vec = tweak.start,
    F.GRAD.FUN = F.GRAD.FUN,
    EVAL.FUN = EVAL.FUN,
    save.all.tweaks = TRUE,
    ...
  )

  # return the cv.object for plotting, and the model with best lambda
  ret <- list(cv.obj = cv.object, best.model = bestModel)
  return(ret)
}
