#' Cross-validation for digital tissue deconvolution
#'
#' Our descent generalized FISTA implementation includes a l1 regularization term.
#' This function performs a k-fold cross validation to find the best fitting regularization parameter.
#'
#' @param lambda.seq numeric vector or NULL: Over this series of lambdas the FISTA will be optimized.
#' If lambda.seq is set to NULL, a generic series of lambdas - depending on the dimensions
#' of the training set -  will be generated
#' @param nfolds integer, number of buckets in the cross validation
#' @param lambda.length integer, how many lambdas will be generated (only used if lambda.seq is NULL)
#' @param train.list list, that can be passed to the GRAD.FUN and EVAL.FUN.
#' Within this list the train/test cross validation will be done.
#' Notice, that the train.list must have an entry named "mixtues". In this entry, the matrix containing the
#' training samples (as columns) and all features (as rows) must be present. (see Vignette for details)
#' @param GRAD.FUN gradient function, see \code{\link{descent_generalized_fista}}
#' @param EVAL.FUN evaluation function, see \code{\link{descent_generalized_fista}}
#' @param tol float, in each cross validation model, the function keeps track how many explaining
#' variables do not contribute (=> equal 0). tol is the limit until a explaining variable is
#' declared as "not contributing".
#' @param ... all parameters that are passed to the \code{\link{descent_generalized_fista}} function.
#' E.g. maxiter, tweak_vec etc ...
#'
#' @return list of length 2. A cross validation matrix, and the model with minimal loss function
#' retrained on the complete dataset
#' @export
#'
#' @examples
#'


DTD_cv_lambda <- function(lambda.seq = NULL,
                          nfolds = 10,
                          lambda.length = 20,
                          train.list = train,
                          GRAD.FUN = DTD.grad.wrapper,
                          EVAL.FUN = DTD.evCor.wrapper,
                          tol = 1e-5,
                          cv.verbose = TRUE,
                          ...
){

  # First, all possible training samples get assigned to a bucket:
  # extract Y:
  train.Y <- train.list$mixtures
  # map every sample to a bucket:
  bucket.indicator <- sample(rep(1:nfolds, each = ceiling(ncol(train.Y)/nfolds)))[1:ncol(train.Y)]
  names(bucket.indicator) <- colnames(train.Y)


  # Next, if no lambda.seq is provided, a generic sequence is created,
  # based on "p" and "n" of the mixture matrix:
  if(is.null(lambda.seq)){
    if(!is.numeric(lambda.length)){lambda.length <- 20}
    p <- nrow(train.Y)
    n <- ncol(train.Y)
    lambda.0 <- 0.5*sqrt(log(p)/n)
    lambda.seq <- lambda.0*2^seq(2, -5, length.out = lambda.length)
  }

  # internal training samples selection function:
  select.fun <- function(list.entry, train.samples){
    if(is.matrix(list.entry)){
      return(list.entry[, train.samples])
    }
    if(!is.null(names(list.entry))){
      return(list.entry[train.samples])
    }
    return(list.entry)
  }

  # Initialise the cv.object:
  cv.object <- data.frame("lambda" = lambda.seq,
                          "cvm" = rep(NA, length(lambda.seq)),
                          "cvsd" = rep(NA, length(lambda.seq)),
                          "cvup" = rep(NA, length(lambda.seq)),
                          "cvlo" = rep(NA, length(lambda.seq)),
                          "nzero" = rep(NA, length(lambda.seq)),
                          "nfoundModels" = rep(NA, length(lambda.seq)))

  # Start of cross validation:
  for(lambda in lambda.seq){
    if(cv.verbose){
      cat("doing lambda: ", lambda, "\n")
      }
    cor.test.vec <- c()
    non_zeros <- c()
    foundMods <- 0
    if(cv.verbose){cat("doing l.fold: ")}
    for(l.fold in 1:nfolds){
      if(cv.verbose){cat(l.fold, "\t")}
      # Split the complete training data into test and train:
      test.samples <- names(which(bucket.indicator == l.fold))
      train.samples <- names(which(bucket.indicator != l.fold))

      # reduce the train to only include the cv train samples ...
      tmp.train.list <- lapply(train.list, select.fun)

      # ... and reset the default values of the gradient and evaluation functions:
      tmp.grad.fun <- function(tmp.tweak, tmp.train = tmp.train.list){
        return(GRAD.FUN(tmp.tweak, train.list = tmp.train))
      }
      tmp.eval.fun <- function(tmp.tweak, tmp.train = tmp.train.list){
        return(EVAL.FUN(tmp.tweak, train.list = tmp.train))
      }

      # Now, try to train a model on the reduced training set:
      catch <- try(descent_generalized_fista(lambda = lambda,
                                             F.GRAD.FUN = tmp.grad.fun,
                                             EVAL.FUN = tmp.eval.fun,
                                             ...)
      )

      # If the regularization parameter lambda is to big,
      # the fista algorithm does not find a model, and throws an error
      if(any(grepl(pattern = "Error", catch))){
        next
      }else{
        # if the fista algorithm did not throw an error => increment the number of found Models ...
        foundMods <- foundMods + 1
        # ... and check how many samples are close to zero
        tmp_non_zero <- sum(
          unlist(
            lapply(catch$Tweak,
                   function(x){isTRUE(all.equal(x, 0, tolerance = tol))})
          )
        )
        non_zeros <- c(non_zeros, tmp_non_zero)
      }
      # Evaluate the reached minimum on the test set:
      cor.test.vec <- c(cor.test.vec, rep(x = EVAL.FUN(catch$Tweak), length(test.samples)))
    }
    if(cv.verbose){cat("\n")}
    # fill the cv.object data.frame:
    pos <- which(cv.object[["lambda"]] == lambda)
    cv.object[["cvm"]][pos] <- suppressWarnings(mean(cor.test.vec, na.rm = TRUE))
    cv.object[["cvsd"]][pos] <- suppressWarnings(sd(cor.test.vec, na.rm = TRUE))
    cv.object[["cvup"]][pos] <- suppressWarnings(max(cor.test.vec, na.rm = TRUE))
    cv.object[["cvlo"]][pos] <- suppressWarnings(min(cor.test.vec, na.rm = TRUE))
    cv.object[["nzero"]][pos] <- suppressWarnings(mean(non_zeros, na.rm = TRUE))
    cv.object[["nfoundModels"]][[pos]] <- foundMods
  }

  # after the cross validation, find the lambda with best evaluation score,
  # and rebuild a model on the complete dataset:

  lmin.pos <- which.min(cv.object$cvm)
  lmin <- cv.object$lambda[lmin.pos]
  bestModel <- descent_generalized_fista(lambda = lmin,
                                         F.GRAD.FUN = GRAD.FUN,
                                         EVAL.FUN = EVAL.FUN,
                                         ...)

  # return the cv.object for plotting, and the model with best lambda
  ret <- list(cv.obj = cv.object, best.model = bestModel)
  return(ret)
}
