#' input tests for cross validation
#'
#' tests common input parameters to the cxx and R implementation of DTD_cv_lambda.
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
#' @param cv.verbose logical, should information about the cv process be
#' printed to the screen?
#' @param warm.start logical, should the solution of a previous model of
#' the cross validation be used as start in the next model.
#' Notice, that the warm.start starts with the most unpenalized model.
#'
#'  @return NULL, or it throws an error
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
                       min = 2,
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
      stop("In 'DTD_cv_lambda': entries of train.data.list must be named 'quantities' and 'mixtures'")
    }else{
      if(!is.matrix(train.data.list$mixtures)){
        stop("In 'DTD_cv_lambda': 'train.data.list$mixtures' is not a matrix")
      }
      if(!is.matrix(train.data.list$quantities)){
        stop("In 'DTD_cv_lambda': 'train.data.list$quantities' is not a matrix")
      }
    }
  }else{
    stop("In 'DTD_cv_lambda': train.data.list must be provided as a list with two entries: 'quantities' and 'mixtures'")
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
  return(NULL)
}
