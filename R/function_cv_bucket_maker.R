#' map every sample to a bucket
#'
#' takes training data and sorts it into buckets, one for each fold.
#'
#' @param train.Y numeric matrix, training data
#' @param folds integer, number of buckets to build.
#'
#' @return named vector
#'
#' @examples
#' some.mat <- matrix(1:20, nrow = 2)
#' colnames(some.mat) <- paste0("some.mat", 1:ncol(some.mat))
#' n.folds <- 5
#' mapped.buckets <- make_buckets(
#'   train.Y = some.mat
#'   , folds = n.folds
#'   )
#' table(mapped.buckets)
make_buckets <- function(train.Y, folds) {
  bucket.indicator <- sample(
    rep(1:folds
        , each = ceiling(ncol(train.Y) / folds)
    )
  )[1:ncol(train.Y)]
  names(bucket.indicator) <- colnames(train.Y)
  return(bucket.indicator)
}
