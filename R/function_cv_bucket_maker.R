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
  bucket.indicator <- sample(
    rep(1:folds
        , each = ceiling(ncol(train.Y) / folds)
    )
  )[1:ncol(train.Y)]
  names(bucket.indicator) <- colnames(train.Y)
  return(bucket.indicator)
}
