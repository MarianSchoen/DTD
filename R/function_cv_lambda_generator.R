#' generates a sequence of lambdas
#'
#' depending on the input, generates a sequence of lambdas or makes it usable (sorts it)
#'
#' @param lambda.seq sequence of lambdas. may be null.
#' @param lambda.length number of lambdas to generate
#' @param train.Y training data, used for heuristics involving the dimensionality of the problem
#'
#' @return sequence of lambdas (numeric vector)
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
    lambda.seq <- sort(lambda.seq, decreasing = FALSE)
  }
  return(lambda.seq)
}
