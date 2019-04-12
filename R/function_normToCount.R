#' normalize_to_count
#'
#' 'normalize_to_count' takes a numeric matrix 'exp.data' with only non-negative
#' entries and normalizes each column (=> sample) to a total number of counts.\cr
#' For every sample i, each feature j gets scaled to
#'  \deqn{exp.data[j, i] = (count * exp.data[j, i]) / sum(exp.data[, i])}
#'
#' @param exp.data positive numeric matrix, with samples as columns and features as rows.
#' Notice, 'normalize_to_count' normalizes the columns of this matrix.
#' @param count float, 0 < 'count', to which every sample of exp.data will be scaled.
#' If 'is.na(count)' (default), count is set to nrow(exp.data)
#'
#' @return "ret", matrix with same dimension as "exp.data"
#' @export
#'
#' @examples
#' library(DTD)
#' some.matrix <- matrix(abs(rnorm(1000 * 5)), ncol = 5, nrow = 1000)
#' # each sample (=column) has different number of total counts:
#' apply(some.matrix, 2, sum)
#'
#' normalized.matrix <- normalize_to_count(some.matrix)
#'
#' # check:
#' apply(normalized.matrix, 2, sum)
normalize_to_count <- function(exp.data,
                               count = NA) {
  if(any(is.na(count))){
    count <- nrow(exp.data)
  }
  # safety check: count
  test <- test_numeric(test.value = count,
                       output.info = c("normalize_to_count", "count"),
                       min = .Machine$double.eps,
                       max = Inf)
  # end -> count

  # saftey check: exp.data
  if(any(is.na(exp.data))){
    stop("In normalize_to_count: exp.data includes NA")
  }
  if(any(exp.data < 0)){
    stop("In normalize_to_count: exp.data includes negative values")
  }
  if(!is.matrix(exp.data)){
    stop("In normalize_to_count: exp.data is not of matrix format")
  }
  # end -> exp.data

  # microbenchmark told me that apply is faster than a loop,
  # therefore normalizing is done via apply:
  norm.data <- apply(exp.data, 2, function(x){count*x/sum(x)})

  return(norm.data)
}
