#' Estimating non negative C
#'
#' Given a reference matrix X, a matrix of bulks Y and a g-vector, "estimate_nn_c" finds
#' the solution of \deqn{arg min || diag(g) (Y - XC) ||_2} with non negative constraint over C,
#' using the package 'nnls' (non-negative least squares)
#'
#' @param X.matrix numeric matrix with cells as columns, and features as rows.
#'  Reference matrix X of the DTD problem. X.matrix can be set to NA (default), if the DTD.model
#'  includes the reference matrix X (default for \code{\link{train_correlatio_model}})
#' @param new.data numeric matrix with samples as columns, and features as rows.
#' In the formula above denoated as Y.
#' @param DTD.model either a numeric vector with length of nrow(X),
#' or a list returned by \code{\link{train_correlatio_model}}, \code{\link{DTD_cv_lambda}},
#' or\code{\link{descent_generalized_fista}}. In the equation above
#'   the DTD.model provides the vector g.
#'
#' @return numeric matrix with ncol(X) rows, and ncol(Y) columns
#'
#' @import nnls
#' @export
#' @examples
#' library(DTD)
#'
#' set.seed(1)
#' # simulate random data:
#' random.data <- generate_random_data(
#'   n.types = 5,
#'   n.samples.per.type = 1,
#'   n.features = 100
#' )
#'
#' # simulate a true c
#' # (this is not used by the estimate_c function, it is only used to show the result!)
#' true.c <- rnorm(n = ncol(random.data), mean = 0.1, sd = 0.5)
#'
#' # calculate bulk y = Xc * some_error
#' bulk <- as.matrix(random.data %*% true.c * rnorm(n = nrow(random.data), mean = 1, sd = 0.01), ncol = 1)
#' colnames(bulk) <- "mixture1"
#'
#' # estimate c
#' estimated.c <- estimate_nn_c(
#'   X.matrix = random.data,
#'   new.data = bulk,
#'   DTD.model = rep(1, nrow(random.data))
#' )
#'
#' # visualize that the estimated c are close to the true c
#' plot(true.c, estimated.c)
estimate_nn_c <- function(X.matrix = NA, new.data, DTD.model){
  if (any(is.na(X.matrix)) && is.list(DTD.model) && "reference.X" %in%
      names(DTD.model)) {
    X <- DTD.model$reference.X
  }
  else {
    X <- X.matrix
  }
  if (is.list(DTD.model)) {
    if ("best.model" %in% names(DTD.model)) {
      gamma.vec <- DTD.model$best.model$Tweak
    }
    else {
      if (!("Tweak" %in% names(DTD.model))) {
        stop("estimate_c: DTD.model does not fit")
      }
      else {
        gamma.vec <- DTD.model$Tweak
      }
    }
  }
  else {
    gamma.vec <- DTD.model
  }
  estimates <- matrix(NA, nrow = ncol(X), ncol = ncol(new.data))
  rownames(estimates) <- colnames(X)
  colnames(estimates) <- colnames(new.data)
  g.X.matrix <- diag(gamma.vec) %*% X
  for(l.mix in colnames(new.data)){
    g.new.data <- diag(gamma.vec) %*% new.data[, l.mix]
    estimates.l.mix <- nnls::nnls(g.X.matrix, g.new.data)
    estimates[, l.mix] <- estimates.l.mix$x
  }
  return(estimates)
}

