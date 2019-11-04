#' Estimating C
#'
#' Given a reference matrix X, a matrix of bulks Y and a g-vector, "estimate_c" finds
#' the solution of \deqn{arg min || diag(g) (Y - XC) ||_2} over all C using
#' direct analytical solution: (see Goertler et al. 2018)
#' \deqn{ C(g) = (X^T \Gamma X )^(-1) X^T \Gamma Y}
#' with \eqn{\Gamma} = diag(g)
#'
#' @param X.matrix numeric matrix with cells as columns, and features as rows.
#'  Reference matrix X of the DTD problem. X.matrix can be set to NA (default), if the DTD.model
#'  includes the reference matrix X (default for \code{\link{train_deconvolution_model}})
#' @param new.data numeric matrix with samples as columns, and features as rows.
#' In the formula above denoated as Y.
#' @param DTD.model either a numeric vector with length of nrow(X),
#' or a list returned by \code{\link{train_deconvolution_model}}, \code{\link{DTD_cv_lambda}},
#' or\code{\link{descent_generalized_fista}}. In the equation above
#'   the DTD.model provides the vector g.
#'
#' @return numeric matrix with ncol(X) rows, and ncol(Y) columns
#'
#' @export
#' @examples
#' library(DTD)
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
#' bulk <- random.data %*% true.c * rnorm(n = nrow(random.data), mean = 1, sd = 0.01)
#'
#' # estimate c
#' estimated.c <- estimate_c(
#'   X.matrix = random.data,
#'   new.data = bulk,
#'   DTD.model = rep(1, nrow(random.data)),
#'   estimate.c.type = "direct"
#' )
#'
#' # visualize that the estimated c are close to the true c
#' plot(true.c, estimated.c)
estimate_c <- function(
  X.matrix = NA,
  new.data,
  DTD.model,
  estimate.c.type = "decide.on.model"
  ) {
  # the reference matrix can either be included in the DTD.model, or has to be past
  # via the X.matrix argument:
  if (any(is.na(X.matrix)) && is.list(DTD.model) && "reference.X" %in% names(DTD.model)) {
    X <- DTD.model$reference.X
  } else {
    X <- X.matrix
  }
  # as a DTD.model either a list, or only the tweak vector can be used:
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
    if ("estimate.c.type" %in% names(DTD.model)){
        model.esti.type <- DTD.model$estimate.c.type
    }
  }
  else {
    gamma.vec <- DTD.model
    model.esti.type <- estimate.c.type
  }

  if(is.list(new.data) && length(new.data) == 2){
    if(!"mixtures" %in%  names(new.data)){
      stop("In estimate.c.type: entry of 'new.data' must be named 'mixtures'.")
    }else{
      if(!is.matrix(new.data$mixtures)){
        stop("In estimate.c.type: 'new.data$mixtures' is not a matrix.")
      }
    }
    Y <- new.data$mixtures
  }else{
    if(!is.matrix(new.data)){
      stop("In estimate.c.type: 'new.data' is not a matrix.")
    }
    Y <- new.data
  }

  if (length(gamma.vec) != nrow(X) || nrow(X) != nrow(Y)) {
    stop("estimate_c: dimension of provided input (X, Y, g) does not match")
  }

  # check if the estimate.c.type of the model, fits the user input:
  # ('direct' model must not be deconvoluting with 'non-negative')
  if(model.esti.type != estimate.c.type){
    # default is 'decide.on.model' => no useless message
    if(estimate.c.type != "decide.on.model"){
      message(
        "in estimate_c:
        parameter 'estimate.c.type' not consistent with 'DTD.model$estimate.c.type'.
        Using 'DTD.model$estimate.c.type'"
        )
    }
  }
  # transform gamma.vec into a diagonal matrix (everything but diagonal is zero)
  Gamma <- diag(as.vector(gamma.vec))

  if(model.esti.type == "direct"){
    sol <- estimate_direct_c(
      X.matrix = X
      , Gamma = Gamma
      , new.data = Y
    )
  } else if (model.esti.type == "non_negative" ||
           model.esti.type == "non-negative" ){
    sol <- estimate_nn_c(
      X.matrix = X
      , new.data = Y
      , Gamma = Gamma
    )
  } else {
    stop(
      paste0(
        "In estimate_c: for 'estimate.c.type': ",
        model.esti.type,
        " there is no estimate_..._c() function"
        )
    )
  }
  return(sol)
}


