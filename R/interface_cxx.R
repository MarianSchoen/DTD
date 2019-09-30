#' check_model
#'
#' @param model list to be checked if it conforms with the model needed in solve_fista_goertler
#'
#' @return nothing if everything is okay
#'
check_model <- function(model) {
  if( ! is.list(model) ) {
    stop("model is not a list.")
  }
  if( ! ( "X" %in% names(model) &&
          "Y" %in% names(model) &&
          "C" %in% names(model) &&
          "tweak" %in% names(model) )) {
    stop("list \"model\" does not contain elements X, Y, C and tweak");
  }
  if( ! ( "normfnid" %in% names(model) &&
          "threshfnid" %in% names(model) &&
          "subspfnid" %in% names(model)) &&
          "coeffestim" %in% names(model)
     ) {
    stop("list \"model\" does not contain function enums (normfnid, threshfnid, subspfnid)");
  }
  # make sure these are integers:
  model$normfnid <- as.integer(model$normfnid)
  model$threshfnid <- as.integer(model$threshfnid)
  model$subspfnid <- as.integer(model$subspfnid)
  model$coeffestim <- as.integer(model$coeffestim)
                                        # size checking
  if( ! (nrow(model$X) == nrow(model$Y) &&
         ncol(model$X) == nrow(model$C) &&
         ncol(model$C) == ncol(model$Y) &&
         nrow(model$X) == length(model$tweak) ) ) {
    stop("input matrices have incompatible sizes.")
  }
}
#' empty_model
#'
#' constructs an empty model.
#'
#' @return the newly constructed model
#'
empty_model <- function() {
  model <- list()
  model$normfnid <- as.integer(0)
  model$threshfnid <- as.integer(0)
  model$subspfnid <- as.integer(0)
  model$coeffestim <- as.integer(0)
  return(model)
}
#' set_model_normfunction
#'
#' @param model input model, as constructed by, e.g., empty_model()
#' @param normfnname the desired value for the norm function. May currently be either "identity" or "norm2"
#'
#' @return the updated model with the norm function set.
#'
set_model_normfunction <- function(model, normfnname) {
  if( normfnname == 'IDENTITY' || normfnname == 'identity') {
    model$normfnid <- as.integer(0)
  } else if( normfnname == 'NORM2' || normfnname == 'norm2' ) {
    model$normfnid <- as.integer(1)
  } else {
    stop("invalid or unimplemented norm function.")
  }
  return(model)
}
#' set_model_subspacefunction
#'
#' @param model input model, as constructed by, e.g., empty_model()
#' @param subspfnname the desired value for the subspace function. May currently be only "positive".
#'
#' @return the updated model with the subspace function set.
#'
set_model_subspacefunction <- function(model, subspfnname) {
  if( subspfnname == 'POSITIVE' || subspfnname == 'positive' || subspfnname == '+' ) {
    model$subspfnid <- as.integer(0)
  } else {
    stop("invalid or unimplemented subsp function.")
  }
  return(model)
}
#' set_model_threshfunction
#'
#' @param model input model, as constructed by, e.g., empty_model()
#' @param threshfnname the desired value for the threshold function. May currently be only "softmax".
#'
#' @return the updated model with the thresh function set.
#'
set_model_threshfunction <- function(model, threshfnname) {
  if( threshfnname == 'SOFTMAX' || threshfnname == 'softmax' ) {
    model$threshfnid <- as.integer(0)
  } else {
    stop("invalid or unimplemented thresh function.")
  }
  return(model)
}
#' set_model_coeff_estimation
#'
#' @param model input model, as constructed by, e.g., empty_model()
#' @param coeffestimname the desired value for the coefficient estimation function. May currently be either "direct" or "nnls" (for non-negative least squares).
#'
#' @return the updated model with the coefficient estimation function set.
#'
set_model_coeff_estimation <- function(model, coeffestimname) {
  if( coeffestimname == 'direct' || coeffestimname == 'DIRECT' ) {
    model$coeffestim <- as.integer(0)
  } else if( coeffestimname == 'nnls' || coeffestimname == 'NNLS' || coeffestimname == 'non_negative' ) {
    model$coeffestim <- as.integer(1)
  } else {
    stop("invalid or unimplemented thresh function.")
  }
  return(model)
}
#' solve_fista_goertler
#'
#' @param model input model, as constructed by, e.g., empty_model()
#' @param lambda regularization parameter lambda
#' @param maxiter maximum number of iterations
#' @param stop.crit.threshold stopping criterion: if loss functions falls by less than this value, stop optimization.
#' @param save.all.tweaks if set to true, keep results after every iteration
#' @param learningrate initial learning rate (if set to NA, learning rate is determined automatically)
#' @param linesearchspeed rate of in- or decrease of learning rate
#' @param cycles number of function evaluations per step.
#' @param restarts if set to true, restart nesterov step width if loss fn is increased.
#'
#' @return list, containing an element "Tweak" with the result vector and, if save.all.tweaks is set to TRUE, "History", containing tweaks at every step.
#' @export
#' @useDynLib DTD, .registration = TRUE
#'
solve_fista_goertler <- function(model, lambda = 0.01, maxiter = 100, stop.crit.threshold = 1e-5, save.all.tweaks = FALSE, learningrate = NA, linesearchspeed = 2.0, cycles = 5, restarts = TRUE) {
  # check input params...
  check_model(model)
  if( ! (is.numeric(lambda) && lambda >= 0.0 )) {
    stop("lambda is not within range (numeric, >= 0)")
  }
  if( ! (is.numeric(maxiter) && maxiter %% 1 == 0 && maxiter >= 2 )) {
    stop("maxiter is not within range (integer, >= 2)")
  }

  result <- .Call(
    "_dtd_solve_fista_goertler",
    as.list(model),
    as.double(lambda),
    as.integer(maxiter),
    as.double(stop.crit.threshold ),
    as.logical(save.all.tweaks),
    learningrate,
    as.double(linesearchspeed),
    as.integer(cycles),
    as.logical(restarts),
    (! is.na(learningrate)),
    PACKAGE="DTD"
  )

  # take over row and colnames from model:
  names(result$Tweak) <- names(model$tweak)
  if( "History" %in% names(result) ) {
    rownames(result$History) <- names(model$tweak)
  }
  return(result)
}
#' cxx_evaluate_model
#'
#' returns loss function, evaluated on the given model
#'
#' @param model
#'
#' @return loss function
#' @export
#' @useDynLib DTD, .registration = TRUE
#'
cxx_evaluate_model <- function(model) {
  check_model(model)
  return(.Call(
    "_dtd_evaluate_model_goertler",
    model,
    PACKAGE = "DTD"
    ))
}
#' cxx_estimate_c
#'
#' returns estimated phenotype composition, given a trained model. Call cxx implementation.
#'
#' @param model
#'
#' @return composition (C)
#' @export
#' @useDynLib DTD, .registration = TRUE
#'
cxx_estimate_c <- function(model) {
  check_model(model)
  return(.Call(
    "_dtd_estimate_c",
    model,
    PACKAGE = "DTD"
    ))
}
#' NNLS
#'
#' @param A
#' @param b
#'
#' @return x
#'
#' @export
#' @useDynLib DTD, .registration = TRUE
#'
cxx_nnls <- function(A, b) {
  m = nrow(A)
  n = ncol(A)
  if( length(b) != m ) {
    stop("b and A have incompatible sizes.")
  }
  result <- .Call("_nnls",
                  as.matrix(A),
                  as.vector(b)
                  )
  return(result)
}
