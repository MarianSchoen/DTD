#' Title
#'
#' @param model
#'
#' @return
#' @export
#'
#' @examples
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
                                        # size checking
  if( ! (nrow(model$X) == nrow(model$Y) &&
         ncol(model$X) == nrow(model$C) &&
         ncol(model$C) == ncol(model$Y) &&
         nrow(model$X) == length(model$tweak) ) ) {
    stop("input matrices have incompatible sizes.")
  }
}
#' Title
#'
#' @param model
#' @param lambda
#' @param maxiter
#' @param save.all.tweaks
#'
#' @return
#' @export
#' @useDynLib DTD, .registration = TRUE
#'
#' @examples
solve_fista_goertler <- function(model, lambda = 0.01, maxiter = 100, save.all.tweaks = FALSE, learningrate = NA, linesearchspeed = 2.0, cycles = 5, restarts = TRUE) {
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
  if( result$Valid == FALSE ) {
    message("result is INVALID!")
  }
  return(result)
}
#' Title
#'
#' @param model
#'
#' @return
#' @export
#'
#' @examples
evaluate_model <- function(model) {
  check_model(model)
  return(.Call(
    "_dtd_evaluate_model_goertler",
    model,
    PACKAGE = "DTD"
    ))
}
