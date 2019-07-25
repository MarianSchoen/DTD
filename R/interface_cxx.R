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
solve_fista_goertler <- function(model, lambda = 0.01, maxiter = 100, save.all.tweaks ) {
  # check input params...
  check_model(model)
  if( ! (is.numeric(lambda) && lambda >= 0.0 )) {
    stop("lambda is not within range (numeric, >= 0)")
  }
  if( ! (is.numeric(maxiter) && maxiter %% 1 == 0 && maxiter >= 2 )) {
    stop("maxiter is not within range (integer, >= 2)")
  }

  result <- .Call("_dtd_solve_fista_goertler", model, lambda, integer(maxiter), save.all.tweaks, PACKAGE="DTD")

  # take over row and colnames from model:
  names(result$Tweak) <- names(model$tweak)
  rownames(result$History) <- rownames(model$tweak)
  return(result)
}
evaluate_model <- function(model) {
  check_model(model)
  return(.Call("_dtd_evaluate_model_goertler", model, PACKAGE = "DTD"))
}
