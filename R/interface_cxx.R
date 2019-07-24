solve_fista_goertler <- function(model, lambda = 0.01, maxiter = 100) {
  # check input params...
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
  if( ! (is.numeric(lambda) && lambda >= 0.0 )) {
    stop("lambda is not within range (numeric, >= 0)")
  }
  if( ! (is.numeric(maxiter) && maxiter %% 1 == 0 && maxiter >= 2 )) {
    stop("maxiter is not within range (integer, >= 2)")
  }

  cat("bla: ", lambda, maxiter)

  result <- .Call("_dtd_solve_fista_goertler", model, lambda, maxiter, PACKAGE="DTD")
  return(result)

  ## names(tweak.vec) <- tweak.names
  ##                                       # build a list to return:
  ## ret <- list("Tweak"=NORM.FUN(tweak.vec),
  ##             "Convergence"=converge_vec,
  ##             "Lambda" = lambda)

                                        # if save.all.tweaks is TRUE add the tweak.history
  # TODO: not implemented
  #if(save.all.tweaks){
  #  ret$History <- tweak.history
  #}
}
