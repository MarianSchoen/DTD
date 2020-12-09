#' Descent generalized FISTA \cr
#' (fast iterative shrinkage thresholding algorithm)
#'
#' descent_generalized_fista takes as input a vector, a gradient function
#' and an evaluation function (and some additional parameters/functions).
#' Then, it iteratively minimizes the tweak vector via FISTA
#' (Beck and Teboulle 2009).
#' Basically,
#' the following equations are used:\cr
#' # prelimary initialization step\cr
#' for k = 2,..., maxit:\cr
#' \itemize{
#'     \item y_vec = NORM.FUN(y_vec)
#'     \item grad = F.GRAD.FUN(y_vec)\cr
#'     \item# find best step size between 0 and learning.rate\cr
#'          for step.size = 0,..., learning.rate:\cr
#'     \itemize{
#'         \item u = ST.FUN(y_vec - step.size * grad, step.size * lambda)\cr
#'         \item eval = EVAL.FUN(u)
#'     }
#'     \item # only keep u with minimal eval\cr
#'     \item tweak.old = tweak.vec\cr
#'     \item # if it was a descent step: \cr tweak.vec = u
#'     \item # if it was not a descent step: \cr tweak.vec = tweak.vec \cr#(and restart as suggested in  O’Donoghue & Candes (2012))
#'     \item # Nesterov extrapolation:
#'     \item nesterov.direction = tweak.vec - tweak.old
#'     \item # find best extrapolation size bewteen 0 and FACTOR.FUN(k):\cr
#'          for ne = 0 ,... FACTOR.FUN(k):\cr
#'     \itemize{
#'         \item y_vec = u_vec + ne * nesterov.direction
#'         \item eval = EVAL.FUN(y_vec)
#'     }
#'     \item # only keep y_vec with minimal eval
#'     \item stop criterion: if the descent in either the gradient
#'     or the nesterov step was below a critical value, then stop.
#'}
#'
#' The gradient and evaluation function both take only one argument,
#' the soft thresholding function two. If your gradient, evaluation or soft
#' thresholding function require more arguments, please write a wrapper.
#'
#' @param tweak.vec numeric vector, starting vector for the DTD algorithm
#' @param lambda non-negative float, regularization factor for ST.FUN function.
#' @param maxit integer, maximum number of iterations for the iterative
#' minimization.
#' @param learning.rate float, step size during optimization. If it is NA,
#' the learning rate will be estimated as published by Barzilai & Borwein 1988.
#' Notice, the algorithm adjusts the learning rate during optimization.
#' @param F.GRAD.FUN function with one parameter: vector with same length as
#'   tweak.vec, and one return argument: vector with same length as tweak.vec
#'   For a given vector, it returns the gradient of the loss-function
#' @param ST.FUN function, with one parameter: vector with same length as
#'   tweak.vec, and one return argument: vector with same length as tweak.vec.
#'   Implementation of the proximal operator for the chosen Loss-function. Here
#'   named soft thresholding function.
#' @param FACTOR.FUN function, with a integer as input, and a float as return
#'   value. This function returns the factor of the nesterov correction
#'   extrapolates.
#' @param EVAL.FUN function, with one parameter: vector with same length as
#'   tweak.vec, and a single float as return. This function evaluates the loss
#'   function.
#' @param NORM.FUN function, with one parameter: After each step the
#' tweak vector will be normed/scaled. If no scaling should be done,
#' set NORM.FUN to identity.
#' @param line.search.speed numeric, factor with which the learning rate
#' changes during optimization. If 'line.search.speed' is, e.g.,
#' 2, the learning rate gets doubled if the highest 'cycle' led to the best
#' eval score. If the 'cycle = 0' led to the best eval score, it would get
#' halved.
#' @param cycles integer, in each iteration one gradient is calculated.
#'  To find the best step size, we do "cycles" steps, and evaluate each of
#'  them to find the best step size.
#' @param save.all.tweaks logical, should all tweak vectores during all
#' iterations be stored.
#' @param use.restart logical, restart the algorithm if the update was
#' not a descent step.
#' @param verbose logical, should information be printed to console?
#' @param NESTEROV.FUN function, applied to the nesterov extrapolation.
#' @param stop.crit.threshold numeric value. The change in either the gradient
#' or nesterov step has to be at least 'stop.crit.threshold', or the algorithm
#' will stop.
#'
#' @return list, including
#' \itemize{
#'    \item "Converge", numeric vector, EVAL.FUN in every step
#'    \item "Tweak" numeric vector, the best distinguishing vector
#'    \item "Lambda", numeric value, used \eqn{\lambda}
#'    \item depending on "save.all.tweaks": "History", numeric matrix, "Tweak" vector of every step
#' }
descent_generalized_fista <- function(tweak.vec,
                                      lambda=0,
                                      maxit=5e2,
                                      learning.rate = NA,
                                      F.GRAD.FUN,
                                      ST.FUN = "softmax",
                                      FACTOR.FUN = "nesterov_factor",
                                      EVAL.FUN,
                                      NORM.FUN = "norm2",
                                      line.search.speed = 2,
                                      cycles=5,
                                      save.all.tweaks=FALSE,
                                      use.restart=TRUE,
                                      verbose=FALSE,
                                      NESTEROV.FUN = "positive",
                                      stop.crit.threshold = 1e-13
                                      ){


  # safety check: tweak.vec
  test <- test_tweak_vec(tweak.vec = tweak.vec,
                         output.info = c("descent_generalized_fista", "tweak.vec"))
  # end -> tweak

  # safety check: F.GRAD.FUN
  # I am going to check if the function returns an error, if called with only 1 parameter
  # I am NOT going to check if the function returns the steepest direction (or at least a descent direction)
  tmp <- try(F.GRAD.FUN(tweak.vec), silent = TRUE)
  if(any(grepl(pattern = "Error ", x = tmp))){
    stop("In descent_generalized_fista: 'F.GRAD.FUN' produces an error, if called with tweak.vec as first, and only argument")
  }
  # end -> F.GRAD.FUN

  # safety check: EVAL.FUN
  # I am going to check if the function returns an error, if called with only 1 parameter
  tmp <- try(EVAL.FUN(tweak.vec), silent = TRUE)
  if(any(grepl(pattern = "Error ", x = tmp))){
    stop("In descent_generalized_fista: 'EVAL.FUN' produces an error, if called with tweak.vec as first, and only argument")
  }
  # end -> EVAL.FUN

  # safety check: ST.FUN
  if( is.character(ST.FUN) && ST.FUN == "softmax" ) {
    ST.FUN = soft_thresholding
  } else if (is.character(ST.FUN)) {
    stop("unknown thresholding function. may currently be just \"softmax\" or some user defined function.")
  }
  # I am going to check if the function returns an error, if called with only 1 parameter
  tmp <- try(ST.FUN(x = tweak.vec,
                    lambda = lambda), silent = TRUE)
  if(any(grepl(pattern = "Error ", x = tmp))){
    stop("In descent_generalized_fista: 'ST.FUN' produces an error, if called with 'tweak.vec' as first, and 'lambda' as second argument")
  }
  # end -> ST.FUN

  # safety check: FACTOR.FUN
  if( is.character(FACTOR.FUN) && FACTOR.FUN == "nesterov_factor" ) {
    FACTOR.FUN = nesterov_factor
  } else if (is.character(FACTOR.FUN)) {
    stop("unknown factor function. may currently be just \"nesterov_factor\" or some user defined function.")
  }
  # I am going to check if the function returns an error, if called with only 1 parameter
  tmp <- try(FACTOR.FUN(maxit), silent = TRUE)
  if(any(grepl(pattern = "Error ", x = tmp))){
    stop("In descent_generalized_fista: 'FACTOR.FUN' produces an error, if called with 'maxit' as first, and only argument")
  }
  # end -> FACTOR.FUN

  # safety check: NORM.FUN
  if( is.character(NORM.FUN) && NORM.FUN == "identity" ) {
    NORM.FUN = identity
  } else if( is.character(NORM.FUN) && NORM.FUN == "norm1") {
    NORM.FUN = n1normed
  } else if( is.character(NORM.FUN) && NORM.FUN == "norm2") {
    NORM.FUN = n2normed
  } else if (is.character(NORM.FUN)) {
    stop("unknown norm function. may currently be \"identity\", \"norm1\", \"norm2\", or some other user defined function.")
  }
  # I am going to check if the function returns an error, if called with only 1 parameter
  tmp <- try(NORM.FUN(tweak.vec), silent = TRUE)
  if(any(grepl(pattern = "Error ", x = tmp))){
    stop("In descent_generalized_fista: 'NORM.FUN' produces an error, if called with 'tweak.vec' as first, and only argument")
  }
  # end -> NORM.FUN

  # safety check: NESTEROV.FUN
  if( is.character(NESTEROV.FUN) && NESTEROV.FUN == "positive" ) {
    NESTEROV.FUN = positive_subspace_pmax
  } else if (is.character(NESTEROV.FUN)) {
    stop("unknown nesterov function. may currently be just \"positive\" or some user defined function.")
  }
  # I am going to check if the function returns an error, if called with only 1 parameter
  tmp <- try(NESTEROV.FUN(tweak.vec), silent = TRUE)
  if(any(grepl(pattern = "Error ", x = tmp))){
    stop("In descent_generalized_fista: 'NESTEROV.FUN' produces an error, if called with 'tweak.vec' as first, and only argument")
  }
  # end -> NORM.FUN

  # safety check: stop.crit.threshold
  test <- test_numeric(test.value = stop.crit.threshold,
                       output.info = c("descent_generalized_fista", "stop.crit.threshold"),
                       min = -Inf,
                       max = Inf)
  # end -> stop.crit.threshold

  # safety check: lambda
  test <- test_numeric(test.value = lambda,
                       output.info = c("descent_generalized_fista", "lambda"),
                       min = 0,
                       max = Inf)
  # end -> lambda

  # safety check: maxit
  test <- test_integer(test.value = maxit,
                       output.info = c("descent_generalized_fista", "maxit"),
                       min = 2,
                       max = Inf)
  # end -> maxit

  # safety check: line.search.speed
  test <- test_numeric(test.value = line.search.speed,
                       output.info = c("descent_generalized_fista", "line.search.speed"),
                       min = .Machine$double.eps,
                       max = Inf)
  # end -> line.search.speed

  # safety check: cycles
  test <- test_integer(test.value = cycles,
                       output.info = c("descent_generalized_fista", "cycles"),
                       min = 1,
                       max = Inf)
  # end -> cycles

  # safety check: save.all.tweaks
  test <- test_logical(test.value = save.all.tweaks,
                       output.info = c("descent_generalized_fista", "save.all.tweaks"))
  # end -> save.all.tweaks

  # safety check: use.restart
  test <- test_logical(test.value = use.restart,
                       output.info = c("descent_generalized_fista", "use.restart"))
  # end -> use.restart

  # safety check: verbose
  test <- test_logical(test.value = verbose,
                       output.info = c("descent_generalized_fista", "verbose"))
  # end -> verbose

  # check if the norm function changes EVAL value:
  if(!isTRUE(all.equal(EVAL.FUN(tweak.vec), EVAL.FUN(NORM.FUN(tweak.vec))))){
    stop("descent_generalized_fista: Norm function changes eval value.")
  }

  # norm the input tweak.vec
  tweak.vec <- NORM.FUN(tweak.vec)

  # If no learning.rate is set, the initial learning rate will be initialized according to:
  # Barzilai & Borwein 1988
  # It estimates via:
  # learning.rate <- <delta(tweak), delta(g)> / <delta(g), delta(g)>
  # here, < , > denotes the scalar product.
  # delta(x) := x_(k) - x_(k-1)
  if(is.na(learning.rate)){
    grad <- F.GRAD.FUN(tweak.vec)
    norm2 <- sqrt(sum(grad**2))
    t <- 1/norm2
    tweak_hat <- tweak.vec - t * grad
    g_hat <- F.GRAD.FUN(tweak_hat)
    learning.rate <- abs((tweak.vec - tweak_hat) %*% (grad - g_hat) / sum((grad - g_hat)**2))
    # learning rate is a matrix => make it a numeric:
    learning.rate <- as.numeric(learning.rate)
    l0 <- learning.rate
  }


  # initialize required variables:
  tweak_old <- tweak.vec
  converge_vec <- EVAL.FUN(tweak.vec)
  y_vec <- tweak.vec
  nesterov.counter <- 2
  factor <- FACTOR.FUN(nesterov.counter)

  if(verbose){
    cat("starting to optimize \n")
  }
  # store the names of tweak, that they can be reset at the end:
  tweak.names <- names(tweak.vec)


  # Notice, if save.all.tweaks is FALSE only the last tweak.vec will be stored,
  # and returned after optimization. if save.all.tweaks is true all tweak.vec will be stored
  # in tweak.history, and returned after optimization:
  if(save.all.tweaks){
    tweak.history <- matrix(NA, nrow = length(tweak.vec), ncol = maxit)
    tweak.history[, 1] <- NORM.FUN(tweak.vec)
    rownames(tweak.history) <- tweak.names
  }

  # Every gradient step is done multiple times with different step sizes larging from 0 to learning.rate
  # sequence holds the different learning rates:
  sequence <- seq(from = 0, to = learning.rate, length.out = cycles)

  # Notice that the for loop starts at 2, due to the extrapolation/correction step of the FISTA algorithm
  for(iter in 2:maxit){
    y_vec <- NORM.FUN(y_vec)

    # calculate gradient at the current position:
    grad <- F.GRAD.FUN(y_vec)

    # The following code junk is hard to understand.
    # In the end, u_mat will be a matrix with new u_vecs in each row.
    # This is done by applying the ST.FUN on the sequence of step sizes.
    # The return of lapply needs to be unlisted and converted to a matrix.
    u_mat <- matrix(
      unlist(
        lapply(
          sequence,
          function(x){ST.FUN(y_vec - x * grad, x*lambda)}),
        use.names = FALSE),
      nrow = cycles,
      byrow = TRUE)

    # every row of u_mat holds a u_vec with another step.size.
    # In order to find the best of them, we use the EVAL.FUN on all of them:
    eval.vec <- apply(u_mat, 1, EVAL.FUN)

    # and find the winner, which is the minimum:
    winner.pos <- which.min(eval.vec)

    # Now we found the best step size between 0 and learning.rate.
    # If the u_vec with step.size = learning.rate (so the last entry) is the winner,
    # the step size needs to be increased.
    if(winner.pos == cycles){
      learning.rate <- learning.rate * line.search.speed
      # calculate new sequence, as learning.rate has changed
      sequence <- seq(from = 0, to = learning.rate, length.out = cycles)
    }
    # If the u_vec with step.size = 0 (so the first entry) is the winner,
    # the step size needs to be decreased.
    if(winner.pos == 1){
      learning.rate <- learning.rate / line.search.speed
      # calculate new sequence, as learning.rate has changed
      sequence <- seq(from = 0, to = learning.rate, length.out = cycles)
    }


    # set the winning u_vec, and eval. Norm the u_vec using the provided function:
    u_vec <- NORM.FUN(u_mat[winner.pos, ])
    eval <- eval.vec[winner.pos]

    # update tweak_old (tweak of last iteration)
    tweak_old <- tweak.vec

    # check if the last step was a descent step:
    if((rev(converge_vec)[1] - eval) > 0){
      nesterov.counter <- nesterov.counter + 1
      # if it was a descent step, update tweak.vec
      tweak.vec <- u_vec
    }else{
      # if not, tweak.vec get's not updated, but eval
      eval <- EVAL.FUN(tweak.vec)
      # according to (O'Donoghue and Candes 2012) restart the nesterov.counter (=> restart)
      if(use.restart){
        nesterov.counter <- 2
      }
    }

    # store the loss-change before nesterov
    change.last.iter.before.nesterov <- rev(converge_vec)[1] - eval

    # and if set, update tweak.history
    if(save.all.tweaks){
      tweak.history[, iter] <- tweak.vec
    }


    # the same line search approach gets applied to the nesterov extrapolation:
    factor <- FACTOR.FUN(nesterov.counter)
    nesterov.sequence <- seq(from=0, to=factor, length.out = cycles)
    nesterov.direction <- tweak.vec - tweak_old
    best.y <- Inf
    for(l.nesterov in nesterov.sequence){
      # extrapolate by "l.nesterov" ...
      # Notice that the NESTEROV.FUN is needed to restrict the solution to a subspace
      # (e.g. for lfl DTD no negativ entries are allowed, therefore NESTEROV.FUN sets the negative entries to 0)
      y_vec.l <- NESTEROV.FUN(tweak.vec + l.nesterov * nesterov.direction)
      # ... check the eval.y value ...
      eval.y <- EVAL.FUN(y_vec.l)
      if(eval.y < best.y){
        # ... and reset the best.y if its a new winner
        y_vec <- y_vec.l
        best.y <- eval.y
        nesterov.winner.pos <- which(nesterov.sequence == l.nesterov)
      }
    }

    # add the last eval to the convergence_vec:
    converge_vec <- c(converge_vec, best.y)

    # find stop criterion
    change.last.iter.after.nesterov <- eval - best.y

    change.last.iter <- max(change.last.iter.after.nesterov, change.last.iter.before.nesterov)
    ###


    ### end stop
    # if verbose = TRUE, print information to the screen
    if(verbose){
      cat("###############################################\n")
      cat("iter: ", iter, "\n")
      cat("Loss: ", rev(converge_vec)[1], "\n")
      cat("learning.rate: ", learning.rate, "\n")
      cat("cycle winner pos:", winner.pos, "\n")
      cat("nesterov extrapolation: ", factor, "\n")
      # plot converge_vec:
      if(iter %% 100 == 0){
        graphics::plot(1:iter, log(-converge_vec))
      }
      # stop crit output
      cat("change of L, only gradient step: ", change.last.iter.before.nesterov, "\n")
      cat("change of L, only nesterov step: ", change.last.iter.after.nesterov, "\n")
      cat("change of L: ", change.last.iter, "\n")
    }


    if(
      change.last.iter <= stop.crit.threshold && # converged
      !change.last.iter == 0
      # learning.rate <= stop.crit.threshold
      ){
      break
    }

  }


  # stop criterion
  names(tweak.vec) <- tweak.names
  # build a list to return:
  ret <- list("Tweak"=NORM.FUN(tweak.vec),
              "Convergence"=converge_vec,
              "Lambda" = lambda)

  # if save.all.tweaks is TRUE add the tweak.history
  if(save.all.tweaks){
    ret$History <- tweak.history[, 1:iter]
  }
  return(ret)
}
#' Descent generalized FISTA \cr
#' (fast iterative shrinkage thresholding algorithm)
#'
#' descent_generalized_fista takes as input a vector, a gradient function
#' and an evaluation function (and some additional parameters/functions).
#' Then, it iteratively minimizes the tweak vector via FISTA
#' (Beck and Teboulle 2009).
#' Basically,
#' the following equations are used:\cr
#' # prelimary initialization step\cr
#' for k = 2,..., maxit:\cr
#' \itemize{
#'     \item y_vec = NORM.FUN(y_vec)
#'     \item grad = F.GRAD.FUN(y_vec)\cr
#'     \item# find best step size between 0 and learning.rate\cr
#'          for step.size = 0,..., learning.rate:\cr
#'     \itemize{
#'         \item u = ST.FUN(y_vec - step.size * grad, step.size * lambda)\cr
#'         \item eval = EVAL.FUN(u)
#'     }
#'     \item # only keep u with minimal eval\cr
#'     \item tweak.old = tweak.vec\cr
#'     \item # if it was a descent step: \cr tweak.vec = u
#'     \item # if it was not a descent step: \cr tweak.vec = tweak.vec \cr#(and restart as suggested in  O’Donoghue & Candes (2012))
#'     \item # Nesterov extrapolation:
#'     \item nesterov.direction = tweak.vec - tweak.old
#'     \item # find best extrapolation size bewteen 0 and FACTOR.FUN(k):\cr
#'          for ne = 0 ,... FACTOR.FUN(k):\cr
#'     \itemize{
#'         \item y_vec = u_vec + ne * nesterov.direction
#'         \item eval = EVAL.FUN(y_vec)
#'     }
#'     \item # only keep y_vec with minimal eval
#'     \item stop criterion: if the descent in either the gradient
#'     or the nesterov step was below a critical value, then stop.
#'}
#'
#' @param model a model as constructed in interface_cxx.R
#' @param lambda non-negative float, regularization factor for ST.FUN function.
#' @param maxit integer, maximum number of iterations for the iterative
#' minimization.
#' @param stop.crit.threshold numeric value. The change in either the gradient
#' or nesterov step has to be at least 'stop.crit.threshold', or the algorithm
#' will stop.
#' @param save.all.tweaks logical, should all tweak vectores during all
#' iterations be stored.
#' @param learning.rate float, step size during optimization. If it is NA,
#' the learning rate will be estimated as published by Barzilai & Borwein 1988.
#' Notice, the algorithm adjusts the learning rate during optimization.
#' @param line.search.speed numeric, factor with which the learning rate
#' changes during optimization. If 'line.search.speed' is, e.g.,
#' 2, the learning rate gets doubled if the highest 'cycle' led to the best
#' eval score. If the 'cycle = 0' led to the best eval score, it would get
#' halved.
#' @param cycles integer, in each iteration one gradient is calculated.
#'  To find the best step size, we do "cycles" steps, and evaluate each of
#'  them to find the best step size.
#' @param use.restart logical, restart the algorithm if the update was
#' not a descent step.
#' @param verbose logical, if set to true, will output information during
#' iteration.
#'
#' @return a list that contains the trained model and its History
descent_generalized_fista_cxx <- function(model,
                                          lambda = 0.01,
                                          maxit = 500,
                                          stop.crit.threshold = 1e-13,
                                          stop.crit.navg = 30,
                                          save.all.tweaks = FALSE,
                                          learning.rate = NA,
                                          line.search.speed = 2.0,
                                          cycles = 5,
                                          use.restart = TRUE,
                                          verbose = FALSE,
                                          ...) {
  return(solve_fista_goertler(model, lambda, maxit, stop.crit.threshold, stop.crit.navg, save.all.tweaks, learning.rate, line.search.speed, cycles, use.restart, verbose))
}
