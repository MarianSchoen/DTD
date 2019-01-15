#' Descent generalized FISTA \cr
#' (fast iterative shrinkage thresholding algorithm)
#'
#' descent_generalized_fista takes as input a vector, a gradient function and an
#' evaluation function (and some additional parameters/functions). Then, it
#' iteratively minimizes the tweak vector via FISTA (Beck and Teboulle 2009). Basically,
#' the following equations are used:\cr
#' # prelimary initialization step\cr
#' for k = 2,..., maxit:\cr
#' \itemize{
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
#'}
#'
#' The gradient and evaluation function both take only one argument, the soft thresholding function two.
#' If your gradient, evaluation or soft thresholding function require more arguments, please write a wrapper.
#' Exemplary wrapper functions can be found in the examples.
#'
#' @param tweak.vec numeric vector, with which the minimization algorithm starts, defaults to NA
#' @param lambda float, regularization factor for ST.FUN function, defaults to 0
#' @param maxit integer, maximum number of iterations for the iterative
#'   minimization, defaults to 1000
#' @param learning.rate float, step size while learning. Defaults to NA. If it is NA, the learning rate will
#'  be estimated as published by Barzilai & Borwein 1988
#' @param F.GRAD.FUN function with one parameter: vector with same length as
#'   tweak.vec, and one return argument: vector with same length as tweak.vec
#'   For a given vector, it returns the gradient of the loss-function
#' @param ST.FUN function, with one parameter: vector with same length as
#'   tweak.vec, and one return argument: vector with same length as tweak.vec
#'   implementation of the proximal operator for the chosen Loss-function. Here
#'   named soft thresholding function. Defaults to soft_thresholding
#' @param FACTOR.FUN function, with a integer as input, and a float as return
#'   value. This function returns the factor the nesterov correction
#'   extrapolates. Defaults to nesterov_faktor
#' @param EVAL.FUN function, with one parameter: vector with same length as
#'   tweak.vec, and a single float as return. This function evaluates the loss
#'   function.
#' @param NORM.FUN function, with one parameter: After each step the tweak vector will be normed/scaled.
#' If no scaling should be done, set NORM.FUN to identity. Defaults to identity.
#' @param line.search.speed numeric, factor with which the learning rate changes,
#'  if the optimium has not been found. Defaults to 2.
#' @param cycles integer, in each iteration one gradient is calculated. To find the
#'  best step size, we do "cycles" steps, and evaluate each of them to find the best step size.
#'  Defaults to 50
#' @param save.all.tweaks boolean, should all tweak vectores during all iterations be stored. Defaults to FALSE
#' @param use.restart boolean, restart the algorithm if the update was not a descent step. Defaults to TRUE
#' @param verbose boolean, should information be printed to console. Defaults to TRUE
#' @param NESTEROV.FUN function, applied to the nesterov extrapolation.
#' In DTD it is all g must be positive. Therefore we set default to a positive_subspace_pmax wrapper function
#'
#' @export
#'
#' @return list, including [["Converge"]], [["Tweak"]] and [["Lambda"]]  (depending on "save_all_tweks": [["History"]])
#'
#' @examples
#' library(DTD)
#' random.data <- generate_random_data(n.types = 10,
#'                                     n.samples.per.type = 150,
#'                                     n.features = 250,
#'                                     sample.type = "Cell",
#'                                     feature.type = "gene")
#'
#' # normalize all samples to the same amount of counts:
#' normalized.data <- normalize_to_count(random.data)
#'
#' # extract indicator list.
#' # This list contains the Type of the sample as value, and the sample name as name
#' indicator.list <- gsub("^Cell[0-9]*\\.", "", colnames(random.data))
#' names(indicator.list) <- colnames(random.data)
#'
#' # extract reference matrix X
#' # First, decide which cells should be deconvoluted.
#' # Notice, in the mixtures there can be more cells than in the reference matrix.
#' include.in.X <- paste0("Type", 2:7)
#'
#' percentage.of.all.cells <- 0.2
#' sample.X <- sample_random_X(included.in.X = include.in.X,
#'                             pheno = indicator.list,
#'                             exp.data = normalized.data,
#'                             percentage.of.all.cells = percentage.of.all.cells)
#' X.matrix <- sample.X$X.matrix
#' samples.to.remove <- sample.X$samples.to.remove
#'
#' # all samples that have been used in the reference matrix, must not be included in
#' # the test/training set
#' remaining.mat <- random.data[, -which(colnames(random.data) %in% samples.to.remove)]
#' train.samples <- sample(x = colnames(remaining.mat),
#'                        size = ceiling(ncol(remaining.mat)/2),
#'                        replace = FALSE)
#' test.samples <- colnames(remaining.mat)[which(!colnames(remaining.mat) %in% train.samples)]
#'
#' train.mat <- remaining.mat[, train.samples]
#' test.mat <- remaining.mat[, test.samples]
#'
#' indicator.train <- indicator.list[names(indicator.list) %in% colnames(train.mat)]
#' training.data <- mix_samples(gene.mat = train.mat,
#'                              pheno = indicator.train,
#'                              included.in.X = include.in.X,
#'                              n.samples = 500,
#'                              n.per.mixture = 100,
#'                              verbose = FALSE)
#'
#' indicator.test <- indicator.list[names(indicator.list) %in% colnames(test.mat)]
#' test.data <-  mix_samples(gene.mat = test.mat,
#'                           pheno = indicator.test,
#'                           included.in.X = include.in.X,
#'                           n.samples = 500,
#'                           n.per.mixture = 100,
#'                           verbose = FALSE)
#'
#' # The descent_generalized_fista optimizer finds the minimum iteratively
#' # using accelareted gradient descent.
#' # Therefore a gradient, and an evaluation function must be provided.
#' # Within the fista implementation the gradient/evalution functions get evoked
#' # with only one parameter.
#' # All other parameters for the gradient (in the following example X, Y and C) must
#' # be set using default parameter
#' # This can be done using wrappers:
#'
#' # wrapper for gradient:
#' DTD.grad.wrapper <- function(tweak,
#'                              X = X.matrix,
#'                              Y = training.data$mixtures,
#'                              C = training.data$quantities){
#'    grad <- gradient_cor_trace(X = X, Y = Y, C = C, tweak = tweak)
#'    return(grad)
#' }
#' # wrapper for evaluate corelation:
#' DTD.evCor.wrapper <- function(tweak,
#'                              X = X.matrix,
#'                              Y = training.data$mixtures,
#'                              C = training.data$quantities){
#'    loss <- evaluate_cor(X = X, Y = Y, C = C, tweak = tweak)/ncol(X)
#'    return(loss)
#' }
#'
#' start.tweak <- rep(1, nrow(X.matrix))
#' start.cor <- DTD.evCor.wrapper(start.tweak)
#' cat("Starting correlation: ", -start.cor, "\n")
#'
#' catch <- descent_generalized_fista(tweak.vec = start.tweak,
#'                                    F.GRAD.FUN = DTD.grad.wrapper,
#'                                    ST.FUN = soft_thresholding,
#'                                    FACTOR.FUN = nesterov_faktor,
#'                                    EVAL.FUN = DTD.evCor.wrapper,
#'                                    line.search.speed = 2,
#'                                    maxit = 250,
#'                                    save.all.tweaks = TRUE,
#'                                    use.restart = TRUE,
#'                                    verbose = FALSE)
#'
#' print(ggplot_convergence(fista.output = catch,
#'                          main="additional title")
#'      )
descent_generalized_fista <- function(tweak.vec = NA,
                                      lambda=0,
                                      maxit=1e3,
                                      learning.rate = NA,
                                      F.GRAD.FUN,
                                      ST.FUN = soft_thresholding,
                                      FACTOR.FUN = nesterov_faktor,
                                      EVAL.FUN,
                                      NORM.FUN = identity,
                                      line.search.speed = 2,
                                      cycles=50,
                                      save.all.tweaks=FALSE,
                                      use.restart=TRUE,
                                      verbose=TRUE,
                                      NESTEROV.FUN = positive_subspace_pmax){
  # safety checks:
  if(any(is.na(tweak.vec))){
    stop("Tweak vector includes NAs")
  }
  if(!(is.numeric(lambda) || is.numeric(maxit) || is.numeric(line.search.speed) || is.numeric(cycles))){
    stop("Set lambda, maxit, line.search.speed and cycles to numeric values")
  }
  if(!(is.logical(save.all.tweaks) || is.logical(verbose) || is.logical(use.restart))){
    stop("Set save.all.tweaks, use.restart and verbose to logicals")
  }

  # check if the norm function changes EVAL value:
  if(!isTRUE(all.equal(EVAL.FUN(tweak.vec), EVAL.FUN(NORM.FUN(tweak.vec))))){
    stop("Norm function changes eval value.")
  }

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


  # initialise required variables:
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
    if(rev(converge_vec)[1] > eval){
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

    # add the last eval to the convergence_vec:
    converge_vec <- c(converge_vec, eval)

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
        nesterov.winner <- l.nesterov
      }
    }

    # if verbose = TRUE, print information to the screen
    if(verbose){
      cat("###############################################\n")
      cat("iter: ", iter, "\n")
      cat("Loss: ", rev(converge_vec)[1], "\n")
      cat("learning.rate: ", learning.rate, "\n")
      cat("cycle winner pos:", winner.pos, "\n")
      cat("factor: ", factor, "\n")
      cat("nesterov rate: ", nesterov.winner, "\n")
      # plot converge_vec:
      if(iter %% 100 == 0){
        graphics::plot(1:iter, converge_vec)
      }
    }
  }

  names(tweak.vec) <- tweak.names
  # build a list to return:
  ret <- list("Tweak"=NORM.FUN(tweak.vec), "Convergence"=converge_vec, "Lambda" = lambda)

  # if save.all.tweaks is TRUE add the tweak.history
  if(save.all.tweaks){
    ret$History <- tweak.history
  }
  return(ret)
}
