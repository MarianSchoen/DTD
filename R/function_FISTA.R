#' Descent generalized FISTA \cr
#' (fast iterative shrinkage thresholding algorithm)
#'
#' descent_generalized_fista takes as input a vector, a gradient function and an
#' evaluation function (and some additional parameters/functions). Then, it
#' iteratively minimizes the tweak vector via the FISTA algorithms. Basically,
#' the following equations are used:
#'
#'
#' # wrapper are needed
#'
#' @param tweak_vec numeric vector, with which the minimization algorithm starts
#' @param lambda float, regularization factor for ST.FUN function
#' @param maxit integer, maximal number of iterations for the iterative
#'   minimization
#' @param learning.rate float, step size while learning
#' @param F.GRAD.FUN function with one parameter: vector with same length as
#'   tweak_vec, and one return argument: vector with same length as tweak_vec
#'   For a given vector, it returns the gradient of the loss-function
#' @param ST.FUN function, with one parameter: vector with same length as
#'   tweak_vec, and one return argument: vector with same length as tweak_vec
#'   implementation of the proximal operator for the chosen Loss-function. Here
#'   named soft thresholding function.
#' @param FACTOR.FUN function, with a integer as input, and a float as return
#'   value. This function returns the factor the nesterov correction
#'   extrapolates.
#' @param EVAL.FUN fucntion, with one parameter: vector with same length as
#'   tweak_vec, and a single float as return. This function evaluates the loss
#'   function.
#' @param line_search_speed ???
#' @param tol ???
#' @param verbose boolean, should information be printed to console
#'
#' @return list, including [["Converge"]] and [["Tweak"]]
#' @export
#'
#' @examples
#' library(DTD)
#' random.data <- generate.random.data(nTypes = 5,
#'                                     nSamples.perType = 10,
#'                                     nFeatures = 100,
#'                                     sample.type = "Cell",
#'                                     feature.type = "gene")
#'
#' # normalize all samples to the same amount of counts:
#' random.data <- normalizeToCount(random.data)
#'
#' # extract indicator list.
#' # This list contains the Type of the sample as value, and the sample name as name
#' indicator.list <- gsub("^Cell([0-9])*.", "", colnames(random.data))
#' names(indicator.list) <- colnames(random.data)
#'
#' # extract reference matrix X
#' # First, decide which cells should be deconvoluted.
#' # Notice, in the mixtures there can be more cells than in the reference matrix.
#' include.in.X <- c("Type2", "Type3", "Type4", "Type5")
#'
#' X.matrix <- matrix(NA, nrow=nrow(random.data), ncol=length(include.in.X))
#' colnames(X.matrix) <- include.in.X
#' rownames(X.matrix) <- rownames(random.data)
#'
#' percentage.of.all.cells <- 0.2
#'
#' # samples that are included in X must not be used in the training set!
#' samples.to.remove <- c()
#'
#' for(l.type in include.in.X){
#'   all.of.type <- names(indicator.list)[which(indicator.list == l.type)]
#'   chosen.for.X <- sample(x = all.of.type,
#'                          size = length(all.of.type) * percentage.of.all.cells,
#'                          replace = FALSE)
#'   samples.to.remove <- c(samples.to.remove, chosen.for.X)
#'
#'   average <- rowSums(random.data[, samples.to.remove])
#'   X.matrix[, l.type] <- average
#' }
#'
#' # here, I declare "Type1" as Tumor cells, and all other as immune cells
#' special.samples <- c("Type1")
#' all.samples <- unique(indicator.list)
#' sample.names <- all.samples[- which(all.samples %in% special.samples)]
#'
#' training.data <- mix.samples.Jitter(sample.names = sample.names,
#'                                      special.samples = special.samples,
#'                                      nMixtures = 1e3,
#'                                      datamatrix = random.data,
#'                                      indicator = indicator.list,
#'                                      singleSpecial = F,
#'                                      add_jitter = T,
#'                                      chosen.mean = 1,
#'                                      chosen.sd = 0.05,
#'                                      min.amount.samples = 1,
#'                                      verbose = FALSE,
#'                                      included.in.X = include.in.X)
#'
#' # wrapper for gradient:
#' DTD.grad.wrapper <- function(tweak){
#'    X <- X.matrix
#'    Y <- training.data$mixtures
#'    C <- training.data$quantities
#'    grad <- Trace.H.gradient(X = X, Y = Y, C = C, gamma.vec = tweak)
#'    return(grad)
#' }
#' # wrapper for evaluate corelation:
#' DTD.evCor.wrapper <- function(tweak){
#'    X <- X.matrix
#'    Y <- training.data$mixtures
#'    C <- training.data$quantities
#'    loss <- evaluate_cor(X = X, Y = Y, C = C, tweak = tweak)
#' }
#'
#' start_tweak <- rep(1, nrow(X.matrix))
#' start_cor <- 1 - DTD.evCor.wrapper(start_tweak)
#' cat("Starting correlation: ", start_cor, "\n")
#'
#' catch <- descent_generalized_fista(tweak_vec = start_tweak,
#'                                    F.GRAD.FUN = DTD.grad.wrapper,
#'                                    ST.FUN = soft_thresholding,
#'                                    FACTOR.FUN = nesterov_faktor,
#'                                    EVAL.FUN = DTD.evCor.wrapper)
#'
#'
#'
#'
descent_generalized_fista <- function(tweak_vec = NA,
                                      lambda=0,
                                      maxit=100,
                                      learning.rate = NA,
                                      F.GRAD.FUN = DTD.grad.wrapper,
                                      ST.FUN = soft_thresholding,
                                      FACTOR.FUN = nesterov_faktor,
                                      EVAL.FUN = DTD.evCor.wrapper,
                                      line_search_speed = 1.01,
                                      tol = 1e-5,
                                      verbose=T){
  # if not set, set learning rate to 1.
  # The idea would be that there are intelligent ways to approximate the best learning rate. (see e.g. apgpy)
  if(is.na(learning.rate)){
    # estimating best step size using barzilai borwein initialization
    # adapted from apgpy !
    grad <- F.GRAD.FUN(tweak_vec)
    norm2 <- sqrt(sum(grad**2))
    t <- 1/norm2
    x_hat <- tweak_vec - t * grad
    g_hat <- F.GRAD.FUN(x_hat)
    learning.rate <- abs((tweak_vec - x_hat) %*% (grad - g_hat) / sum((grad - g_hat)**2))
  }

  # safety check.
  # Actually, the idea was to check if the user predefined tweak_vec vector does not fit the data. But therefore the algorithm needs to know the f-function. (E.g. dimension c in ||y-Xc|| and dimension of g of -sum(cor(c, c_g)) are not the same)
  if(any(is.na(tweak_vec))){
    stop("Tweak vector includes NAs")
  }
  lsspeed_safe <- line_search_speed
  # In each iteration the updated tweak_vec will be added to tweak.history. This matrix will be returned after all iterations.
  # tweak.history is being initialised with two replicates of the input tweak_vec, because of the nesterov step.
  tweak_old <- tweak_vec
  converge_vec <- EVAL.FUN(tweak_vec)
  y_vec <- tweak_vec
  factor <- FACTOR.FUN(2)
  cat("starting to optimize \n")
  for(iter in 2:maxit){
    grad <- F.GRAD.FUN(y_vec)

    u_vec <- ST.FUN(y_vec - learning.rate* grad, learning.rate*lambda)
    eval <- EVAL.FUN(u_vec)
    i <- 1
    switched <- TRUE
    while(converge_vec[iter - 1] < eval){ # did not decrease
      l.rate <- learning.rate/line_search_speed^i
      u_vec <- ST.FUN(y_vec - l.rate*grad, l.rate*lambda)
      summary(u_vec)
      eval.old <- eval
      eval <- EVAL.FUN(u_vec)
      if(i %% 100 == 0){
        cat("line searching step: ", i, " in ", iter, " \n")
      }
      i <- i + 1
      if(abs(eval.old - eval) < tol && switched && i > 100){ # means that decreasing learning rate does not change the result
        line_search_speed <- 1/line_search_speed
        i <- 2
        cat("switched\n")
        switched <- FALSE
      }
      if( l.rate > 1e50){ # honestly, here i don't know what to do ...
        cat("Can't increase Loss anymore\n")
        ret <- list("Tweak"=tweak_vec, "Convergence"=converge_vec)
        return(ret)
      }
    }
    line_search_speed <- lsspeed_safe
    converge_vec <- c(converge_vec, eval)
    tweak_old <- tweak_vec

    tweak_vec <- u_vec

    v_vec <- tweak_old + 1/factor * (tweak_vec - tweak_old)

    # now factor ---> iter + 1
    factor <- FACTOR.FUN(iter)
    y_vec <- (1-factor) * tweak_vec + factor * v_vec

    # user can set if the algorithm prints out information in each step
    if(verbose){
      cat("###############################################\n")
      cat("iter: ", iter, "\n")
      cat("Loss: ", rev(converge_vec)[1], "\n")
      cat("factor: ", factor, "\n")
      if(iter %% 100 == 0){
        plot(1:iter, converge_vec)
      }
    }
  }
  ret <- list("Tweak"=tweak_vec, "Convergence"=converge_vec)
  return(ret)
}
