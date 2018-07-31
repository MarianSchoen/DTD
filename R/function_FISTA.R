#' Descent generalized FISTA \cr
#' (fast iterative shrinkage thresholding algorithm)
#'
#' descent_generalized_fista takes as input a vector, a gradient function and an
#' evaluation function (and some additional parameters/functions). Then, it
#' iteratively minimizes the tweak vector via the FISTA algorithms. Basically,
#' the following equations are used:
#'
#' The gradient and evaluation function both take only one argument, the soft thresholding function two.
#' If your gradient, evaluation or soft thresholding function require more arguments, please write a wrapper.
#' Exemplary wrapper functions can be found in the examples.
#'
#' @param tweak_vec numeric vector, with which the minimization algorithm starts
#' @param lambda float, regularization factor for ST.FUN function
#' @param maxit integer, maximum number of iterations for the iterative
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
#' @param EVAL.FUN function, with one parameter: vector with same length as
#'   tweak_vec, and a single float as return. This function evaluates the loss
#'   function.
#' @param line_search_speed numeric, factor with which the learning rate changes,
#'  if the optimium has not been found
#' @param verbose boolean, should information be printed to console
#' @param cycles integer, in each iteration one gradient is calculated. To find the
#'  best step size, with "cycles" different step sizes
#' @param save_all_tweaks boolean, should all tweak vectores during all iterations be stored
#'
#' @export
#'
#' @return list, including [["Converge"]] and [["Tweak"]]
#'
#' @examples
#' library(DTD)
#' random.data <- generate.random.data(nTypes = 5,
#'                                     nSamples.perType = 20,
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
#'                          size = ceiling(length(all.of.type) * percentage.of.all.cells),
#'                          replace = FALSE)
#'   samples.to.remove <- c(samples.to.remove, chosen.for.X)
#'
#'   average <- rowSums(random.data[, chosen.for.X, drop = F])
#'   X.matrix[, l.type] <- average
#' }
#'
#' # here, I declare "Type1" as Tumor cells, and all other as immune cells
#' special.samples <- c("Type1")
#' all.samples <- unique(indicator.list)
#' sample.names <- all.samples[- which(all.samples %in% special.samples)]
#'
#' # all samples that have been used in the reference matrix, must not be included in the test/training set
#' reduced.mat <- random.data[, -which(colnames(random.data) %in% samples.to.remove)]
#' # all remaining samples will be equally split into test and training:
#' test.samples <- sample(x = colnames(reduced.mat), size = ceiling(ncol(reduced.mat)/2), replace = FALSE)
#' test.mat <- random.data[, test.samples]
#' indicator.list.test <- indicator.list[-which(!names(indicator.list) %in% colnames(test.mat))]
#'
#' train.samples <- colnames(reduced.mat)[which(!colnames(reduced.mat) %in% test.samples)]
#' train.mat <- random.data[, train.samples]
#' indicator.list.train <- indicator.list[-which(!names(indicator.list) %in% colnames(train.mat))]
#'
#'
#' training.data <- mix.samples.jitter(sample.names = sample.names,
#'                                      special.samples = special.samples,
#'                                      nSamples = 1e3,
#'                                      datamatrix = train.mat,
#'                                      pheno = indicator.list.train,
#'                                      singleSpecial = F,
#'                                      add_jitter = T,
#'                                      chosen.mean = 1,
#'                                      chosen.sd = 0.05,
#'                                      min.amount.samples = 1,
#'                                      verbose = FALSE,
#'                                      included.in.X = include.in.X)
#'
#' test.data <-  mix.samples.jitter(sample.names = sample.names,
#'                                      special.samples = special.samples,
#'                                      nSamples = 1e3,
#'                                      datamatrix = test.mat,
#'                                      pheno = indicator.list.test,
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
#'                                    EVAL.FUN = DTD.evCor.wrapper,
#'                                    line_search_speed = 2,
#'                                    maxit = 1e3,
#'                                    save_all_tweaks = T)
#'
#' print(ggplot_correlation(fista.output = catch,
#'                          test.set = test.data,
#'                          X.matrix = X.matrix))
#'
#'
descent_generalized_fista <- function(tweak_vec = NA,
                                      lambda=0,
                                      maxit=1e3,
                                      learning.rate = NA,
                                      F.GRAD.FUN = DTD.grad.wrapper,
                                      ST.FUN = soft_thresholding,
                                      FACTOR.FUN = nesterov_faktor,
                                      EVAL.FUN = DTD.evCor.wrapper,
                                      line_search_speed = 2,
                                      cycles=50,
                                      save_all_tweaks=F,
                                      verbose=T){


  if(is.na(learning.rate)){
    # estimating best step size using barzilai borwein initialization
    # adapted from apgpy !
    grad <- F.GRAD.FUN(tweak_vec)
    norm2 <- sqrt(sum(grad**2))
    t <- 1/norm2
    x_hat <- tweak_vec - t * grad
    g_hat <- F.GRAD.FUN(x_hat)
    learning.rate <- abs((tweak_vec - x_hat) %*% (grad - g_hat) / sum((grad - g_hat)**2))
    # learning rate is a matrix => make it a numeric:
    learning.rate <- as.numeric(learning.rate)
  }



  # safety check.
  if(any(is.na(tweak_vec))){
    stop("Tweak vector includes NAs")
  }


  tweak_old <- tweak_vec
  converge_vec <- EVAL.FUN(tweak_vec)
  y_vec <- tweak_vec
  factor <- FACTOR.FUN(2)
  if(verbose){
    cat("starting to optimize \n")
  }

  if(save_all_tweaks){
    tweak.history <- matrix(NA, nrow = length(tweak_vec), ncol = maxit)
    tweak.history[, 1] <- tweak_vec
  }

  # step size sequence:
  sequence <- seq(from = 0, to = learning.rate, length.out = cycles)


  for(iter in 2:maxit){
    grad <- F.GRAD.FUN(y_vec)


    # u_mat will be a matrix with new u_vecs in each row
    u_mat <- matrix(unlist(lapply(sequence,
                           function(x){ST.FUN(y_vec - x * grad, x*lambda)}),
                    use.names = FALSE), nrow = cycles, byrow = T )

    # use evaluation function on each row => eval.vec is a vector with Loss-function values
    eval.vec <- apply(u_mat, 1, EVAL.FUN)

    winner.pos <- which.min(eval.vec)

    if(winner.pos == cycles){
      learning.rate <- learning.rate * line_search_speed
      sequence <- seq(from = 0, to = learning.rate, length.out = cycles)
    }
    if(winner.pos == 1){
      learning.rate <- learning.rate / line_search_speed
      sequence <- seq(from = 0, to = learning.rate, length.out = cycles)
    }

    u_vec <- u_mat[winner.pos, ]
    eval <- eval.vec[winner.pos]

    converge_vec <- c(converge_vec, eval)
    tweak_old <- tweak_vec

    tweak_vec <- u_vec

    if(save_all_tweaks){
      tweak.history[, iter] <- tweak_vec
    }

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
  if(save_all_tweaks){
    ret$History <- tweak.history
  }
  return(ret)
}
