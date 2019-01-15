#' Evaluate correlation
#'
#' The loss-function learning digital tissue deconvolution finds a vector g which minimizes the Loss function L\cr
#' \deqn{L(g) = - \sum cor(true_C, estimatd_C(g))}
#' The evaluate_cor function returns the value of the Loss function.
#' The evaluate_cor function takes 4 arguments. However, we wrote the FISTA implementation in a way that both the
#' gradient and the evaluation function only take one argument, which is the g vector. In order to use this function
#' within the `descent_generalized_fista` function, you need a wrapper function (see \code{\link{descent_generalized_fista}})
#'
#' @param X numeric matrix with cells as columns, and features as rows. Reference matrix X of the DTD problem.
#' @param Y numeric matrix with samples as columns, and features as rows. Each sample in Y is a bulk measurement,
#' for which the quantity of the cells in X are known (and saved in C)
#' @param C numeric matrix with cells as rows, and mixtures as columns.
#' Each row of C holds the distribution of the cell over all mixtures.
#' @param tweak numeric vector with length of nrow(X).
#'
#' @return float, value of the Loss function
#'
#' @export
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
evaluate_cor <- function(X, Y, C, tweak){
  # estimate C using reference matrix, bulks and the tweak vector:
  esti.cs <- est_cs(X, Y, tweak)

  if(any(dim(esti.cs) != dim(C))){
    stop("evaluate_cor: dimension of estimated C do not fit")
  }

  # initialise loss:
  loss <- 0
  for(l1 in 1:nrow(C)){
    # calculate the correlation per Type and add them up
    loss <- loss + stats::cor(C[l1, ], esti.cs[l1, ])
  }

  # Our loss function is set to be a minimization problem, therefore we invert the sign:
  loss <- -loss
  return(loss)
}

