#' Evaluate correlation
#'
#' The loss-function learning digital tissue deconvolution finds a vector g which minimizes the Loss function L\cr
#' \deqn{L(g) = 1 - \sum cor(true_C, estimatd_C(g))}
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
#' training.data <- mix.samples.jitter(sample.names = sample.names,
#'                                      special.samples = special.samples,
#'                                      nSamples = 1e3,
#'                                      datamatrix = random.data,
#'                                      pheno = indicator.list,
#'                                      singleSpecial = FALSE,
#'                                      add_jitter = TRUE,
#'                                      chosen.mean = 1,
#'                                      chosen.sd = 0.05,
#'                                      min.amount.samples = 1,
#'                                      verbose = FALSE,
#'                                      included.in.X = include.in.X)
#'
#' training.data$mixtures <- normalizeToCount(training.data$mixtures)
#'
#' loss <- evaluate_cor(X = X.matrix,
#'                      Y = training.data$mixtures,
#'                      C = training.data$quantities,
#'                      tweak = rep(1, nrow(X.matrix)))
#'
#' cor <- 1 - loss
#'
evaluate_cor <- function(X, Y, C, tweak){
  # estimate C using reference matrix, bulks and the tweak vector:
  esti.cs <- est.cs(X, Y, tweak)

  if(any(dim(esti.cs) != dim(C))){
    stop("evaluate_cor: dimension of estimated C do not fit")
  }

  # initialise cor_per_cType:
  cor_per_cType <- rep(NA, nrow(C))
  for(l1 in 1:nrow(C)){
    # calculate the correlation per Type:
    cor_per_cType[l1] <- stats::cor(C[l1, ], esti.cs[l1, ])
  }

  # The value of loss-function is 1 - the averaged correlation:
  loss <- 1 - mean(cor_per_cType)
  return(loss)
}
