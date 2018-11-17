#' Sample random X
#'
#' The loss-function learning digital tissue deconvolution approach published by Goertler et al 2018
#' estimates cell compositions for a given reference matrix X. Basically, there are two methods to
#'  specify the reference profiles in X. Either they are selected using external knowledge
#' (e.g. additional measurements) or they are randomly selected out of the complete data set.
#' The sample.random.X function is an implementation for the second method.
#'
#' @param included.in.X vector of strings, which cell types should be included in X?
#' @param pheno named vector of strings, names have to match exp.data.
#' Information about the cell type (value of vector) for each sample (name of vector)
#' @param exp.data numeric matrix with features as rows, and samples as columns
#' @param percentage.of.all.cells float, which percentage of all possible cells should
#'  be use to generate a cell type profile?
#'
#' @return numeric matrix with as many rows as exp.data, and as many columns as length(included.in.X)
#' @export
sample.random.X <- function(included.in.X,
                            pheno,
                            exp.data,
                            percentage.of.all.cells = 0.1){


  # initialise empty matrix:
  X.mat <- matrix(NA,
                  nrow = nrow(exp.data),
                  ncol = length(included.in.X))
  colnames(X.mat) <- included.in.X
  rownames(X.mat) <- rownames(exp.data)

  samples.to.remove <- c()

  for(l.type in included.in.X){
    # get sample names of all cells of type "l.type"
    all.of.type <- names(pheno)[which(pheno == l.type)]

    # randomly sample some cells
    chosen.for.X <- sample(x = all.of.type,
                           size = ceiling(length(all.of.type) * percentage.of.all.cells),
                           replace = FALSE)

    # Add those cells which will be included in X to the list of samples.to.remove
    samples.to.remove <- c(samples.to.remove, chosen.for.X)

    # for each gene average over the selected
    average <- rowSums(exp.data[, chosen.for.X, drop = FALSE])
    X.mat[, l.type] <- average
  }

  # list X.matrix and samples to remove
  ret <- list("X.matrix" = X.mat, "samples.to.remove" = samples.to.remove)
  return(ret)
}
