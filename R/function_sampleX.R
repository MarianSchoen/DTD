#' Sample random X
#'
#' The loss-function learning digital tissue deconvolution approach published by Goertler et al 2018
#' estimates cell compositions for a given reference matrix X (supervised deconvolution).
#' Basically, there are two methods to specify the reference profiles in X.
#' Either they are selected using external knowledge (e.g. additional measurements) or they are
#' randomly selected out of the complete data set.
#' The sample_random_X function is an implementation for the second method.
#'
#' For each entry of 'included.in.X', 'percentage.of.all.cells' are randomly selected.
#' Then, the reference profile is built by adding up all selected profiles of a type.
#' Afterwards, the reference profiles are normalized to a total number of counts.
#'
#' For examples see the DTD vignette: browseVignettes("DTD")
#'
#' @param included.in.X vector of strings, which cell types should be included in X?
#' @param pheno named vector of strings, names have to match 'colnames(expr.data)'.
#' Information about the cell type (values of vector) for each sample (names of vector)
#' @param expr.data numeric matrix with features as rows, and samples as columns
#' @param percentage.of.all.cells 0 < float < 1, which percentage of all possible cells should
#'  be use to generate a cell type profile? Defaults to 0.1
#' @param normalize.to.count logical, normalize the reference profiles? Defaults to TRUE
#'
#' @return list with two entries: (1) X.matrix: numeric matrix with as many rows as expr.data,
#' and as many columns as length(included.in.X)
#' (2) samples.to.remove: vector of strings, all samples that have been used in generating X.
#' @export
sample_random_X <- function(included.in.X,
                            pheno,
                            expr.data,
                            percentage.of.all.cells = 0.1,
                            normalize.to.count = TRUE) {

  # Safety checks:
  if(any(is.numeric(percentage.of.all.cells)) && length(percentage.of.all.cells) == 1){
    if(percentage.of.all.cells <= 0 || percentage.of.all.cells >= 1){
      stop("in sample_random_X: 'percentage.of.all.cells' must be above 0 and below 1")
    }
  }else{
    stop("in sample_random_X: 'percentage.of.all.cells' is not a single numeric")
  }

  if(!any(included.in.X %in% pheno)){
    stop("in sample_random_X: no cell type in 'included.in.X' fits 'pheno'")
  }

  if(!(all(names(pheno) %in% colnames(expr.data)) && length(pheno) == ncol(expr.data))){
    stop("in sample_random_X: 'names(pheno)' do not fit 'colnames(expr.data)'. For every entry of 'colnames(expr.data)' there has to be an entry in 'pheno'")
  }

  if(!is.matrix(expr.data)){
    stop("in sample_random_X: 'expr.data' is no matrix")
  }

  # test: normalize.to.count:
  test <- test_logical(
    test.value = normalize.to.count,
    output.info = c("sample_random_X", "normalize.to.count")
               )
  # end -> normalize.to.count
  ############################
  # initialise empty matrix:
  X.mat <- matrix(NA,
    nrow = nrow(expr.data),
    ncol = length(included.in.X)
  )
  colnames(X.mat) <- included.in.X
  rownames(X.mat) <- rownames(expr.data)
  # Keep track of all samples that have been used while generating X,
  # these have to be removed from the training set afterwards
  samples.to.remove <- c()

  for (l.type in included.in.X) {
    # get sample names of all cells of type "l.type"
    all.of.type <- names(pheno)[which(pheno == l.type)]

    # randomly sample some cells
    chosen.for.X <- sample(
      x = all.of.type,
      size = ceiling(length(all.of.type) * percentage.of.all.cells),
      replace = FALSE
    )

    # Add those cells which will be included in X to the list of samples.to.remove
    samples.to.remove <- c(samples.to.remove, chosen.for.X)

    # for each gene average over the selected
    average <- rowSums(expr.data[, chosen.for.X, drop = FALSE])
    X.mat[, l.type] <- average
  }
  # normalize to common number of counts:
  if(normalize.to.count){
    X.mat <- normalize_to_count(X.mat)
  }
  # list X.matrix and samples to remove
  ret <- list("X.matrix" = X.mat, "samples.to.remove" = samples.to.remove)
  return(ret)
}
