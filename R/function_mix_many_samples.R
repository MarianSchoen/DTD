#' Mix samples for loss-function learning DTD
#'
#' mix.samples takes a gene expresssion matrix, and pheno information.
#' It then mixes the samples with known quantities such that it can be
#' used for loss-function learning digital tissue deconvolution.
#'
#' @param gene.mat numeric matrix, with features as rows and samples as columns
#' @param pheno matrix, with pheno information for each sample in gene.mat
#' @param nSamples integer, numbers of samples to be drawn
#' @param nPerMixture integer, how many samples should be included per mixutre
#' @param included.in.X list of strings, types that are in the reference matrix.
#' Only those, and sorted like that, must be included in the quantity matrix
#' @param verbose boolean, should information be printed
#'
#' @return list with random profiles, and their quantity matrix
#'
#' @export
#'
#' @examples
#' library(DTD)
#' random.data <- generate.random.data(nTypes = 10,
#'                                     nSamples.perType = 1e2,
#'                                     nFeatures = 500,
#'                                     sample.type = "Cell",
#'                                     feature.type = "gene")
#'
#' # normalize all samples to the same amount of counts:
#' random.data <- normalizeToCount(random.data)
#'
#' # extract indicator list.
#' # This list contains the type of the sample as value, and the sample name as names
#' indicator.list <- gsub("^Cell[0-9]*\\.", "", colnames(random.data))
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
#' train.mat <- random.data[, -which(colnames(random.data) %in% samples.to.remove)]
#' indicator.train <- indicator.list[names(indicator.list) %in% colnames(train.mat)]
#'
#' training.data <- mix.samples(gene.mat = train.mat,
#'                              pheno = indicator.train,
#'                              included.in.X = include.in.X)
#'
mix.samples <- function(gene.mat,
                       pheno,
                       included.in.X,
                       nSamples = 1e3,
                       nPerMixture = 100,
                       verbose = FALSE){

  if(any(!names(pheno) %in% colnames(gene.mat))){
    stop("pheno does not match colnames of gene mat")
  }

  # initialise geneExpression matrix (here, the mixtures will be stored)
  geneExpression <- matrix(NA, nrow=nrow(gene.mat), ncol=nSamples)
  rownames(geneExpression) <- rownames(gene.mat)
  colnames(geneExpression) <- paste0("mixtures", 1:nSamples)

  # which types are within the pheno?
  types <- unique(pheno)
  # how many types?
  nTypes <- length(types)

  # initialise quantities matrix (here, the quantities for each cell type in every mixture will be stored)
  quantities <- matrix(NA, nrow = nTypes, ncol=nSamples)
  rownames(quantities) <- types
  colnames(quantities) <- colnames(geneExpression)

  for(lsample in 1:nSamples){
    # randomly select 'nPerMixture' samples (without replacing!)
    randomSamples <- sample(x = colnames(gene.mat), size = nPerMixture, replace = FALSE)
    # average over the selected samples ...
    avg <- rowSums(gene.mat[, randomSamples])
    # ... store the result in the geneExpression matrix:
    geneExpression[,lsample] <- avg

    # Next, we extract how often each cell type has been selected in this mixture:
    chosen.pheno <- pheno[which(names(pheno) %in% randomSamples)]
    table.chosen.pheno <- table(chosen.pheno)

    # due to the randomness some types may not be included in this sample:
    missing.types <- rownames(quantities)[which(!rownames(quantities) %in% names(table.chosen.pheno))]
    # if there are any missing types, add those with quantity 0
    if(length(missing.types) != 0){
      add_0 <- rep(0, length(missing.types))
      names(add_0) <- missing.types
      table.chosen.pheno <- c(table.chosen.pheno, add_0)
    }

    # divide quants by the number of used samples:
    quants <- table.chosen.pheno/nPerMixture
    # order quants by the rownames of quantities ...
    quants <- quants[rownames(quantities)]
    # ... and add them to the quantities matrix
    quantities[, lsample] <- quants
    if(verbose){
      cat("done ", 100 * lsample/nSamples, " %\n")
    }
  }

  # only the information about  cells "included in X" should be stored in quantities:
  quantities <- quantities[included.in.X, ]

  # normalize the expression matrix
  geneExpression <- normalizeToCount(geneExpression)
  # and return both matrices as a list:
  ret <- list("mixtures"= geneExpression, "quantities" = quantities)
  return(ret)
}

