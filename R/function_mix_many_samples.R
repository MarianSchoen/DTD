#' Mix samples for loss-function learning DTD
#'
#' mix_samples takes a gene expresssion matrix, and pheno information.
#' It then mixes the samples with known quantities such that it can be
#' used for loss-function learning digital tissue deconvolution.
#' For mixture it randomly selects "n.samples" samples from "gene.mat", and averages over them.
#' Using the information stored in pheno, it can get the quantities per cell in each mixture.
#'
#' @param gene.mat numeric matrix, with features as rows and samples as columns
#' @param pheno named vector, with pheno information for each sample in gene.mat
#' @param n.samples integer, numbers of samples to be drawn (defaults to 1000)
#' @param n.per.mixture integer, how many samples should be included per mixutre. (Default: 100)
#' @param included.in.X named vector of strings, indicating types that are in the reference matrix.
#' Only those types, and sorted in that order, will be included in the quantity matrix
#' @param verbose boolean, should information be printed (Default: FALSE)
#'
#' @return list with random profiles, and their quantity matrix
#'
#' @export
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
mix_samples <- function(gene.mat,
                       pheno,
                       included.in.X,
                       n.samples = 1e3,
                       n.per.mixture = 100,
                       verbose = FALSE){

  if(any(!names(pheno) %in% colnames(gene.mat))){
    stop("pheno does not match colnames of gene mat")
  }

  # initialise geneExpression matrix (here, the mixtures will be stored)
  geneExpression <- matrix(NA, nrow=nrow(gene.mat), ncol=n.samples)
  rownames(geneExpression) <- rownames(gene.mat)
  colnames(geneExpression) <- paste0("mixtures", 1:n.samples)

  # which types are within the pheno?
  types <- unique(pheno)
  # how many types?
  nTypes <- length(types)

  # initialise quantities matrix (here, the quantities for each cell type in every mixture will be stored)
  quantities <- matrix(NA, nrow = nTypes, ncol=n.samples)
  rownames(quantities) <- types
  colnames(quantities) <- colnames(geneExpression)

  for(lsample in 1:n.samples){
    # randomly select 'n.per.mixture' samples (without replacing!)
    randomSamples <- sample(x = colnames(gene.mat),
                            size = n.per.mixture,
                            replace = FALSE)
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
    quants <- table.chosen.pheno/n.per.mixture
    # order quants by the rownames of quantities ...
    quants <- quants[rownames(quantities)]
    # ... and add them to the quantities matrix
    quantities[, lsample] <- quants
    if(verbose){
      cat("done ", 100 * lsample/n.samples, " %\n")
    }
  }

  # only the information about  cells "included in X" should be stored in quantities:
  quantities <- quantities[included.in.X, ]

  # normalize the expression matrix
  geneExpression <- normalize_to_count(geneExpression)
  # and return both matrices as a list:
  ret <- list("mixtures"= geneExpression, "quantities" = quantities)
  return(ret)
}

