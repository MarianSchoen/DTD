#' Mix samples for loss-function learning DTD
#'
#' mix_samples takes a gene expresssion matrix ('exp.data'),
#' and 'pheno' information.
#' It then mixes the samples with known quantities such that it can be
#' used for loss-function learning digital tissue deconvolution.
#' For a mixture it randomly selects "n.samples" samples from "exp.data", and averages over them.
#' Using the information stored in pheno, it can get the quantities per cell in each mixture.
#'
#' @param exp.data numeric matrix, with features as rows and samples as columns
#' @param pheno named vector of strings, with pheno information ('pheno') for each sample ('name(pheno)') in exp.data
#' @param n.samples integer above 0, numbers of samples to be drawn (defaults to 1000)
#' @param n.per.mixture integer above 0, below ncol(exp.data),
#'  how many samples should be included per mixutre. (Default: 100)
#' @param included.in.X vector of strings, indicating types that are in the reference matrix.
#' Only those types, and sorted in that order, will be included in the quantity matrix.
#' Notice, every profile of 'exp.data' might be included in the mixture.
#' But the quantity matrix only reports quantity information for the cell types in 'included.in.X'.
#' @param verbose logical, should information be printed? (Default: FALSE)
#' @param normalize.to.count logical, normalize each mixture? Defaults to TRUE
#'
#' @return list with two entries: 'mixtures' and 'quantities'.
#'
#' @export
#'
#' @examples
#' library(DTD)
#' random.data <- generate_random_data(
#'       n.types = 10,
#'       n.samples.per.type = 150,
#'       n.features = 250,
#'       sample.type = "Cell",
#'       feature.type = "gene"
#'       )
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
#' sample.X <- sample_random_X(
#'       included.in.X = include.in.X,
#'       pheno = indicator.list,
#'       exp.data = normalized.data,
#'       percentage.of.all.cells = percentage.of.all.cells
#'       )
#' X.matrix <- sample.X$X.matrix
#' samples.to.remove <- sample.X$samples.to.remove
#'
#' # all samples that have been used in the reference matrix, must not be included in
#' # the test/training set
#' remaining.mat <- random.data[, -which(colnames(random.data) %in% samples.to.remove)]
#' train.samples <- sample(
#'       x = colnames(remaining.mat),
#'       size = ceiling(ncol(remaining.mat)/2),
#'       replace = FALSE
#'       )
#' test.samples <- colnames(remaining.mat)[which(!colnames(remaining.mat) %in% train.samples)]
#'
#' train.mat <- remaining.mat[, train.samples]
#' test.mat <- remaining.mat[, test.samples]
#'
#' indicator.train <- indicator.list[names(indicator.list) %in% colnames(train.mat)]
#' training.data <- mix_samples(
#'       exp.data = train.mat,
#'       pheno = indicator.train,
#'       included.in.X = include.in.X,
#'       n.samples = 500,
#'       n.per.mixture = 100,
#'       verbose = FALSE
#'       )
#'
#' indicator.test <- indicator.list[names(indicator.list) %in% colnames(test.mat)]
#' test.data <-  mix_samples(
#'       exp.data = test.mat,
#'       pheno = indicator.test,
#'       included.in.X = include.in.X,
#'       n.samples = 500,
#'       n.per.mixture = 100,
#'       verbose = FALSE
#'       )
mix_samples <- function(exp.data,
                       pheno,
                       included.in.X,
                       n.samples = 1e3,
                       n.per.mixture = 100,
                       verbose = FALSE,
                       normalize.to.count = TRUE){

  # Safety checks
  if(!is.vector(included.in.X)){
    stop("in mix_samples: 'included.in.X' is not provided as vector")
  }
  if(!is.vector(pheno)){
    stop("in mix_samples: 'pheno' is not provided as vector")
  }
  if(!is.matrix(exp.data)){
    stop("in mix_samples: 'exp.data' is not a numeric matrix")
  }

  if(!any(included.in.X %in% pheno)){
    stop("in mix_samples: no cell type in 'included.in.X' fits 'pheno'")
  }

  if(!(all(names(pheno) %in% colnames(exp.data)) && length(pheno) == ncol(exp.data))){
    stop("in mix_samples: 'names(pheno)' do not fit 'colnames(exp.data)'.
         For every entry of 'colnames(exp.data)' there has to be a entry in 'pheno'")
  }
  # test: normalize.to.count:
  test <- test_logical(
    test.value = normalize.to.count,
    output.info = c("mix_samples", "normalize.to.count")
  )
  # end -> normalize.to.count

  # safety checks: n.samples
  test <- test_integer(test.value = n.samples,
                       output.info = c("mix_samples", "n.samples"),
                       min = 1,
                       max = Inf)
  # end -> n.samples

  # safety checks: n.per.mixture
  test <- test_integer(test.value = n.per.mixture,
                       output.info = c("mix_samples", "n.per.mixutre"),
                       min = 1,
                       max = ncol(exp.data))

  # end -> n.per.mixture

  # safety checks: verbose
  test <- test_logical(test.value = verbose,
                       output.info = c("mix_samples", "verbose")
                       )
  # end -> verbose

  # initialise geneExpression matrix (here, the mixtures will be stored)
  geneExpression <- matrix(NA, nrow=nrow(exp.data), ncol=n.samples)
  rownames(geneExpression) <- rownames(exp.data)
  colnames(geneExpression) <- paste0("mixtures", 1:n.samples)

  # which types are within the pheno (or 'include.in.X')?
  # if there are only a few samples, it might happen that
  # there is no sample for a cell type.
  types <- unique(pheno)
  if(!all(included.in.X %in% types)){
    useable.types <- unique(
      c(
        pheno,
        included.in.X
      )
    )
    stop(
      paste(
        "In 'mix_samples': for some types in 'included.in.X' there are no samples in your data. Provide e.g. the following vector as 'included.in.X': \n",
        paste0(
          useable.types
          , collapse = ", "
        )
      )
    )
  }

  # how many types?
  nTypes <- length(types)

  # initialise quantities matrix (here, the quantities for each cell type in every mixture will be stored)
  quantities <- matrix(NA, nrow = nTypes, ncol=n.samples)
  rownames(quantities) <- types
  colnames(quantities) <- colnames(geneExpression)

  for(lsample in 1:n.samples){
    # randomly select 'n.per.mixture' samples (without replacing!)
    randomSamples <- sample(x = colnames(exp.data),
                            size = n.per.mixture,
                            replace = FALSE)
    # average over the selected samples ...
    avg <- rowSums(exp.data[, randomSamples])
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
  if(normalize.to.count){
    geneExpression <- normalize_to_count(geneExpression)
  }
  # and return both matrices as a list:
  ret <- list("mixtures"= geneExpression, "quantities" = quantities)
  return(ret)
}

