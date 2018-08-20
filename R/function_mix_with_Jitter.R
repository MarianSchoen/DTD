#' Mix samples with Jitter
#'
#' "mix.samples.jitter" takes pheno information (sample.names, special.names) and a expression matrix.
#' Its output is a expression and a quantity matrix.
#' Each sample of the output is a mixture of input samples, multiplied with Jitter vector.
#' "mix.samples.Jitter" mixes samples in a way that they look similar to biological data.
#' For example, in a tumor tissue there are several cells included, but we expect that
#' most of them are tumor cells. Therefore, the mix.samples.Jitter function expects a
#' list of samples which represent immune cells (and occur in minor fractions) and a
#' list of special samples (which occur in major fractions)
#'
#' @param sample.names vector of strings, have to match colnames(datamatrix).
#' @param special.samples vector of strings, have to match colnames(datamatrix), these samples will
#' occur with higher quantites within the mix data.
#' @param nSamples integer, how many mixtures should be drawn
#' @param datamatrix numeric matrix, with samples as columns, and features as rows
#' @param singleSpecial boolean, should all special names be used? or only a single one
#' @param add_jitter boolean, should the mixtures be mulitplied with a vector of normally distributed numbers? (JITTER)
#' @param verbose boolean, should function tell about progression?
#' @param pheno named list of characters, indicates which of the samples in datamatrix belongs to which type in sample.names/special.names
#' @param chosen.mean numeric, mean of jitter
#' @param chosen.sd numeric, standard deviation of jitter
#' @param min.amount.samples integer, how many samples have to be present such that it averages over them,
#'  instead of taking only 1
#' @param included.in.X list of strings, which cell types are included in the reference matrix.
#' The mix function will mix all samples in the data set, but will only return the quantity matrix
#' for the types included in X (and in the right ordering)
#' @param per.type integer, how many samples "per type" should be used for each mixture
#'
#' @return list with two entries. "quantities" matrix (nrow = ncol(datamatrix), ncol = nMixtures) and "mixture"
#' matrix (nrow = nrow(datamatrix), ncol = nMixtures)
#'
#' @export
#'
#' @examples
#' library(DTD)
#' random.data <- generate.random.data(nTypes = 10,
#'                                     nSamples.perType = 10,
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
#' # here, I declare "Type1" as Tumor cells, and all other as immune cells
#' special.samples <- c("Type1")
#' all.samples <- unique(indicator.list)
#' sample.names <- all.samples[- which(all.samples %in% special.samples)]
#'
#' training.data <- random.data[ , -which(colnames(random.data) %in% samples.to.remove)]
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
#'                                      min.amount.samples = 1)
mix.samples.jitter <- function(sample.names,
                               special.samples,
                               nSamples = 1e3,
                               datamatrix,
                               pheno,
                               verbose = FALSE,
                               singleSpecial = FALSE,
                               add_jitter = FALSE,
                               chosen.mean = 1,
                               chosen.sd = 0.05,
                               min.amount.samples = 1,
                               per.type = 1,
                               included.in.X){
  # Safety check
  if(any(!names(pheno) %in% colnames(datamatrix))){
    stop("not all names provided are within the data matrix")
  }
  if(any(!c(sample.names, special.samples) %in% pheno)){
    stop("sample.names and special.samples do not match pheno")
  }
  len <- length(c(sample.names, special.samples))

  # initialise return variables:
  mix.matrix <- matrix(nrow=nrow(datamatrix), ncol=0)
  rownames(mix.matrix) <- rownames(datamatrix)

  quant.matrix <- matrix(nrow=len, ncol=0)
  rownames(quant.matrix) <- c(sample.names, special.samples)


  for(lrun in 1:nSamples){
    if(verbose){
      cat("doing ", lrun, " or ", 100*lrun/nSamples, "%\n" )
    }
    # get quantities for the special samples:
    # The special samples (e.g. malignant cells) dominate the cells (=> high quants)
    quant.special  <- stats::runif(n = 1, min = 0.5, max = 0.7)

    # singleSpecial indicates if there is only one special sample present in the mixtures:
    if(singleSpecial){
      special <- c(quant.special, rep(0, length(special.samples)-1))
    }else{
      # if there are more than one special samples, the quant.special gets split,
      # such that the sum over all special samples equals quant.special:
      remaining <- quant.special
      count <- 1
      special <- c()
      while(count <= length(special.samples)){
        tmp.special <- stats::runif(n=1, min=0, max=remaining)
        special <- c(special, tmp.special)
        remaining <- remaining - tmp.special
        count <- count + 1
      }
      # fill last one, such that sum(tumors) equals quant.tumors
      special[count-1] <- special[count-1] + quant.special - sum(special)
    }
    # shuffle it, so every special sample can be the biggest one
    special <- sample(special)
    names(special) <- special.samples

    # next, the normal samples get random quantities:
    remaining <- 1 - quant.special
    count <- 1
    other <- c()
    while(count <= length(sample.names)){
      tmp <- stats::runif(n=1, min=0, max=remaining)
      other <- c(other, tmp)
      remaining <- remaining - tmp
      count <- count + 1
    }
    # sum of quant.special + sum of other = 1:
    other[count-1] <- other[count-1] + (1 - quant.special) - sum(other)
    # sample "others" to be not sorted:
    other <- sample(other)
    names(other) <- sample.names

    # add quantities to quant.matrix:
    quant.matrix <- cbind(quant.matrix, c(other, special))
    colnames(quant.matrix)[lrun] <- paste0("mix", lrun)

    # Now, using the sampled quantities, calculated the expression profiles of the mixtures.
    # Initialise a numerc vector with zeros:
    mixture <- rep(0, nrow(datamatrix))
    # Go through all samples:
    for(lsample in c(sample.names, special.samples)){
      # Find all cells in the pheno which match to "lsample"
      potentialCells <- names(pheno[pheno == lsample])
      # if there are many samples of type "lsample" in the set:
      if(length(potentialCells) >= min.amount.samples){
        # sample a subset of the potentiallCells:
        chosen.sample <-  sample(x = potentialCells, size = per.type)
        # and multiple these with their quantity:
        next_quant <- quant.matrix[lsample, lrun] * rowSums(datamatrix[, chosen.sample, drop = FALSE])
      }else{
        # fi there are only a few samples present, we do not sample some of the.
        # The expression is the rowSum over all potentialCells, which gets multiplied with its quantity:
        next_quant <- quant.matrix[lsample, lrun] * rowSums(datamatrix[, potentialCells, drop = FALSE])
      }
      # jitter means to multply every entry of every "next_quant" with a random number close to 1:
      if(add_jitter){
        # get a vector of random numbers:
        factors <- stats::rnorm(n=length(next_quant), mean=chosen.mean, sd=chosen.sd)
        # elementwise multiply next_quant with it
        next_quant <- next_quant * factors
      }
      # add the new "next_quant" vector to the mixture:
      mixture <- mixture + next_quant
    }
    # add the mixture to the mix.matrix, and add the colname:
    mix.matrix <- cbind(mix.matrix, mixture)
    colnames(mix.matrix)[lrun] <- paste0("mix", lrun)
  }

  # only keep the quantity information of cells that are included in X
  quant.matrix <- quant.matrix[included.in.X, ]

  # build a return list, and add both matrices:
  rets <- vector(mode="list")
  rets[["mixtures"]] <- mix.matrix
  rets[["quantities"]] <- quant.matrix
  return(rets)
}

