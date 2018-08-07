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
#'
#' @return list with two entries. "quantities" matrix (nrow = ncol(datamatrix), ncol = nMixtures) and "mixture"
#' matrix (nrow = nrow(datamatrix), ncol = nMixtures)
#'
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
#' training.data <- random.data[ , -which(colnames(random.data) %in% samples.to.remove)]
#'
#' training.data <- mix.samples.jitter(sample.names = sample.names,
#'                                      special.samples = special.samples,
#'                                      nSamples = 1e3,
#'                                      datamatrix = random.data,
#'                                      pheno = indicator.list,
#'                                      singleSpecial = F,
#'                                      add_jitter = T,
#'                                      chosen.mean = 1,
#'                                      chosen.sd = 0.05,
#'                                      min.amount.samples = 1)
mix.samples.jitter <- function(sample.names,
                               special.samples,
                               nSamples=1e3,
                               datamatrix,
                               pheno,
                               verbose=F,
                               singleSpecial = F,
                               add_jitter = F,
                               chosen.mean=1,
                               chosen.sd=0.05,
                               min.amount.samples = 1,
                               included.in.X){
  # Safety check
  if(any(!names(pheno) %in% colnames(datamatrix))){
    stop("not all names provided are within the data matrix")
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
    quant.special  <- runif(n = 1, min = 0.5, max = 0.7)
    if(singleSpecial){
      special <- c(quant.special, rep(0, length(special.samples)-1))
    }else{
      remaining <- quant.special
      count <- 1
      special <- c()
      while(count <= length(special.samples)){
        tmp.special <- runif(n=1, min=0, max=remaining)
        special <- c(special, tmp.special)
        remaining <- remaining - tmp.special
        count <- count + 1
      }
      # fill last one, such that sum(tumors) equals quant.tumors
      special[count-1] <- special[count-1] + quant.special - sum(special)
    }
    special <- sample(special) # shuffle it, so every special sample can be the biggest one
    names(special) <- special.samples

    ### mixture quantities:
    remaining <- 1 - quant.special
    count <- 1
    other <- c()
    while(count <= length(sample.names)){
      tmp <- runif(n=1, min=0, max=remaining)
      other <- c(other, tmp)
      remaining <- remaining - tmp
      count <- count + 1
    }
    other[count-1] <- other[count-1] + (1 - quant.special) - sum(other)
    other <- sample(other)
    names(other) <- sample.names

    quant.matrix <- cbind(quant.matrix, c(other, special))
    colnames(quant.matrix)[lrun] <- paste0("mix", lrun)

    ### mixture: initialise with zero for each gene
    mixture <- rep(0, nrow(datamatrix))
    for(l1 in c(sample.names, special.samples)){
      #### get one of the samples
      potentialCells <- names(pheno[pheno == l1])
      if(length(potentialCells) == 0){# this cell type can not be found in the mixture ...
        # i am totally unsure what to do now ...
        next_quant <- quant.matrix[l1, lrun] *
          apply(datamatrix, 1, mean) *
          abs(rnorm(nrow(datamatrix), mean=chosen.mean, sd=chosen.sd))
      }
      if(length(potentialCells) < min.amount.samples){ # only a view cells are within the pool => average over all of them
        next_quant <- rowSums(datamatrix[, potentialCells, drop = FALSE])
      }else{ # there are enough cells in here, such that sampling is possible
        chosen.sample <-  sample(x = names(pheno[pheno == l1]), size = 1)
        next_quant <- (quant.matrix[l1, lrun] * datamatrix[, chosen.sample])
      }
      if(add_jitter){
        factors <- rnorm(n=length(next_quant), mean=chosen.mean, sd=chosen.sd)
        next_quant <- next_quant * factors
      }
      mixture <- mixture + next_quant
    }
    mix.matrix <- cbind(mix.matrix, mixture)
    colnames(mix.matrix)[lrun] <- paste0("mix", lrun)
  }

  # only keep the quantity information of cells that are included in X
  quant.matrix <- quant.matrix[included.in.X, ]

  rets <- vector(mode="list")
  rets[["mixtures"]] <- mix.matrix
  rets[["quantities"]] <- quant.matrix
  return(rets)
}
