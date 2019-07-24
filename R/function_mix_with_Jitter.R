#' Mix samples with Jitter
#'
#' "mix_samples_with_jitter" takes pheno information ('minor.fraction.type', 'special.names')
#' and a expression matrix ('exp.data').
#' Its output is a list of two matrices. A matrix holding in-silicio 'mixtures' and a quantity matrix,
#' holding the corresponding compositions.
#' Each profile of the 'mixture' matrix is the result of a linear combination of input profiles ('exp.data'),
#' multiplied with a random jitter vector.
#' "mix_samples_with_jitter" mixes samples in a way that they look similar to biological data.
#' For example, in a tumor tissue there are several cells included, but we expect that
#' most of them are tumor cells. Therefore, the mix_samples_with_jitter function expects a
#' vector of samples which represent cells, that occur in minor fractions (e.g. immune cells) and a
#' vector of special samples, which occur in major fractions.
#'
#' @param minor.fraction.type vector of strings, have to match names(pheno).
#' @param major.fraction.type vector of strings, have to match names(pheno), these samples will
#' occur with higher quantites in the in-silicio mixtures.
#' @param n.samples integer above 0, numbers of samples to be drawn (defaults to 1000)
#' @param exp.data non-negative numeric matrix, with features as rows and samples as columns
#' @param single.special logical, in a mixture, should all special names be used? (=> TRUE) Or, should a mixture only include one special sample? (=> FALSE)(Defaults to FALSE)
#' @param add.jitter logical, should the mixtures be multiplied with a vector of normally distributed numbers? (JITTER)
#' @param verbose logical, should the function print about progress to the terminal? (Defaults to FALSE)
#' @param pheno named vector of strings, with pheno information ('pheno') for each sample ('names(pheno)') in exp.data
#' @param chosen.mean float, mean of jitter (Default: 1)
#' @param chosen.sd float, standard deviation of jitter (Default: 0.05)
#' @param included.in.X vector of strings, indicating types that are in the reference matrix.
#' Only those types, and sorted in that order, will be included in the quantity matrix.
#' Notice, every profile of 'exp.data' might be included in the mixture.
#' But the quantity matrix only reports quantity information for the cell types in 'included.in.X'.
#' @param n.per.mixture integer, 1 <= 'n.per.mixture', how many samples per type should be used for each mixture (Default: 1)
#' @param min.major float, 0 <= 'min.major' <= 'max.major', minimal fraction of 'major.fraction.type' in each mixture
#' @param max.major float, 'min.major' <= 'max.major' <= 1, maximal fraciton of 'major.fraction.type' in each mixture
#'
#' @return list with two entries. "quantities" matrix (nrow = ncol(exp.data), ncol = nMixtures) and "mixture"
#' matrix (nrow = nrow(exp.data), ncol = nMixtures)
#'
#' @export
#'
#' @examples
#' library(DTD)
#' random.data <- generate_random_data(n.types = 10,
#'                                     n.samples.per.type = 10,
#'                                     n.features = 500,
#'                                     sample.type = "Cell",
#'                                     feature.type = "gene")
#'
#' # normalize all samples to the same amount of counts:
#' random.data <- normalize_to_count(random.data)
#'
#' # extract indicator list.
#' # This list contains the type of the sample as value, and the sample name as names
#' indicator.list <- gsub("^Cell[0-9]*\\.", "", colnames(random.data))
#' names(indicator.list) <- colnames(random.data)
#'
#' # First, decide which cells should be deconvoluted.
#' include.in.X <- c("Type2", "Type3", "Type4", "Type5")
#'
#' # here, I declare "Type1" as Tumor cells, and all other as immune cells
#' major.fraction.type <- c("Type1")
#' all.samples <- unique(indicator.list)
#' minor.fraction.type <- all.samples[- which(all.samples %in% major.fraction.type)]
#'
#' training.data <- mix_samples_with_jitter(
#'     minor.fraction.type = minor.fraction.type,
#'     major.fraction.type = major.fraction.type,
#'     included.in.X = include.in.X,
#'     n.samples = 1e3,
#'     exp.data = random.data,
#'     pheno = indicator.list,
#'     single.special = FALSE,
#'     add.jitter = TRUE,
#'     chosen.mean = 1,
#'     chosen.sd = 0.05,
#'     min.major = 0.5,
#'     max.major = 0.5
#' )
#'
mix_samples_with_jitter <- function(minor.fraction.type, # tested
                               major.fraction.type, # tested
                               n.samples, # tested
                               exp.data, # tested
                               pheno, # tested
                               verbose = FALSE, # tested
                               single.special = FALSE, # tested
                               add.jitter = FALSE, # tested
                               chosen.mean = 1, # tested
                               chosen.sd = 0.05, # tested
                               n.per.mixture = 1, # tested
                               included.in.X, # tested
                               min.major = 0.5, # tested
                               max.major = 0.7
                               ){

  # Safety check: pheno, exp.data
  if(!is.vector(included.in.X)){
    stop("in mix_samples_with_jitter: 'included.in.X' is not provided as vector")
  }
  if(!is.vector(pheno)){
    stop("in mix_samples_with_jitter: 'pheno' is not provided as vector")
  }
  if(!any(included.in.X %in% pheno)){
    stop("in mix_samples_with_jitter: no cell type in 'included.in.X' fits 'pheno'")
  }

  if(!(is.matrix(exp.data) && is.numeric(exp.data))){
    stop("in mix_samples_with_jitter: 'exp.data' is not a numeric matrix")
  }

  if(!(all(names(pheno) %in% colnames(exp.data)) && length(pheno) == ncol(exp.data))){
    stop("in mix_samples_with_jitter: 'names(pheno)' do not fit 'colnames(exp.data)'.
         For every entry of 'colnames(exp.data)' there has to be a entry in 'pheno'")
  }


  # end -> pheno, exp.data

  # Safety check: major.fraction.type
  if(!is.vector(major.fraction.type)){
    stop("in mix_samples_with_jitter: 'major.fraction.type' is not provided as vector")
  }
  if(!all(major.fraction.type %in% pheno)){
    stop("in mix_samples_with_jitter: There is type in 'major.fraction.type' that can not be matched to 'pheno'")
  }
  # end -> major.fraction.type
  # Safety check: minor.fraciton.type
  if(!is.vector(minor.fraction.type)){
    stop("in mix_samples_with_jitter: 'minor.fraction.type' is not provided as vector")
  }
  if(!all(minor.fraction.type %in% pheno)){
    stop("in mix_samples_with_jitter: There is type in 'minor.fraciton.type' that can not be matched to 'pheno'")
  }
  # end -> minor.fraciton.type

  len <- length(c(minor.fraction.type, major.fraction.type))
  # Safety check: verbose
  test <- test_logical(test.value = verbose,
                       output.info = c("mix_samples_with_jitter", "verbose"))
  # end -> verbose
  # Safety check: single.special
  test <- test_logical(test.value = single.special,
                       output.info = c("mix_samples_with_jitter", "single.special"))
  # end -> single.special
  # Safety check: add.jitter
  test <- test_logical(test.value = add.jitter,
                       output.info = c("mix_samples_with_jitter", "add.jitter"))
  # end -> add.jitter
  # Safety check: chosen.mean
  test <- test_numeric(test.value = chosen.mean,
                       output.info = c("mix_samples_with_jitter", "chosen.mean"),
                       min = -Inf,
                       max = Inf)
  # end -> chosen.mean
  # Safety check: chosen.sd
  test <- test_numeric(test.value = chosen.sd,
                       output.info = c("mix_samples_with_jitter", "chosen.sd"),
                       min = -Inf,
                       max = Inf)
  # end -> chosen.sd
  # Safety check: n.per.mixture
  test <- test_integer(test.value = n.per.mixture,
                       output.info = c("mix_samples_with_jitter", "n.per.mixture"),
                       min = 0,
                       max = Inf)
  # end -> chosen.sd

  # Safety check: min.major
  test <- test_numeric(test.value = min.major,
                       output.info = c("mix_samples_with_jitter", "min.major"),
                       min = 0,
                       max = max.major)
  # end -> chosen.sd
  # Safety check: chosen.sd
  test <- test_numeric(test.value = max.major,
                       output.info = c("mix_samples_with_jitter", "max.major"),
                       min = min.major,
                       max = 1)
  # end -> chosen.sd

  # Safety check: n.samples
  test <- test_integer(test.value = n.samples,
                       output.info = c("mix_samples_with_jitter", "n.samples"),
                       min = 1,
                       max = Inf)
  # end -> n.samples

  # initialise return variables:
  mix.matrix <- matrix(nrow=nrow(exp.data), ncol=0)
  rownames(mix.matrix) <- rownames(exp.data)

  quant.matrix <- matrix(nrow=len, ncol=0)
  rownames(quant.matrix) <- c(minor.fraction.type, major.fraction.type)


  for(lrun in 1:n.samples){
    if(verbose){
      cat("doing ", lrun, " or ", 100*lrun/n.samples, "%\n" )
    }
    # get quantities for the special samples:
    # The special samples (e.g. malignant cells) dominate the cells (=> high quants)
    quant.special  <- stats::runif(n = 1,
                                   min = min.major,
                                   max = max.major)

    # single.special indicates if there is only one special sample present in the mixtures:
    if(single.special){
      special <- c(quant.special, rep(0, length(major.fraction.type)-1))
    }else{
      # if there are more than one special samples, the quant.special gets split,
      # such that the sum over all special samples equals quant.special:
      remaining <- quant.special
      counter <- 1
      special <- c()
      while(counter <= length(major.fraction.type)){
        tmp.special <- stats::runif(n=1, min=0, max=remaining)
        special <- c(special, tmp.special)
        remaining <- remaining - tmp.special
        counter <- counter + 1
      }
      # fill last one, such that sum(tumors) equals quant.tumors
      special[counter-1] <- special[counter-1] + quant.special - sum(special)
    }
    # shuffle it, so every special sample can be the biggest one
    special <- sample(special)
    names(special) <- major.fraction.type

    # next, the normal samples get random quantities:
    remaining <- 1 - quant.special
    counter <- 1
    other <- c()
    while(counter <= length(minor.fraction.type)){
      tmp <- stats::runif(n=1, min=0, max=remaining)
      other <- c(other, tmp)
      remaining <- remaining - tmp
      counter <- counter + 1
    }
    # sum of quant.special + sum of other = 1:
    other[counter-1] <- other[counter-1] + (1 - quant.special) - sum(other)
    # sample "others" to be not sorted:
    other <- sample(other)
    names(other) <- minor.fraction.type

    # add quantities to quant.matrix:
    quant.matrix <- cbind(quant.matrix, c(other, special))
    colnames(quant.matrix)[lrun] <- paste0("mix", lrun)

    # Now, using the sampled quantities, calculated the expression profiles of the mixtures.
    # Initialise a numerc vector with zeros:
    mixture <- rep(0, nrow(exp.data))
    # Go through all samples:
    for(lsample in c(minor.fraction.type, major.fraction.type)){
      # Find all cells in the pheno which match to "lsample"
      potentialCells <- names(pheno[pheno == lsample])
      # sample a subset of the potentiallCells:
      chosen.sample <- sample(x = potentialCells,
                              size = n.per.mixture,
                              replace = TRUE)
      # and multiple these with their quantity:
      next_quant <- quant.matrix[lsample, lrun] * rowSums(exp.data[, chosen.sample, drop = FALSE])
      # jitter means to multply every entry of every "next_quant" with a random number close to 1:
      if(add.jitter){
        # get a vector of random numbers:
        factors <- stats::rnorm(
          n=length(next_quant),
          mean=chosen.mean,
          sd=chosen.sd)
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
  if(!any(is.na(included.in.X))){
    quant.matrix <- quant.matrix[included.in.X, ]
  }
  # build a return list, and add both matrices:
  rets <- vector(mode="list")
  rets[["mixtures"]] <- mix.matrix
  rets[["quantities"]] <- quant.matrix
  return(rets)
}

