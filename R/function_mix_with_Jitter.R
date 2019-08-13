#' Mix samples with jitter
#'
#' "mix_samples_with_jitter" takes 'pheno' information and a expression matrix ('exp.data').
#' Its output is a list of two matrices. A matrix holding in-silicio 'mixtures' and a quantity matrix,
#' holding the corresponding compositions.
#' Each profile of the 'mixture' matrix is the result of a linear combination of input expression,
#' multiplied with a random jitter vector.
#'
#' In the DTD package there are 2 ways to generate 'in-silicio' mixtures:
#'  - take 'n.per.mixture' random profiles from the expression matrix, and average over them.
#'  Here, the quantities (or probabilities) for each cell type is set by their relative
#'  frequency in the expression matrix. (DTD::mix_samples)
#'  - for each cell type, pick a random number, interpret this number as the cell types quantity.
#'  Multiply the profile with its quantity, and average over the resulting mixture. (DTD::mix_samples_with_jitter)
#'
#' @param n.samples integer above 0, numbers of samples to be drawn (defaults to 1000)
#' @param prob.each numeric vector with same length as 'include.in.X.' For each cell type in 'include.in.X'
#' 'prob.each' holds the average quantity in the mixtures.
#' @param exp.data non-negative numeric matrix, with features as rows and samples as columns
#' @param add.jitter logical, should the mixtures be multiplied with a vector of normally distributed numbers? (JITTER)
#' @param verbose logical, should the function print about progress to the terminal? (Defaults to FALSE)
#' @param pheno named vector of strings, with pheno information ('pheno') for each sample ('names(pheno)') in exp.data
#' @param chosen.mean float, mean of jitter (Default: 1)
#' @param chosen.sd float, standard deviation of jitter (Default: 0.05)
#' @param include.in.X vector of strings, indicating types that are in the reference matrix.
#' Only those types, and sorted in that order, will be included in the quantity matrix.
#' Notice, every profile of 'exp.data' might be included in the mixture.
#' But the quantity matrix only reports quantity information for the cell types in 'include.in.X'.
#' @param n.per.mixture integer, 1 <= 'n.per.mixture', how many samples per type should be used for each mixture (Default: 1)
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
#' training.data <- mix_samples_with_jitter(
#'     include.in.X = include.in.X
#'     , prob.each = c(2,1,1,1)
#'     , n.samples = 1e3
#'     , exp.data = random.data
#'     , pheno = indicator.list
#'     , add.jitter = TRUE
#'     , chosen.mean = 1
#'     , chosen.sd = 0.05
#' )
#'
mix_samples_with_jitter <- function(
  include.in.X, # tested
  prob.each = NA, # tested
  n.samples, # tested
  exp.data, # tested
  pheno, # tested
  verbose = FALSE, # tested
  add.jitter = FALSE, # tested
  chosen.mean = 1, # tested
  chosen.sd = 0.05, # tested
  n.per.mixture = 1 # tested
){
  # Safety check: pheno, exp.data
  if(!is.vector(include.in.X)){
    stop("in mix_samples_with_jitter: 'include.in.X' is not provided as vector")
  }
  if(!is.vector(pheno)){
    stop("in mix_samples_with_jitter: 'pheno' is not provided as vector")
  }
  if(!any(include.in.X %in% pheno)){
    stop("in mix_samples_with_jitter: no cell type in 'include.in.X' fits 'pheno'")
  }

  if(!(is.matrix(exp.data) && is.numeric(exp.data))){
    stop("in mix_samples_with_jitter: 'exp.data' is not a numeric matrix")
  }

  if(!(all(names(pheno) %in% colnames(exp.data)) && length(pheno) == ncol(exp.data))){
    stop("in mix_samples_with_jitter: 'names(pheno)' do not fit 'colnames(exp.data)'.
         For every entry of 'colnames(exp.data)' there has to be a entry in 'pheno'")
  }
  # end -> pheno, exp.data

  if(any(is.na(prob.each))){
    prob.each <- rep(1, length(include.in.X))
  }else{
    if(any(!is.numeric(prob.each))){
      stop("in mix_samples_with_jitter: there are non numeric entries in 'prob.each'.")
    }
  }
  prob.each <- prob.each/sum(prob.each)

  len <- length(include.in.X)

  # Safety check: verbose
  test <- test_logical(test.value = verbose,
                       output.info = c("mix_samples_with_jitter", "verbose"))
  # end -> verbose
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
  rownames(quant.matrix) <- include.in.X

  for(lrun in 1:n.samples){
    if(verbose){
      cat("doing ", lrun, " or ", 100*lrun/n.samples, "%\n" )
    }

    quantities <- sapply(prob.each, function(x){
      return(
        stats::runif(
          n = 1
          , min = 0
          , max = x
        )
      )
    })

    # add quantities to quant.matrix:
    quant.matrix <- cbind(quant.matrix, quantities)
    colnames(quant.matrix)[lrun] <- paste0("mix", lrun)

    # Now, using the sampled quantities, calculated the expression profiles of the mixtures.
    # Initialise a numerc vector with zeros:
    mixture <- rep(0, nrow(exp.data))
    # Go through all samples:
    for(lsample in include.in.X){
      # Find all cells in the pheno which match to "lsample"
      potentialCells <- names(pheno[pheno == lsample])
      # sample a subset of the potentiallCells:
      chosen.sample <- sample(
        x = potentialCells
        , size = n.per.mixture
        , replace = TRUE
      )
      # and multiple these with their quantity:
      next.mix <- quant.matrix[lsample, lrun] * rowSums(exp.data[, chosen.sample, drop = FALSE])
      # jitter means to multply every entry of every "next.mix" with a random number close to 1:
      if(add.jitter){
        # get a vector of random numbers:
        factors <- stats::rnorm(
          n=length(next.mix),
          mean=chosen.mean,
          sd=chosen.sd
        )
        # elementwise multiply next.mix with it
        next.mix <- next.mix * factors
      }
      # add the new "next.mix" vector to the mixture:
      mixture <- mixture + next.mix
    }
    # add the mixture to the mix.matrix, and add the colname:
    mix.matrix <- cbind(mix.matrix, mixture)
    colnames(mix.matrix)[lrun] <- paste0("mix", lrun)
  }

  # only keep the quantity information of cells that are included in X
  if(!any(is.na(include.in.X))){
    quant.matrix <- quant.matrix[include.in.X, ]
  }
  # build a return list, and add both matrices:
  rets <- vector(mode="list")
  rets[["mixtures"]] <- mix.matrix
  rets[["quantities"]] <- quant.matrix
  return(rets)
  }
