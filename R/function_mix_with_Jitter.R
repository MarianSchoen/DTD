#' Mix samples with jitter
#'
#' 'mix_samples_with_jitter' takes a expression matrix ('expr.data')
#' and 'pheno' information.
#' It then mixes the samples with known quantities, such that it can be used for
#' loss-function learning digital tissue deconvolution.
#' For a mixture it randomly selects a quantity for each cell type. Then, it
#' randomly selects profiles from 'expr.data' for each cell type, multiplies it
#' with the respective quantity, and averages it to a simulated bulk profile.
#' Notice, in the mixtures, the frequency of a cell type is reflected by the
#' 'prob.each' vector, and not by their occurrence in 'pheno'. Alternatively,
#' there is a 'mix_samples' function, that reflects the underlying composition
#' of the data set: \code{\link{mix_samples}}
#'
#' @param n.samples integer above 0, numbers of samples to be drawn
#' @param prob.each numeric vector with same length as 'included.in.X.'
#' For each cell type in 'included.in.X', prob.each' holds the expected
#' average quantity in the mixtures.
#' @param expr.data numeric matrix, with features as rows and samples as columns
#' @param add.jitter logical, should each mixture be multiplied with a vector
#' of normally distributed numbers? (JITTER)
#' @param verbose logical, should information be printed to console?
#' @param pheno named vector of strings, with pheno information ('pheno')
#' for each sample in 'expr.data'. names(pheno)' must all be in
#' 'colnames(expr.data)'
#' @param chosen.mean float, mean of jitter
#' @param chosen.sd float, standard deviation of jitter
#' @param included.in.X vector of strings, indicating types that are in the
#' reference matrix. Only those types, and sorted in that order, will be
#' included in the quantity matrix.
#' @param n.per.mixture integer above 0, below ncol(expr.data),
#'  how many samples should be included per mixutre
#' @param normalize.to.count logical, normalize each mixture?
#'
#' @return list with two entries. "quantities": matrix (nrow = ncol(expr.data), ncol = n.samples)
#' and "mixtures": matrix (nrow = nrow(expr.data), ncol = n.samples)
#'
#' @export
#'
#' @examples
#' library(DTD)
#' random.data <- generate_random_data(
#'       n.types = 10,
#'       n.samples.per.type = 10,
#'       n.features = 500,
#'       sample.type = "Cell",
#'       feature.type = "gene"
#'       )
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
#' included.in.X <- c("Type2", "Type3", "Type4", "Type5")
#'
#' training.data <- mix_samples_with_jitter(
#'     included.in.X = included.in.X
#'     , prob.each = c(2,1,1,1)
#'     , n.samples = 1e3
#'     , expr.data = random.data
#'     , pheno = indicator.list
#'     , add.jitter = TRUE
#'     , chosen.mean = 1
#'     , chosen.sd = 0.05
#' )
#'
#' # see the effect of "Type2" having higher 'prob.each' entry:
#' apply(training.data$quantities, 1, mean)
#'
mix_samples_with_jitter <- function(
  included.in.X, # tested
  prob.each = NA, # tested
  n.samples, # tested
  expr.data, # tested
  pheno, # tested
  verbose = FALSE, # tested
  add.jitter = FALSE, # tested
  chosen.mean = 1, # tested
  chosen.sd = 0.05, # tested
  n.per.mixture = 1, # tested
  normalize.to.count = TRUE
){
  # Safety check: pheno, expr.data
  if(!is.vector(included.in.X)){
    stop("in mix_samples_with_jitter: 'included.in.X' is not provided as vector")
  }
  if(!is.vector(pheno)){
    stop("in mix_samples_with_jitter: 'pheno' is not provided as vector")
  }
  if(!any(included.in.X %in% pheno)){
    stop("in mix_samples_with_jitter: no cell type in 'included.in.X' fits 'pheno'")
  }

  if(!(is.matrix(expr.data) && is.numeric(expr.data))){
    stop("in mix_samples_with_jitter: 'expr.data' is not a numeric matrix")
  }

  if(!(all(names(pheno) %in% colnames(expr.data)) && length(pheno) == ncol(expr.data))){
    stop("in mix_samples_with_jitter: 'names(pheno)' do not fit 'colnames(expr.data)'.
         For every entry of 'colnames(expr.data)' there has to be a entry in 'pheno'")
  }
  test <- test_logical(
    test.value = normalize.to.count,
    output.info = c("mix_samples_with_jitter", "normalize.to.count")
  )
  # end -> pheno, expr.data

  if(any(is.na(prob.each))){
    warning("in mix_samples_with_jitter: There are 'NA' in 'prob.each', therefore set it to a all 1 vector")
    prob.each <- rep(1, length(included.in.X))
  }else{
    if(any(!is.numeric(prob.each))){
      stop("in mix_samples_with_jitter: there are non numeric entries in 'prob.each'.")
    }
    if(length(prob.each) != length(included.in.X)){
      stop("in mix_samples_with_jitter: 'length(prob.each)' unequal to 'length(included.in.X)'.")
    }

  }
  prob.each <- prob.each/sum(prob.each)

  len <- length(included.in.X)

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
        "In 'mix_samples_with_jitter': for some types in 'included.in.X' there are no samples in your data. Provide e.g. the following vector as 'included.in.X': \n",
        paste0(
          useable.types
          , collapse = ", "
        )
      )
    )
  }

  # initialise return variables:
  mix.matrix <- matrix(nrow=nrow(expr.data), ncol=0)
  rownames(mix.matrix) <- rownames(expr.data)

  quant.matrix <- matrix(nrow=len, ncol=0)
  rownames(quant.matrix) <- included.in.X

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
    mixture <- rep(0, nrow(expr.data))
    # Go through all samples:
    for(lsample in included.in.X){
      # Find all cells in the pheno which match to "lsample"
      potentialCells <- names(pheno[which(pheno == lsample)])
      # sample a subset of the potentiallCells:
      chosen.sample <- sample(
        x = potentialCells
        , size = n.per.mixture
        , replace = TRUE
      )
      # and multiple these with their quantity:
      next.mix <- quant.matrix[lsample, lrun] * rowSums(expr.data[, chosen.sample, drop = FALSE])
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
  if(!any(is.na(included.in.X))){
    quant.matrix <- quant.matrix[included.in.X, ]
  }

  if(normalize.to.count){
    mix.matrix <- normalize_to_count(mix.matrix)
  }
  # build a return list, and add both matrices:
  rets <- vector(mode="list")
  rets[["mixtures"]] <- mix.matrix
  rets[["quantities"]] <- quant.matrix
  return(rets)
  }
