#' Generating random data
#'
#' generate.random.data simulates data which can be used for digital tissue deconvolution.
#' It will generate a numeric matrix with "nFeatures" rows, and "nTypes" * "nSamples.perType"
#' columns. Each column represents a sample of special type. The function will generate "nTypes",
#' and for each type "nSamples.perType". Notice, that "nSamples.perType" will be randomised,
#' such that the number of sampler per type differ between the types.
#'
#' @param nTypes integer, how many different types should be included in the data set
#' @param nSamples.perType integer, how many samples should be generated per type.
#' Notice, for each type, the number of samples will be randomized a bit.
#' @param nFeatures integer, how many features should be included
#' @param sample.type string, name of samples
#' @param feature.type string, name of features
#' @param seed integer, to which the seed will be set
#'
#' @return matrix with ~ (nTypes * nSamples.perType) columns, and nFeatures rows
#' @export
#'
#' @examples
#' #' library(DTD)
#' random.data <- generate.random.data(nTypes = 5,
#'                                     nSamples.perType = 10,
#'                                     nFeatures = 100,
#'                                     sample.type = "Cell",
#'                                     feature.type = "gene")
#'
#' # normalize all samples to the same amount of counts:
#' random.data <- normalizeToCount(random.data)
generate.random.data <- function(nTypes = 5,
                                 nSamples.perType = 10,
                                 nFeatures = 1000,
                                 sample.type = "Cell",
                                 feature.type = "gene",
                                 seed = 1310){

  # Check if input is valid:
  if(!is.numeric(nTypes) || !is.numeric(nSamples.perType) || !is.numeric(nFeatures)){
    stop("Set input parameters (nTypes, nSamples.perType, nFeatures) to valid values")
  }

  if(!is.character(sample.type)){
    sample.type <- "Cell"
  }
  if(!is.character(feature.type)){
    feature.type <- "gene"
  }

  # set seed, notice that if the provided seed is not numeric, seed will be set to 1310
  if(is.numeric(seed)){
    set.seed(seed)
  }else{
    cat("The provided seed could not be used! Therefore, set it to default \n")
    set.seed(1310)
  }



  # fluctuate number of samples per Type
  cells.perType <- abs(round(rnorm(n=nTypes, mean = nSamples.perType, sd = 1)))

  # calculate total number of samples
  n.totalSamples <- sum(cells.perType)

  # initalising the expression matrix with names:
  expression.matrix <- matrix(NA, nrow = nFeatures, ncol=n.totalSamples)

  # calculate column names ...
  type.names <- c()
  for(l.name in 1:nTypes){
    type.names <- c(type.names, rep(paste0("Type", l.name), cells.perType[l.name]))
  }
  sample.names <- paste0(sample.type, 1:n.totalSamples, ".", type.names)

  # ... and set them:
  colnames(expression.matrix) <- sample.names
  rownames(expression.matrix) <- paste0(feature.type, 1:nFeatures)


  # Data is generated via rnorm.
  # For each cell type the mean and standard deviation changes.
  # Therefore, for each cell type a mean and a sd is randomly generated
  means <- abs(rnorm(n = nTypes, mean = 5, sd=25))
  sds <- abs(rnorm(n = nTypes, mean = 5, sd=50))

  # means and sds are a list of length nTypes.
  # To have the same mean/sds for each sample of the same type each entry
  # of means will be replicated as often as there are cells.perType
  sample.means <- c()
  sample.sds <- c()
  for(l.sample in 1:nTypes){
    sample.means <- c(sample.means, rep(means[l.sample], cells.perType[l.sample]))
    sample.sds <- c(sample.sds, rep(sds[l.sample], cells.perType[l.sample]))
  }

  # Here, the actual data generation starts.
  # Each sample will be simulated iteratively:
  for(l.sample in 1:n.totalSamples){
    # expression will be a positive value --> abs
    # all values are drawn from a normal distribution with already generated means and sd:
    simulated.expression <- abs(
                            rnorm(n = nFeatures,
                                  mean=sample.means[l.sample],
                                  sd = sample.sds[l.sample]
                                  )
                                )
    # replace column of expression matrix with the simulated.expression
    expression.matrix[, l.sample] <- simulated.expression
  }
  # return expression matrix
  return(expression.matrix)
}
