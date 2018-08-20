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
#'  library(DTD)
#' random.data <- generate.random.data(nTypes = 5,
#'                                     nSamples.perType = 10,
#'                                     nFeatures = 100,
#'                                     sample.type = "Cell",
#'                                     feature.type = "gene")
#'
#' # undo log transformation to have counts:
#'
#' random.data <- (2^random.data) - 1
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
    stop("Set input parameters (nTypes, nSamples.perType, nFeatures) to valid numeric values")
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

  n.totalSamples <- nTypes * nSamples.perType
  expression.matrix <- matrix(NA, nrow = nFeatures, ncol=n.totalSamples)

  # get column names ...
  type.names <- c()
  for(l.name in 1:nTypes){
    type.names <- c(type.names, rep(paste0("Type", l.name), nSamples.perType))
  }
  sample.names <- paste0(sample.type, 1:n.totalSamples, ".", type.names)

  # ... and set them:
  colnames(expression.matrix) <- sample.names
  rownames(expression.matrix) <- paste0(feature.type, 1:nFeatures)


  # Expression will be generated randomly using a poisson distribution.
  # Each gene in every type is poisson distributed with an unique lambda.
  # Here, I generate the lambda for each gene as a normal distributed random variable:
  lambda.mat <- matrix(data = abs(stats::rnorm(nTypes*nFeatures, mean = 5, sd = 2)),
                       nrow = nFeatures,
                       ncol = nTypes)
  rownames(lambda.mat) <- rownames(expression.matrix)
  colnames(lambda.mat) <- paste0("Type", 1:nTypes)


  # Now, I loop over every type ...
  for(l.type in colnames(lambda.mat)){
    # ... and every gene ...
    for(l.gene in rownames(lambda.mat)){
      # (There are many samples per type: )
      pos <- which(grepl(x = colnames(expression.matrix),
                         pattern = l.type))

      # ... and sample "nSamples.perType" expression values, with the previously generated lambda:
      expression.matrix[l.gene, pos] <- stats::rpois(n = nSamples.perType,
                                                     lambda = lambda.mat[l.gene, l.type])
    }
  }

  return(expression.matrix)
}
