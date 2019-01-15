#' Generating random data
#'
#' generate_random_data simulates data which can be used exemplary for digital tissue deconvolution.
#' It will generate a numeric matrix with "n.features" rows, and "n.types" * "n.samples.per.type"
#' columns. Each column represents a sample of special type. The function will generate "n.types",
#' and for each type "n.samples.per.type".
#' Mathematically, the function randomly generates a lambda for each feature in each type.
#' Then it generates multiple samples per type. Each feature in every sample will be drawn from a poisson
#' distribution with the previously sampled lambda.
#'
#' @param n.types integer, how many different types should be included in the data set (default 5)
#' @param n.samples.per.type integer, how many samples should be generated per type (default 10)
#' Notice, for each type, the number of samples will be randomized a bit
#' @param n.features integer, how many features should be included (default 1000)
#' @param sample.type string, name of samples, defaults to "Cell"
#' @param feature.type string, name of features, defaults to "gene"
#' @param seed integer, to which the seed will be set, defaults to 1310
#'
#' @return matrix with (n.types * n.samples.per.type) columns, and n.features rows
#' @export
#'
#' @examples
#' library(DTD)
#' random.data <- generate_random_data(n.types = 5,
#'                                     n.samples.per.type = 10,
#'                                     n.features = 100,
#'                                     sample.type = "Cell",
#'                                     feature.type = "gene")
#'
#' # normalize all samples to the same amount of counts:
#' random.data <- normalize_to_count(random.data)
generate_random_data <- function(n.types = 5,
                                 n.samples.per.type = 10,
                                 n.features = 1000,
                                 sample.type = "Cell",
                                 feature.type = "gene",
                                 seed = 1310){

  # Check if input is valid:
  if(!is.numeric(n.types) || !is.numeric(n.samples.per.type) || !is.numeric(n.features)){
    stop("Set input parameters (n.types, n.samples.per.type, n.features) to valid numeric values")
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

  n.totalSamples <- n.types * n.samples.per.type
  expression.matrix <- matrix(NA, nrow = n.features, ncol=n.totalSamples)

  # get column names ...
  type.names <- c()
  for(l.name in 1:n.types){
    type.names <- c(type.names, rep(paste0("Type", l.name), n.samples.per.type))
  }
  sample.names <- paste0(sample.type, 1:n.totalSamples, ".", type.names)

  # ... and set them:
  colnames(expression.matrix) <- sample.names
  rownames(expression.matrix) <- paste0(feature.type, 1:n.features)


  # Expression will be generated randomly using a poisson distribution.
  # Each gene in every type is poisson distributed with an unique lambda.
  # Here, I generate the lambda for each gene as a normal distributed random variable:
  lambda.mat <- matrix(data = abs(stats::rnorm(n.types*n.features, mean = 10, sd = 10)),
                       nrow = n.features,
                       ncol = n.types)
  rownames(lambda.mat) <- rownames(expression.matrix)
  colnames(lambda.mat) <- paste0("Type", 1:n.types)


  # Now, I loop over every type ...
  for(l.type in colnames(lambda.mat)){
    # ... and every gene ...
    for(l.gene in rownames(lambda.mat)){
      # (There are many samples per type: )
      pos <- which(grepl(x = colnames(expression.matrix),
                         pattern = l.type))

      # ... and sample "n.samples.per.type" expression values, with the previously generated lambda:
      expression.matrix[l.gene, pos] <- stats::rpois(n = n.samples.per.type,
                                                     lambda = lambda.mat[l.gene, l.type])
    }
  }

  return(expression.matrix)
}
