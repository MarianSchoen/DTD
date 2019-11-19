#' Generating random data
#'
#' 'generate_random_data' simulates data which can be used exemplary for
#' digital tissue deconvolution. It will generate a numeric matrix with '
#' n.features' rows, and ('n.types' * 'n.samples.per.type') columns.
#' Each column represents a sample of special type. The function will generate
#' 'n.types', and for each type 'n.samples.per.type'.\cr
#' Mathematically, each feature is drawn from a poisson distribution.
#' For each feature in every cell type, a lambda is drawn randomly.
#' Then it generates multiple samples per type. This ensures that samples
#' from the same cell type have similar counts for the same feature.
#'
#' @param n.types integer, 2 <= 'n.types', how many different types should be
#' included in the data set
#' @param n.samples.per.type integer 1 <= 'n.samples.per.type',
#' how many samples should be generated per type
#' @param n.features integer, 1 <= 'n.features', how many features should be
#' included
#' @param sample.type string, name of samples
#' @param feature.type string, name of features
#' @param seed integer, will be passed to "set_seed"
#'
#' @return matrix with ('n.types' * 'n.samples.per.type') columns,
#' and 'n.features' rows
#' @export
#'
#' @examples
#' library(DTD)
#' random.data <- generate_random_data(
#'   n.types = 5,
#'   n.samples.per.type = 10,
#'   n.features = 100,
#'   sample.type = "Cell",
#'   feature.type = "gene"
#' )
#'
#' # normalize all samples to the same amount of counts:
#' random.data <- normalize_to_count(random.data)
generate_random_data <- function(n.types = 5,
                                 n.samples.per.type = 10,
                                 n.features = 1000,
                                 sample.type = "Cell",
                                 feature.type = "gene",
                                 seed = 1310) {

  # safety check: n.types
  test <- test_integer(test.value = n.types,
                       output.info = c("generate_random_data", "n.types"),
                       min = 2,
                       max = Inf)
  # end -> n.types

  # safety check: n.samples.per.type
  test <- test_integer(test.value = n.samples.per.type,
                       output.info = c("generate_random_data", "n.samples.per.type"),
                       min = 1,
                       max = Inf)
  # end -> n.samples.per.type

  # safety check: n.features
  test <- test_integer(test.value = n.features,
                       output.info = c("generate_random_data", "n.features"),
                       min = 1,
                       max = Inf)
  # end -> n.features

  if (!is.character(sample.type) || length(sample.type) != 1) {
    message("sample.type is not a single 'character', therefore set to 'Cell'\n")
    sample.type <- "Cell"
  }
  if (!is.character(feature.type) || length(feature.type) != 1) {
    message("feature.type  is not a single 'character', therefore set to 'gene'\n")
    feature.type <- "gene"
  }

  # set seed, notice that if the provided seed is not numeric, seed will be set to 1310
  if(is.numeric(seed)){
    if (round(seed) == seed || length(seed) == 1) {
      set.seed(seed)
    } else {
      message("The provided seed could not be used! Therefore, set it to default \n")
      seed <- 1310
    }
  }else{
    message("The provided seed could not be used! Therefore, set it to default \n")
    seed <- 1310
  }
  set.seed(seed)

  n.totalSamples <- n.types * n.samples.per.type
  expression.matrix <- matrix(NA, nrow = n.features, ncol = n.totalSamples)

  # get column names ...
  type.names <- c()
  for (l.name in 1:n.types) {
    type.names <- c(type.names, rep(paste0("Type", l.name), n.samples.per.type))
  }
  sample.names <- paste0(sample.type, 1:n.totalSamples, ".", type.names)

  # ... and set them:
  colnames(expression.matrix) <- sample.names
  rownames(expression.matrix) <- paste0(feature.type, 1:n.features)


  # Expression will be generated randomly using a poisson distribution.
  # Each gene in every type is poisson distributed with an unique lambda.
  # Here, I generate the lambda for each gene as a normal distributed random variable:
  lambda.mat <- matrix(
    data = abs(stats::rnorm(n.types * n.features, mean = 10, sd = 10)),
    nrow = n.features,
    ncol = n.types
  )
  rownames(lambda.mat) <- rownames(expression.matrix)
  colnames(lambda.mat) <- paste0("Type", 1:n.types)


  # Now, I loop over every type ...
  for (l.type in colnames(lambda.mat)) {
    # ... and every gene ...
    for (l.gene in rownames(lambda.mat)) {
      # (There are many samples per type: )
      pos <- which(grepl(
        x = colnames(expression.matrix),
        pattern = l.type
      ))

      # ... and sample "n.samples.per.type" expression values, with the previously generated lambda:
      expression.matrix[l.gene, pos] <- stats::rpois(
        n = n.samples.per.type,
        lambda = lambda.mat[l.gene, l.type]
      )
    }
  }

  return(expression.matrix)
}
