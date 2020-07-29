# written by Marian Schön
rm(list = ls())
gc()
if ("DTD" %in% (.packages())) {
  detach("package:DTD", unload = TRUE)
}
# devtools::load_all()
library(DTD)
library(testthat)

context("Test generate_random_data ")

test_that("generate_random_data ", {
  expect_error(
    generate_random_data(n.types = "25"),
    "In generate_random_data: 'n.types' is not a single integer"
  )
  expect_error(
    generate_random_data(n.types = 2.2),
    "In generate_random_data: 'n.types' is not an integer"
  )
  expect_error(
    generate_random_data(n.types = c(1, 2)),
    "In generate_random_data: 'n.types' is not a single integer"
  )
  expect_error(
    generate_random_data(n.types = 0),
    "In generate_random_data: 'n.types' is below minimal value"
  )
  expect_message(
    generate_random_data(sample.type = c("Cell", "gene")),
    "sample.type is not a single 'character', therefore set to 'Cell'"
  )
  expect_message(
    generate_random_data(feature.type = c("Cell", "gene")),
    "feature.type  is not a single 'character', therefore set to 'gene'"
  )
  expect_message(
    generate_random_data(seed = c("iwas2", 2)),
    "The provided seed could not be used! Therefore, set it to default"
  )
})
{
  random.data <- generate_random_data(
    n.types = 25,
    n.samples.per.type = 100,
    n.features = 100,
    seed = 102
  )
}
context("Test normalize_to_count ")
test_that("normalize_to_count ", {
  expect_error(
    normalize_to_count(expr.data = random.data, count = "iwas"),
    "In normalize_to_count: 'count' is not a single integer"
  )
  expect_error(
    normalize_to_count(expr.data = random.data, count = c(1, 2, 3)),
    "In normalize_to_count: 'count' is not a single integer"
  )
  expect_error(
    normalize_to_count(expr.data = 1, count = 1),
    "In normalize_to_count: expr.data is not of matrix format"
  )
  expect_error(
    normalize_to_count(expr.data = matrix(-9:-1, nrow = 3)),
    "In normalize_to_count: expr.data includes negative values"
  )
  expect_error(
    normalize_to_count(expr.data = matrix(c(1:8, NA), nrow = 3)),
    "In normalize_to_count: expr.data includes NA"
  )
})
{
  normalized.data <- normalize_to_count(random.data)

  indicator.list <- gsub("^Cell[0-9]*\\.", "", colnames(normalized.data))
  names(indicator.list) <- colnames(normalized.data)
  include.in.X <- paste0("Type", 1:5)
  perc.of.all.cells <- 0.1
}
context("Test sample_random_X ")
test_that("sample_random_X", {
  expect_error(
    sample_random_X(
      included.in.X = 1,
      pheno = indicator.list,
      expr.data = normalized.data,
      percentage.of.all.cells = perc.of.all.cells
    ),
    "in sample_random_X: no cell type in 'included.in.X' fits 'pheno'"
  )
  expect_error(
    sample_random_X(
      included.in.X = c("1", "2", "3"),
      pheno = indicator.list,
      expr.data = normalized.data,
      percentage.of.all.cells = perc.of.all.cells
    ),
    "in sample_random_X: no cell type in 'included.in.X' fits 'pheno'"
  )
  expect_error(
    sample_random_X(
      included.in.X = include.in.X,
      pheno = 1,
      expr.data = normalized.data,
      percentage.of.all.cells = perc.of.all.cells
    ),
    "in sample_random_X: no cell type in 'included.in.X' fits 'pheno'"
  )
  expect_error(
    sample_random_X(
      included.in.X = include.in.X,
      pheno = indicator.list[1:10],
      expr.data = normalized.data,
      percentage.of.all.cells = perc.of.all.cells
    ),
    "in sample_random_X: 'names(pheno)' do not fit 'colnames(expr.data)'. For every entry of 'colnames(expr.data)' there has to be an entry in 'pheno'",
    fixed = TRUE
  )
  expect_error(
    sample_random_X(
      included.in.X = include.in.X,
      pheno = indicator.list,
      expr.data = 1,
      percentage.of.all.cells = perc.of.all.cells
    ),
    "in sample_random_X: 'names(pheno)' do not fit 'colnames(expr.data)'. For every entry of 'colnames(expr.data)' there has to be an entry in 'pheno'",
    fixed = TRUE
  )
  expect_error(
    sample_random_X(
      included.in.X = include.in.X,
      pheno = indicator.list,
      expr.data = normalized.data,
      percentage.of.all.cells = "1"
    ),
    "in sample_random_X: 'percentage.of.all.cells' is not a single numeric"
  )
  expect_error(
    bol.sam.ram.X <- sample_random_X(
      included.in.X = include.in.X,
      pheno = indicator.list,
      expr.data = normalized.data,
      percentage.of.all.cells = c(1, 2, 3)
    ),
    "in sample_random_X: 'percentage.of.all.cells' is not a single numeric"
  )
  expect_error(
    sample_random_X(
      included.in.X = include.in.X,
      pheno = indicator.list,
      expr.data = normalized.data,
      percentage.of.all.cells = matrix(1:4, nrow = 2)
    ),
    "in sample_random_X: 'percentage.of.all.cells' is not a single numeric"
  )
  expect_error(
    sample_random_X(
      included.in.X = include.in.X,
      pheno = indicator.list,
      expr.data = normalized.data,
      percentage.of.all.cells = -1
    ),
    "in sample_random_X: 'percentage.of.all.cells' must be above 0 and below 1"
  )
  expect_error(
    sample_random_X(
      included.in.X = include.in.X,
      pheno = indicator.list,
      expr.data = normalized.data,
      percentage.of.all.cells = 2
    ),
    "in sample_random_X: 'percentage.of.all.cells' must be above 0 and below 1"
  )
})
{
  sample.X <- sample_random_X(
    included.in.X = include.in.X,
    pheno = indicator.list,
    expr.data = normalized.data,
    percentage.of.all.cells = perc.of.all.cells
  )

  # the following are calls independent of DTD, therefore no tests
  X.matrix <- sample.X$X.matrix
  genes <- rownames(X.matrix)


  # here, I sample 75% of the features, making them useless for deconvolution:
  X.matrix[1:(0.75 * nrow(X.matrix)), ] <- sample(X.matrix[1:(0.75 * nrow(X.matrix)), ])

  samples.to.remove <- sample.X$samples.to.remove

  remaining.mat <- normalized.data[, -which(colnames(normalized.data) %in% samples.to.remove)]

  # sampling training samples: (notice, that train test seperation is 50:50)
  train.samples <- sample(
    x = colnames(remaining.mat),
    size = ceiling(ncol(remaining.mat) / 2),
    replace = FALSE
  )
  # selecting test samples:
  test.samples <- colnames(remaining.mat)[which(!colnames(remaining.mat) %in% train.samples)]

  # extract data matrices for training and testing:
  train.mat <- remaining.mat[, train.samples]
  test.mat <- remaining.mat[, test.samples]

  indicator.train <- indicator.list[names(indicator.list) %in% colnames(train.mat)]
}
context("Test mix_with_jitter")
{
  major.fraction.type <- c("Type1")
  all.samples <- unique(indicator.list)
  minor.fraction.type <- all.samples[-which(all.samples %in% major.fraction.type)]
}
test_that("mix_with_jitter", {
  expect_error(
    mix_samples_with_jitter(
      included.in.X = "include.in.X",
      prob.each = rep(1, length(include.in.X)),
      n.per.mixture = 200,
      normalize.to.count = FALSE,
      n.samples = 1e3,
      expr.data = random.data,
      pheno = indicator.list,
      add.jitter = TRUE,
      chosen.mean = 1,
      chosen.sd = 0.05,
    ),
    "in mix_samples_with_jitter: no cell type in 'included.in.X' fits 'pheno'"
  )
  expect_error(
    mix_samples_with_jitter(
      included.in.X = include.in.X,
      prob.each = rep(1, length(include.in.X)),
      n.per.mixture = 200,
      normalize.to.count = FALSE,
      n.samples = -1e3,
      expr.data = random.data,
      pheno = indicator.list,
      add.jitter = TRUE,
      chosen.mean = 1,
      chosen.sd = 0.05,
    ),
    "In mix_samples_with_jitter: 'n.samples' is below minimal value"
  )
  expect_error(
    mix_samples_with_jitter(
      included.in.X = include.in.X,
      prob.each = rep(1, length(include.in.X)),
      n.per.mixture = 200,
      normalize.to.count = FALSE,
      n.samples = 1e3,
      expr.data = "random.data",
      pheno = indicator.list,
      add.jitter = TRUE,
      chosen.mean = 1,
      chosen.sd = 0.05,
    ),
    "in mix_samples_with_jitter: 'expr.data' is not a numeric matrix"
  )
  expect_error(
    mix_samples_with_jitter(
      included.in.X = include.in.X,
      prob.each = rep(1, length(include.in.X)),
      n.per.mixture = 200,
      normalize.to.count = FALSE,
      n.samples = 1e3,
      pheno = indicator.list,
      add.jitter = TRUE,
      chosen.mean = 1,
      chosen.sd = 0.05,
      expr.data = cbind(random.data, as.character(random.data[, 1]))
    ),
    "in mix_samples_with_jitter: 'expr.data' is not a numeric matrix"
  )
  expect_error(
    mix_samples_with_jitter(
      included.in.X = include.in.X,
      prob.each = rep(1, length(include.in.X)),
      n.per.mixture = 200,
      normalize.to.count = FALSE,
      n.samples = 1e3,
      expr.data = random.data,
      pheno = "indicator.list",
      add.jitter = TRUE,
      chosen.mean = 1,
      chosen.sd = 0.05,
    ),
    "in mix_samples_with_jitter: no cell type in 'included.in.X' fits 'pheno'"
  )
})
context("Test mix_samples ")
test_that("mix_samples ", {
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = "0",
      n.per.mixture = 100,
      verbose = FALSE
    ),
    "In mix_samples: 'n.samples' is not a single integer"
  )
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = 0,
      n.per.mixture = 100,
      verbose = FALSE
    ),
    "In mix_samples: 'n.samples' is below minimal value"
  )
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = 3.1415,
      n.per.mixture = 100,
      verbose = FALSE
    ),
    "In mix_samples: 'n.samples' is not an integer"
  )
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = c(1, 2, 3, -3.1415),
      n.per.mixture = 100,
      verbose = FALSE
    ),
    "In mix_samples: 'n.samples' is not a single integer"
  )
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = length(genes),
      n.per.mixture = "10",
      verbose = FALSE
    ),
    "In mix_samples: 'n.per.mixutre' is not a single integer"
  )
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = length(genes),
      n.per.mixture = -10,
      verbose = FALSE
    ),
    "In mix_samples: 'n.per.mixutre' is below minimal value"
  )
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = length(genes),
      n.per.mixture = c(1, 2, 3),
      verbose = FALSE
    ),
    "In mix_samples: 'n.per.mixutre' is not a single integer"
  )
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = length(genes),
      n.per.mixture = ncol(train.mat) + 1,
      verbose = FALSE
    ),
    "In mix_samples: 'n.per.mixutre' is above maximal value"
  )
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = length(genes),
      n.per.mixture = 100,
      verbose = 1
    ),
    "In mix_samples: 'verbose' must be a single value, either 'TRUE' or 'FALSE' (not 1 or 0)",
    fixed = TRUE
  )
  expect_error(
    mix_samples(
      expr.data = train.mat,
      pheno = indicator.train,
      included.in.X = include.in.X,
      n.samples = length(genes),
      n.per.mixture = 100,
      verbose = c(TRUE, FALSE, FALSE)
    ),
    "In mix_samples: 'verbose' must be a single value, either 'TRUE' or 'FALSE' (not 1 or 0)",
    fixed = TRUE
  )
})
{
  training.data <- mix_samples(
    expr.data = train.mat,
    pheno = indicator.train,
    included.in.X = include.in.X,
    n.samples = length(genes),
    n.per.mixture = 100,
    verbose = FALSE
  )

  indicator.test <- indicator.list[names(indicator.list) %in% colnames(test.mat)]
  test.data <- mix_samples(
    expr.data = test.mat,
    pheno = indicator.test,
    included.in.X = include.in.X,
    n.samples = length(genes),
    n.per.mixture = 100,
    verbose = FALSE
  )

  start.tweak <- rep(1, nrow(X.matrix))
  names(start.tweak) <- rownames(X.matrix)
}
context("Test evaluate_cor ")
test_that("evaluate_cor ", {
  expect_error(
    evaluate_cor(
      X.matrix = NA,
      new.data = training.data$mixtures,
      true.compositions = training.data$quantities,
      DTD.model = start.tweak,
      estimate.c.type = "direct"
    ),
    "In evaluate_cor: 'X.matrix' must be provided either as the 'X.matrix' argument, or within the DTD.model"
  )
  expect_error(
    evaluate_cor(
      X.matrix = X.matrix,
      new.data = training.data$mixtures,
      true.compositions = training.data$quantities,
      DTD.model = c(start.tweak, "ichbinkeinenumeric"),
      estimate.c.type = "direct"
    ),
    "In evaluate_cor: 'DTD.model' does not fit 'X.matrix'. Check if the provided model includes as many 'Tweak' entries as there are features (=> rows) in 'X.matrix'",
    fixed = TRUE
  )
  expect_error(
    evaluate_cor(
      X.matrix = X.matrix,
      new.data = training.data$mixtures,
      true.compositions = training.data$quantities,
      DTD.model = as.character(start.tweak),
      estimate.c.type = "direct"
    ),
    "In evaluate_cor: Used 'Tweak' (from 'DTD.model') is not numeric",
    fixed = TRUE
  )
  expect_error(
    evaluate_cor(
      X.matrix = X.matrix,
      new.data = training.data$mixtures,
      true.compositions = training.data$quantities,
      DTD.model = start.tweak,
      estimate.c.type = "iwasanders"
    )
  )

  expect_equal(
    evaluate_cor(
      X.matrix = X.matrix[sample(1:nrow(X.matrix)), sample(1:ncol(X.matrix))],
      new.data = training.data$mixtures[
        sample(1:nrow(training.data$mixtures)),
        sample(1:ncol(training.data$mixtures))
      ],
      true.compositions = training.data$quantities[
        sample(1:nrow(training.data$quantities)),
        sample(1:ncol(training.data$quantities))
      ],
      DTD.model = sample(start.tweak),
      estimate.c.type = "direct"
    ) / ncol(X.matrix),
    evaluate_cor(
      X.matrix = X.matrix[sample(1:nrow(X.matrix)), sample(1:ncol(X.matrix))],
      new.data = training.data$mixtures[
        sample(1:nrow(training.data$mixtures)),
        sample(1:ncol(training.data$mixtures))
      ],
      true.compositions = training.data$quantities[
        sample(1:nrow(training.data$quantities)),
        sample(1:ncol(training.data$quantities))
      ],
      DTD.model = sample(start.tweak),
      estimate.c.type = "direct"
    ) / ncol(X.matrix)
  )
})

context("Test train_deconvolution_model ")
test_that("train_deconvolution_model ", {
  expect_error(
    train_deconvolution_model(
      tweak = "hallo",
      X.matrix = X.matrix,
      train.data.list = training.data,
      test.data.list = test.data,
      estimate.c.type = "direct"
    ),
    "In train_deconvolution_model: 'tweak' does not fit 'X.matrix'. 'length(tweak)' must be 'nrow(X.matrix)'",
    fixed = TRUE
  )
  expect_error(
    train_deconvolution_model(
      tweak = 1:10,
      X.matrix = X.matrix,
      train.data.list = training.data,
      test.data.list = test.data,
      estimate.c.type = "direct"
      ,
    ),
    "In train_deconvolution_model: 'tweak' does not fit 'X.matrix'. 'length(tweak)' must be 'nrow(X.matrix)'",
    fixed = TRUE
  )
  tmp <- start.tweak
  names(tmp) <- paste0(1:length(tmp))
  expect_error(
    train_deconvolution_model(
      tweak = tmp,
      X.matrix = X.matrix,
      train.data.list = training.data,
      test.data.list = test.data,
      estimate.c.type = "direct"
    ),
    "In train_deconvolution_model: There are features within 'X.matrix' where no entry in 'tweak' can be found"
  )
  tmp <- start.tweak
  tmp[1] <- NA
  expect_error(
    train_deconvolution_model(
      tweak = tmp,
      X.matrix = X.matrix,
      train.data.list = training.data,
      test.data.list = test.data,
      estimate.c.type = "direct"
    ),
    "In train_deconvolution_model: 'tweak' includes NA"
  )

  expect_error(
    train_deconvolution_model(
      tweak = start.tweak,
      X.matrix = X.matrix,
      train.data.list = NA,
      test.data.list = test.data,
      estimate.c.type = "direct"
    ),
    "In train_deconvolution_model: train.data.list must be provided as a list with two entries: 'quantities' and 'mixtures'"
  )
  tmp <- training.data
  tmp$quantities <- NULL
  expect_error(
    train_deconvolution_model(
      tweak = start.tweak,
      X.matrix = X.matrix,
      train.data.list = tmp,
      test.data.list = test.data,
      estimate.c.type = "direct"
    ),
    "In train_deconvolution_model: train.data.list must be provided as a list with two entries: 'quantities' and 'mixtures'"
  )
  tmp <- training.data
  tmp$quantities <- rep(1:10)
  expect_error(
    train_deconvolution_model(
      tweak = start.tweak,
      X.matrix = X.matrix,
      train.data.list = tmp,
      test.data.list = test.data,
      estimate.c.type = "direct"
    ),
    "In train_deconvolution_model: 'train.data.list$quantities' is not a matrix",
    fixed = TRUE
  )
  tmp <- training.data
  tmp$quantities <- NULL
  tmp$wasanders <- 1
  expect_error(
    train_deconvolution_model(
      tweak = start.tweak,
      X.matrix = X.matrix,
      train.data.list = tmp,
      test.data.list = test.data,
      estimate.c.type = "direct"
    ),
    "In train_deconvolution_model: entries of train.data.list must be named 'quantities' and 'mixtures'"
  )
  expect_error(
    train_deconvolution_model(
      tweak = start.tweak,
      X.matrix = X.matrix,
      train.data.list = training.data,
      test.data.list = test.data,
      estimate.c.type = "schmuh"
    ),
    "In train_deconvolution_model: 'estimate.c.type' does not match 'non_negative' or 'direct."
  )
  expect_error(
    train_deconvolution_model(
      tweak = start.tweak,
      X.matrix = X.matrix,
      train.data.list = training.data,
      test.data.list = test.data,
      estimate.c.type = "direct",
      use.implementation = "nichtRundnichtC"
    ),
    "cannot use implementation \" nichtRundnichtC \": not implemented."
  )
})

catch <- train_deconvolution_model(
  tweak = start.tweak,
  X.matrix = X.matrix,
  train.data.list = training.data,
  test.data.list = test.data,
  estimate.c.type = "direct"#, verbose = TRUE
)


context("Test DTD_cv_lambda_cxx ")
test_that("DTD_cv_lambda_cxx ", {
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = "1",
      n.folds = 10,
      lambda.length = 10,
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = FALSE,
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: 'tweak.start' is not numeric"
    , fixed = TRUE
  )
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = c(1, 2, 3),
      n.folds = 10,
      lambda.length = 10,
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = FALSE,
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: 'nrow(train.data.list$mixtures)' does not match 'length(tweak.start)'",
    fixed = TRUE
  )
  tmp <- start.tweak
  names(tmp)[1] <- "iwas anderes"
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = tmp,
      n.folds = 10,
      lambda.length = 10,
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = FALSE,
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: there are features in 'train.data.list$mixtures' that don't match with 'names(tweak.start)'",
    fixed = TRUE
  )
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = start.tweak,
      n.folds = "10",
      lambda.length = 10,
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = FALSE,
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: 'n.folds' is not a single integer"
  )
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = start.tweak,
      n.folds = c(1, 2, 3, "4"),
      lambda.length = 10,
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = FALSE,
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: 'n.folds' is not a single integer"
  )
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = start.tweak,
      n.folds = 2.72,
      lambda.length = 10,
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = FALSE,
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: 'n.folds' is not an integer"
  )
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = start.tweak,
      n.folds = 10,
      lambda.length = "10",
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = FALSE,
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: 'lambda.length' is not a single integer"
  )
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = start.tweak,
      n.folds = 10,
      lambda.length = c(1, 2, 3, "10"),
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = FALSE,
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: 'lambda.length' is not a single integer"
  )
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = start.tweak,
      n.folds = 10,
      lambda.length = 22.02,
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = FALSE,
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: 'lambda.length' is not an integer"
  )
  expect_error(
    DTD_cv_lambda_cxx(
      lambda.seq = NULL,
      tweak.start = start.tweak,
      n.folds = 10,
      lambda.length = 22,
      train.data.list = training.data,
      F.GRAD.FUN = NA,
      EVAL.FUN = NA,
      cv.verbose = "TRUE",
      warm.start = FALSE
    ),
    "In DTD_cv_lambda: 'cv.verbose' must be a single value, either 'TRUE' or 'FALSE' (not 1 or 0)",
    fixed = TRUE
  )
})


context("Test desecent_generalized_fista ")
# most of the cases are tested via "standard test functions", therefore they
# will not be retested ...
# empty, i am going to fill this with test as soon as errors occur
# test_that("descent_generalized_fista ", {
#   }
# )

# the following include tests, but I am not going to test them here:
# nesterov_faktor
# nesterov_supspace
# nesterov_faktor
# gradient
# no test, as I think, if the arguments are invalid, they have been caught before:
# identiy, n2normed?


context("Test ggplot_gpath ")
test_that("ggplot_gpath ", {
  expect_error(
    ggplot_gpath(
      DTD.model = catch,
      number.pics = -1
      #               G.TRANSFORM.FUN = ,
      #               ITER.TRANSFORM.FUN = ,
      #               y.lab = , x.lab = ,
      #              subset = ,
      #              main = ,
      #              plot.legend =
    ),
    "In ggplot_gpath: 'number.pics' is below minimal value"
  )
  expect_error(
    ggplot_gpath(
      DTD.model = catch,
      #               number.pics = -1,
      G.TRANSFORM.FUN = as.character
      #               ITER.TRANSFORM.FUN = ,
      #               y.lab = , x.lab = ,
      #              subset = ,
      #              main = ,
      #              plot.legend =
    ),
    "In ggplot_gpath: 'G.TRANSFORM.FUN' does not return numeric vector."
  )
  expect_error(
    ggplot_gpath(
      DTD.model = catch,
      #               number.pics = -1,
      #              G.TRANSFORM.FUN = as.character
      ITER.TRANSFORM.FUN = as.character
      #               y.lab = , x.lab = ,
      #              subset = ,
      #              main = ,
      #              plot.legend =
    ),
    "In ggplot_gpath: 'ITER.TRANSFORM.FUN' does not return numeric vector."
  )
})

context("Test ggplot_convergence")
test_that("ggplot_convergence ", {
  expect_message(
    ggplot_convergence(
      DTD.model = catch,
      X.matrix = matrix(1:9, nrow = 3),
      test.data = test.data,
      estimate.c.type = "direct",
      title = ""
    ),
    "In ggplot_convergence: rownames('X.matrix') does not fit 'DTD.model' (rownames of History entry). Therefore convergence can not be shown on 'test.data'.",
    fixed = TRUE
  )
  expect_error(
    ggplot_convergence(DTD.model = NULL),
    "In ggplot_convergence: DTD.model is not a list"
  )
  expect_message(
    ggplot_convergence(
      DTD.model = catch,
      estimate.c.type = "direct",
      test.data = NA
    ),
    "In ggplot_convergence: 'test.data' must be provided as a list with two entries: 'quantities' and 'mixtures'. Therefore, 'test.data' can not be used. Plotting only training convergence",
    fixed = TRUE
  )
  expect_message(
    ggplot_convergence(
      DTD.model = catch,
      estimate.c.type = "direct",
      test.data = test.data,
      title = log2
    ),
    "In ggplot_convergence: 'title' can not be used 'as.character'"
  )
})

context("Test ggplot_true_vs_esti")
test_that("ggplot_true_vs_esti", {
  expect_error(
    ggplot_true_vs_esti(
      DTD.model = 1:100,
      X.matrix = NULL,
      test.data = test.data,
      estimate.c.type = "direct"
    ),
    "In ggplot_true_vs_esti: can't use 'X.matrix'. Neither via X.matrix argument nor included in the DTD.model."
  )
  expect_message(
    ggplot_true_vs_esti(
      DTD.model = catch,
      X.matrix = -1:10,
      test.data = test.data,
      estimate.c.type = "direct"
    ),
    "In ggplot_true_vs_esti: provided 'X.matrix' could not be used, therefore used: 'DTD.model$reference.X'",
    fixed = TRUE
  )
  expect_error(
    ggplot_true_vs_esti(
      DTD.model = catch,
      X.matrix = -1:10,
      test.data = NULL,
      estimate.c.type = "direct"
    ),
    "In ggplot_true_vs_esti: 'test.data' must be provided as a list with two entries: 'quantities' and 'mixtures'."
  )
  expect_error(
    ggplot_true_vs_esti(
      DTD.model = catch,
      X.matrix = NA,
      test.data = NULL,
      estimate.c.type = "direct"
    ),
    "In ggplot_true_vs_esti: 'test.data' must be provided as a list with two entries: 'quantities' and 'mixtures'."
  )
  expect_error(
    ggplot_true_vs_esti(
      DTD.model = catch,
      X.matrix = -1:10,
      test.data = NA,
      estimate.c.type = "direct"
    ),
    "In ggplot_true_vs_esti: 'test.data' must be provided as a list with two entries: 'quantities' and 'mixtures'."
  )
})
context("Test ggplot_ghistorgram ")
test_that("ggplot_ghistorgram ", {
  expect_error(
    ggplot_ghistogram(
      DTD.model = "NA",
      n.bins = 10,
      G.TRANSFORM.FUN = log2,
      title = "",
      x.lab = ""
    ),
    "In ggplot_ghistogram: 'DTD.model' does not fit"
  )
  expect_error(
    ggplot_ghistogram(
      DTD.model = catch,
      n.bins = -7,
      G.TRANSFORM.FUN = log2,
      title = "",
      x.lab = ""
    ),
    "In ggplot_ghistogram: 'n.bins' is below minimal value"
  )
  expect_error(
    ggplot_ghistogram(
      DTD.model = catch,
      n.bins = 2,
      G.TRANSFORM.FUN = "log2",
      title = "",
      x.lab = ""
    ),
    "In ggplot_ghistogram: 'G.TRANSFORM.FUN' is not a function"
  )
  # rest is checked via default functions
})

context("Test ggplot_gpath ")
test_that("ggplot_gpath ", {
  expect_error(
    ggplot_gpath(
      DTD.model = list(),
      number.pics = 3,
      G.TRANSFORM.FUN = identity,
      ITER.TRANSFORM.FUN = identity,
      y.lab = "",
      x.lab = "",
      subset = NA,
      title = "",
      show.legend = FALSE
    )
  )
  expect_error(
    ggplot_gpath(
      DTD.model = catch,
      number.pics = -3,
      G.TRANSFORM.FUN = identity,
      ITER.TRANSFORM.FUN = identity,
      y.lab = "",
      x.lab = "",
      subset = NA,
      title = "",
      show.legend = FALSE
    ),
    "In ggplot_gpath: 'number.pics' is below minimal value"
  )
  expect_error(
    ggplot_gpath(
      DTD.model = catch,
      number.pics = -3,
      G.TRANSFORM.FUN = identity,
      ITER.TRANSFORM.FUN = identity,
      y.lab = "",
      x.lab = "",
      subset = NA,
      title = "",
      show.legend = FALSE
    ),
    "In ggplot_gpath: 'number.pics' is below minimal value"
  )
  expect_error(
    ggplot_gpath(
      DTD.model = catch,
      number.pics = c(1, 2, 3, -3),
      G.TRANSFORM.FUN = identity,
      ITER.TRANSFORM.FUN = identity,
      y.lab = "",
      x.lab = "",
      subset = NA,
      title = "",
      show.legend = FALSE
    ),
    "In ggplot_gpath: 'number.pics' is not a single integer"
  )
  expect_error(
    ggplot_gpath(
      DTD.model = catch,
      number.pics = 4,
      G.TRANSFORM.FUN = "identity",
      ITER.TRANSFORM.FUN = "identity",
      y.lab = "",
      x.lab = "",
      subset = NA,
      title = "",
      show.legend = FALSE
    ),
    "In ggplot_gpath: 'G.TRANSFORM.FUN' does not return numeric vector"
  )
  expect_error(
    ggplot_gpath(
      DTD.model = catch,
      number.pics = 4,
      G.TRANSFORM.FUN = identity,
      ITER.TRANSFORM.FUN = "identity",
      y.lab = "",
      x.lab = "",
      subset = NA,
      title = "",
      show.legend = FALSE
    ),
    "In ggplot_gpath: 'ITER.TRANSFORM.FUN' does not return numeric vector"
  )
  expect_message(
    ggplot_gpath(
      DTD.model = catch,
      number.pics = 4,
      G.TRANSFORM.FUN = identity,
      ITER.TRANSFORM.FUN = identity,
      y.lab = identity,
      x.lab = "",
      subset = NA,
      title = "",
      show.legend = FALSE
    ),
    "In ggplot_gpath: 'y.lab' can not be used 'as.character'"
  )
  expect_message(
    ggplot_gpath(
      DTD.model = catch,
      number.pics = 4,
      G.TRANSFORM.FUN = identity,
      ITER.TRANSFORM.FUN = identity,
      y.lab = "",
      x.lab = "",
      subset = c(1, 2, 3),
      title = "",
      show.legend = FALSE
    ),
    "In ggplot_gpath: subset could not be used, therefore complete tweak, and history will be used"
  )

  expect_error(
    ggplot_gpath(
      DTD.model = catch,
      number.pics = 4,
      G.TRANSFORM.FUN = identity,
      ITER.TRANSFORM.FUN = identity,
      y.lab = "",
      x.lab = "",
      subset = c(1, 2, 3),
      title = "",
      show.legend = "FALSE"
    )
  )
})

context("Test ggplot_cv ")
test_that("ggplot_cv ", {
  catch <- train_deconvolution_model(
    tweak = start.tweak,
    X.matrix = X.matrix,
    train.data.list = training.data,
    test.data.list = test.data,
    lambda.seq = NULL,
    cv.verbose = FALSE,
    estimate.c.type = "direct"#, verbose = TRUE
  )

  expect_error(
    ggplot_cv(
      DTD.model = catch,
      title = "iwas",
      LAMBDA.TRANS.FUN = as.character(),
      x.lab = "1",
      y.lab = "2"
    ),
    "In ggplot_cv: provided 'LAMBDA.TRANS.FUN' does not return a numeric, if called with '2'"
  )
  expect_error(
    ggplot_cv(
      DTD.model = list(),
      title = "",
      LAMBDA.TRANS.FUN = log2,
      upper.x.axis.info = "geometric-mean",
      x.lab = "",
      y.lab = ""
    ),
    "In ggplot_cv: 'DTD.model' does not fit"
  )
  expect_error(
    ggplot_cv(
      DTD.model = catch,
      title = "",
      LAMBDA.TRANS.FUN = log2(),
      upper.x.axis.info = "geometric-mean",
      x.lab = "",
      y.lab = ""
    ),
    "In ggplot_cv: provided 'LAMBDA.TRANS.FUN' does not return a numeric, if called with '2'."
  )
  expect_error(
    ggplot_cv(
      DTD.model = catch,
      title = "",
      LAMBDA.TRANS.FUN = log2(),
      upper.x.axis.info = "some_not_set_thing",
      x.lab = "",
      y.lab = ""
    ),
    "In ggplot_cv: 'upper.x.axis.info' must either be 'non-zero' or 'geometric-mean'"
  )
})


context("Test ggplot_heatmap ")
test_that("ggplot_heatmap ", {
  expect_error(
    ggplot_heatmap(
      DTD.model = "catch",
      X.matrix = X.matrix,
      test.data = test.data,
      estimate.c.type = "direct",
      title = "",
      feature.subset = 1.0
    ),
    "In ggplot_heatmap: DTD.model is neither a list nor a numeric vector"
  )
  expect_error(
    ggplot_heatmap(
      DTD.model = list(),
      X.matrix = X.matrix,
      test.data = test.data,
      estimate.c.type = "direct",
      title = "",
      feature.subset = 1.0
    ),
    "In ggplot_heatmap: There is no Tweak entry in the 'DTD.model'"
  )
})

context("Comparing R and cpp implemenation: ")
test_that("r and cpp: ",
  {
  list.of.models <- list()
  for(implemenation in c("cpp", "R")){
    for(runs in 1:5){
      name <- paste0(implemenation, runs)
      set.seed(1)
      list.of.models[[name]] <- train_deconvolution_model(
        tweak = start.tweak,
        X.matrix = X.matrix,
        use.implementation = implemenation,
        train.data.list = training.data,
        test.data.list = test.data,
        estimate.c.type = "direct"
        , verbose = FALSE
        , cv.verbose = FALSE
        , maxit = 100
      )
    }
  }

  tmp <- matrix(
    unlist(
      lapply(list.of.models, function(model){return(model$best.model$Tweak)})
      ),
    ncol = length(list.of.models)
  )
  colnames(tmp) <- names(list.of.models)

  tmp.cpp <- tmp[, which(grepl(pattern = "cpp", colnames(tmp)))]
  tmp.R <- tmp[, which(grepl(pattern = "R", colnames(tmp)))]

  expect_equal(
    object = apply(
      tmp.cpp, 1, function(row){
        r <- range(row)
        return(r[1] - r[2])
        }
      ),
    expected = rep(0, nrow(tmp.cpp))
    )


  expect_equal(
    object = apply(
      tmp.R, 1, function(row){
        r <- range(row)
        return(r[1] - r[2])
      }
    ),
    expected = rep(0, nrow(tmp.R))
  )

  expect_equal(
    object = tmp.R[, 1]
    , expected = tmp.cpp[, 1]
  )
  cor.on.test <- unlist(
    lapply(
      list.of.models,
      function(model){
        DTD::evaluate_cor(
          new.data = test.data$mixtures,
          true.compositions = test.data$quantities,
          DTD.model = model,
          estimate.c.type = "decide.on.model"
        )
        }
      ),
    use.names = TRUE
  )

  cpp.cot <- cor.on.test[grepl(pattern = "cpp", names(cor.on.test))]
  R.cot <- cor.on.test[grepl(pattern = "R", names(cor.on.test))]
  # print(cpp.cot)
  # print(R.cot)
  expect_equal(
    expected = length(unique(cpp.cot)),
    object = 1
  )
  expect_equal(
    expected = length(unique(R.cot)),
    object = 1
  )
  names(cpp.cot) <- NULL
  names(R.cot) <- NULL
  expect_equal(
    expected = cpp.cot,
    object = R.cot,
  )
  }
)
