# clean environment
rm(list = ls()); gc()

# I'd like to start with 'training the g vector'. Therefore I need a lot of stuff ...
library(DTD)
library(tictoc)
# library(devtools)
# load_all(".")

# fast testing:
lambda.len <- 4
n.folds <- 2
###
number.types <- 5
n.features <- 200
n.per.type <- 20
n.per.mixtures <- 10
maxit <- 25
n.samples <- n.features
random.data <- generate_random_data(
  n.types = number.types,
  n.samples.per.type = n.per.type,
  n.features = n.features
)
normalized.data <- normalize_to_count(random.data)
indicator.vector <- gsub("^Cell[0-9]*\\.", "", colnames(normalized.data))
names(indicator.vector) <- colnames(normalized.data)
perc.of.all.cells <- 0.1
include.in.X <- paste0("Type", 1:5)
sample.X <- sample_random_X(
  included.in.X = include.in.X,
  pheno = indicator.vector,
  expr.data = normalized.data,
  percentage.of.all.cells = perc.of.all.cells
)
X.matrix <- sample.X$X.matrix
samples.to.remove <- sample.X$samples.to.remove
remaining.mat <- normalized.data[, -which(colnames(normalized.data) %in% samples.to.remove)]
train.samples <- sample(
  x = colnames(remaining.mat),
  size = ceiling(ncol(remaining.mat) / 2),
  replace = FALSE
)
test.samples <- colnames(remaining.mat)[which(!colnames(remaining.mat) %in% train.samples)]

train.mat <- remaining.mat[, train.samples]
test.mat <- remaining.mat[, test.samples]
indicator.train <- indicator.vector[names(indicator.vector) %in% colnames(train.mat)]
training.data <- mix_samples(
  expr.data = train.mat,
  pheno = indicator.train,
  included.in.X = include.in.X,
  n.samples = n.samples,
  n.per.mixture = n.per.mixtures,
  verbose = FALSE
)
indicator.test <- indicator.vector[names(indicator.vector) %in% colnames(test.mat)]
test.data <- mix_samples(
  expr.data = test.mat,
  pheno = indicator.test,
  included.in.X = include.in.X,
  n.samples = n.samples,
  n.per.mixture = n.per.mixtures,
  verbose = FALSE
)

start.tweak <- rep(1, nrow(X.matrix))
names(start.tweak) <- rownames(X.matrix)


# # ab hier c++
# library(microbenchmark)
#
# time.table <- microbenchmark(
#   "cpp" = train_deconvolution_model(
#     tweak = start.tweak,
#     X.matrix = X.matrix,
#     train.data.list = training.data,
#     test.data.list = test.data,
#     estimate.c.type = "direct",
#     maxit = maxit,
#     n.folds = n.folds,
#     lambda.len = lambda.len,
#     cv.verbose = FALSE,
#     verbose = FALSE,
#     useImplementation = "cpp"
#   ),
#   "R" = train_deconvolution_model(
#     tweak = start.tweak,
#     X.matrix = X.matrix,
#     train.data.list = training.data,
#     test.data.list = test.data,
#     estimate.c.type = "direct",
#     maxit = maxit,
#     n.folds = n.folds,
#     lambda.len = lambda.len,
#     cv.verbose = FALSE,
#     verbose = FALSE,
#     useImplementation = "R"
#   )
# )

# "non_negative", 
for(esti.type in c("direct")){
  tic(esti.type)
  model <- train_deconvolution_model(
    tweak = start.tweak,
    X.matrix = X.matrix,
    train.data.list = training.data,
    test.data.list = test.data,
    estimate.c.type = esti.type,
    maxit = maxit,
    n.folds = 10,
    lambda.len = 10,
    cv.verbose = TRUE,
    verbose = TRUE,
    NORM.FUN = "identity",
    use.implementation = "R"
  )
  toc()
}
