## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----include=FALSE-------------------------------------------------------
  # for fast testing 
  maxit <- 250
  nSamples <- 500

## ------------------------------------------------------------------------
  library(DTD)

## ------------------------------------------------------------------------
  number.types <- 10
  random.data <- generate.random.data(nTypes = number.types, 
                                      nSamples.perType = 150, 
                                      nFeatures = 250)

## ------------------------------------------------------------------------
  print(random.data[1:5, 1:5])
  
  # for further exemplary visualization 
  example.index <- c(1, 151, 301, 451)

## ------------------------------------------------------------------------
  print(colnames(random.data)[example.index])

## ------------------------------------------------------------------------
    random.data <- 2^random.data - 1

## ------------------------------------------------------------------------
  # In random.data the number of counts differs over all samples: 
  apply(random.data, 2, sum)[example.index]
  normalized.data <- normalizeToCount(random.data)
  # In normalized.data all samples share the same number of counts: 
  apply(normalized.data, 2, sum)[example.index]

## ------------------------------------------------------------------------
 indicator.list <- gsub("^Cell[0-9]*\\.", "", colnames(normalized.data))
 names(indicator.list) <- colnames(normalized.data)
 print(indicator.list[example.index])

## ------------------------------------------------------------------------
  include.in.X <- paste0("Type", 2:7)
  print(include.in.X)

## ------------------------------------------------------------------------
  X.matrix <- matrix(NA, nrow=nrow(normalized.data), ncol=length(include.in.X))
  colnames(X.matrix) <- include.in.X
  rownames(X.matrix) <- rownames(normalized.data)

## ------------------------------------------------------------------------
  percentage.of.all.cells <- 0.2
  samples.to.remove <- c()
  for(l.type in include.in.X){
    # get sample names of all cells of type "l.type" 
    all.of.type <- names(indicator.list)[which(indicator.list == l.type)]
    
    # randomly sample some cells
    chosen.for.X <- sample(x = all.of.type,
                           size = ceiling(length(all.of.type) * percentage.of.all.cells),
                           replace = FALSE)
    
    # Add those cells which will be included in X to the list of samples.to.remove 
    samples.to.remove <- c(samples.to.remove, chosen.for.X)

    # for each gene average over the selected 
    average <- rowSums(normalized.data[, chosen.for.X, drop = FALSE])
    X.matrix[, l.type] <- average
 }

## ------------------------------------------------------------------------
  remaining.mat <- normalized.data[, -which(colnames(normalized.data) %in% samples.to.remove)]
  train.samples <- sample(x = colnames(remaining.mat), 
                          size = ceiling(ncol(remaining.mat)/2), 
                          replace = FALSE)
  test.samples <- colnames(remaining.mat)[which(!colnames(remaining.mat) %in% train.samples)]
  
  train.mat <- remaining.mat[, train.samples]
  test.mat <- remaining.mat[, test.samples]

## ------------------------------------------------------------------------
  indicator.train <- indicator.list[names(indicator.list) %in% colnames(train.mat)]
  training.data <- mix.samples(gene.mat = train.mat,
                               pheno = indicator.train,
                               included.in.X = include.in.X, 
                               nSamples = nSamples, 
                               nPerMixture = 100, 
                               verbose = F)
  str(training.data)

## ------------------------------------------------------------------------
  indicator.test <- indicator.list[names(indicator.list) %in% colnames(test.mat)]
  test.data <- mix.samples(gene.mat = test.mat,
                           pheno = indicator.test,
                           included.in.X = include.in.X, 
                           nSamples = nSamples, 
                           nPerMixture = 100, 
                           verbose = F)

## ------------------------------------------------------------------------
  # Here, we set "Type1" to be special:
  special.samples <- c("Type1")
  # and all other to be normal: 
  all.samples <- unique(indicator.list)
  sample.names <- all.samples[- which(all.samples %in% special.samples)]
  
  # reduce indicator list to those samples included in training: 
  indicator.list <- indicator.list[names(indicator.list) %in% colnames(train.mat)]
  
  training.data.jitter <- mix.samples.jitter(sample.names = sample.names,
                                             special.samples = special.samples, 
                                             nSamples = nSamples, 
                                             datamatrix = train.mat, 
                                             pheno = indicator.list, 
                                             verbose = F, 
                                             add_jitter = T, 
                                             included.in.X = include.in.X)

## ------------------------------------------------------------------------
   # wrapper for gradient:
   DTD.grad.wrapper <- function(tweak){
      X <- X.matrix
      Y <- training.data$mixtures
      C <- training.data$quantities
      grad <- Trace.H.gradient(X = X, Y = Y, C = C, tweak = tweak)
      return(grad)
   }
  # wrapper for evaluation:
  DTD.evCor.wrapper <- function(tweak){
      X <- X.matrix
      Y <- training.data$mixtures
      C <- training.data$quantities
      loss <- evaluate_cor(X = X, Y = Y, C = C, tweak = tweak)
      return(loss)
  }

## ------------------------------------------------------------------------
  start_tweak <- rep(1, nrow(X.matrix))
  catch <- descent_generalized_fista(tweak_vec = start_tweak,
                                     F.GRAD.FUN = DTD.grad.wrapper,
                                     ST.FUN = soft_thresholding,
                                     FACTOR.FUN = nesterov_faktor,
                                     EVAL.FUN = DTD.evCor.wrapper,
                                     line_search_speed = 2,
                                     maxit = maxit,
                                     save_all_tweaks = T, 
                                     verbose = F)
  str(catch)

## ---- fig.width=5--------------------------------------------------------
  print(ggplot_correlation(fista.output = catch, 
                           test.set = test.data, 
                           X.matrix = X.matrix, 
                           main = "DTD Vignette"))

