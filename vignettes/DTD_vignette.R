## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----include=FALSE-------------------------------------------------------
  # for fast testing 
  maxit <- 2
  nSamples <- 5

## ------------------------------------------------------------------------
  library(DTD)

## ------------------------------------------------------------------------
  number.types <- 5
  number.samples.per.type <- 10
  random.data <- generate.random.data(nTypes = number.types, 
                                      nSamples.perType = number.samples.per.type, 
                                      nFeatures = 100)
  print(random.data[1:5, 1:5])
  
  # for further exemplary visualization 
  example.index <- c(1, 12, 23, 34, 44)

## ------------------------------------------------------------------------
  print(colnames(random.data)[example.index])

## ------------------------------------------------------------------------
  if(max(random.data) < 25){ 
    random.data <- 2^random.data
  }

## ------------------------------------------------------------------------
  # In random.data the number of counts differs over all samples: 
  apply(random.data, 2, sum)[example.index]
  normalized.data <- normalizeToCount(random.data)
  # In normalized.data all samples share the same number of counts: 
  apply(normalized.data, 2, sum)[example.index]

## ------------------------------------------------------------------------
 indicator.list <- gsub("^Cell([0-9])*.", "", colnames(normalized.data))
 names(indicator.list) <- colnames(normalized.data)
 print(indicator.list[example.index])

## ------------------------------------------------------------------------
  include.in.X <- c("Type2", "Type3", "Type4", "Type5")

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
    average <- rowSums(normalized.data[, chosen.for.X, drop = F])
    X.matrix[, l.type] <- average
 }

## ------------------------------------------------------------------------
  train.mat <- random.data[, -which(colnames(random.data) %in% samples.to.remove)]

## ------------------------------------------------------------------------
  training.data <- mix.samples(gene.mat = train.mat,
                               pheno = indicator.list,
                               included.in.X = include.in.X, 
                               nSamples = nSamples, 
                               nPerMixture = 5, 
                               verbose = F)
  str(training.data)

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
  DTD.grad.wrapper <- function(tweak){
      X <- X.matrix
      Y <- training.data$mixtures
      C <- training.data$quantities
      grad <- Trace.H.gradient(X = X, Y = Y, C = C, gamma.vec = tweak)
      return(grad)
  }
  DTD.evCor.wrapper <- function(tweak){
      X <- X.matrix
      Y <- training.data$mixtures
      C <- training.data$quantities
      loss <- evaluate_cor(X = X, Y = Y, C = C, tweak = tweak)
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

