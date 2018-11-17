## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.align = "center"
)

## ----include=FALSE-------------------------------------------------------
  # for fast testing ...
  maxit <- 75
  nSamples <- 500
  do.cv <- TRUE

## ----include=FALSE-------------------------------------------------------
  # I'd like to start with 'training the g vector'. Therefore I need a lot of stuff ...
  library(DTD)
  number.types <- 20
  n.features <- 250 
  n.per.type <- 200
  nPerMixture <- 100
  random.data <- generate.random.data(nTypes = number.types, 
                                      nSamples.perType = n.per.type, 
                                      nFeatures = n.features)
  normalized.data <- normalizeToCount(random.data)
  indicator.list <- gsub("^Cell[0-9]*\\.", "", colnames(normalized.data))
  names(indicator.list) <- colnames(normalized.data)
  perc.of.all.cells <- 0.1
  include.in.X <- paste0("Type", 2:5)
  sample.X <- sample.random.X(included.in.X = include.in.X, 
                              pheno = indicator.list, 
                              exp.data = normalized.data, 
                              percentage.of.all.cells = perc.of.all.cells)
  X.matrix <- sample.X$X.matrix
  samples.to.remove <- sample.X$samples.to.remove
  remaining.mat <- normalized.data[, -which(colnames(normalized.data) %in% samples.to.remove)]
  train.samples <- sample(x = colnames(remaining.mat), 
                          size = ceiling(ncol(remaining.mat)/2), 
                          replace = FALSE)
  test.samples <- colnames(remaining.mat)[which(!colnames(remaining.mat) %in% train.samples)]
  
  train.mat <- remaining.mat[, train.samples]
  test.mat <- remaining.mat[, test.samples]
  indicator.train <- indicator.list[names(indicator.list) %in% colnames(train.mat)]
  training.data <- mix.samples(gene.mat = train.mat,
                               pheno = indicator.train,
                               included.in.X = include.in.X, 
                               nSamples = nSamples, 
                               nPerMixture = nPerMixture, 
                               verbose = F)
  indicator.test <- indicator.list[names(indicator.list) %in% colnames(test.mat)]
  test.data <- mix.samples(gene.mat = test.mat,
                           pheno = indicator.test,
                           included.in.X = include.in.X, 
                           nSamples = nSamples, 
                           nPerMixture = nPerMixture, 
                           verbose = F)
  # wrapper for gradient:
   DTD.grad.wrapper <- function(tweak, 
                                X = X.matrix,
                                train.list = training.data){
      Y <- train.list$mixtures
      C <- train.list$quantities
      grad <- Trace.H.gradient(X = X, Y = Y, C = C, tweak = tweak)
      grad[grad > 0] <- 0
      return(grad)
   }
  # wrapper for evaluation:
  DTD.evCor.wrapper <- function(tweak, 
                                X = X.matrix, 
                                train.list = training.data){
      Y <- train.list$mixtures
      C <- train.list$quantities
      loss <- evaluate_cor(X = X, Y = Y, C = C, tweak = tweak)/ncol(X.matrix)
      return(loss)
  }
  start_tweak <- rep(1, nrow(X.matrix)); names(start_tweak) <- rownames(X.matrix)
  start.loss <- DTD.evCor.wrapper(start_tweak)
  DTD.evCor.wrapper.test <- function(tweak, 
                              X = X.matrix, 
                              yc.list = test.data){
    Y <- yc.list$mixtures
    C <- yc.list$quantities
    loss <- evaluate_cor(X = X, Y = Y, C = C, tweak = tweak)/ncol(X.matrix)
    return(loss)
  }

## ---- fig.width = 7, fig.align="center"----------------------------------
  catch <- descent_generalized_fista(tweak_vec = start_tweak,
                                   F.GRAD.FUN = DTD.grad.wrapper,
                                   ST.FUN = soft_thresholding,
                                   FACTOR.FUN = nesterov_faktor,
                                   EVAL.FUN = DTD.evCor.wrapper,
                                   NORM.FUN = n2normed, 
                                   line_search_speed = 2,
                                   maxit = maxit,
                                   save_all_tweaks = TRUE, 
                                   verbose = F)
  str(catch)
  
  print(ggplot_convergence(fista.output = catch, 
                         EVAL.FUN = DTD.evCor.wrapper.test,
                         main = "DTD Vignette"))

## ------------------------------------------------------------------------
  library(DTD)

## ----echo=FALSE, results = "asis"----------------------------------------
  cat("```\n")
  cat(" number.types <- " , number.types, "\n",
      "n.features <- ", n.features, "\n",
      "n.per.type <- ", n.per.type, "\n")
  cat("```\n")

## ------------------------------------------------------------------------
  random.data <- generate.random.data(nTypes = number.types, 
                                      nSamples.perType = n.per.type, 
                                      nFeatures = n.features)

## ------------------------------------------------------------------------
  print(random.data[1:5, 1:5])

## ----include=FALSE-------------------------------------------------------
  # for further exemplary visualization 
  example.index <- c(1, n.per.type + 1, 2*n.per.type + 1, 3 * n.per.type+1)

## ------------------------------------------------------------------------
  print(colnames(random.data)[example.index])

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

## ----echo=FALSE, results = "asis"----------------------------------------
  cat("```\n")
  cat(" include.in.X <- " , paste0("c(\"", paste(include.in.X ,collapse = "\", \""), "\")"), "\n")
  cat(" perc.of.all.cells <- ", perc.of.all.cells,  "\n")
  cat("```\n")

## ------------------------------------------------------------------------
  sample.X <- sample.random.X(included.in.X = include.in.X, 
                              pheno = indicator.list, 
                              exp.data = normalized.data, 
                              percentage.of.all.cells = perc.of.all.cells)
  X.matrix <- sample.X$X.matrix
  samples.to.remove <- sample.X$samples.to.remove

## ------------------------------------------------------------------------
  # removing samples that have been used in X.matrix
  remaining.mat <- normalized.data[, -which(colnames(normalized.data) %in% samples.to.remove)]

  # sampling training samples: (notice, that train test seperation is 50:50)
  train.samples <- sample(x = colnames(remaining.mat), 
                          size = ceiling(ncol(remaining.mat)/2), 
                          replace = FALSE)
  # selecting test samples: 
  test.samples <- colnames(remaining.mat)[which(!colnames(remaining.mat) %in% train.samples)]
  
  # extract data matrices for training and testing: 
  train.mat <- remaining.mat[, train.samples]
  test.mat <- remaining.mat[, test.samples]

## ----echo=FALSE, results = "asis"----------------------------------------
  cat("```\n")
  cat(" nSamples <- ", nSamples,  "\n",
      "nPerMixture <- ", nPerMixture,  "\n")
  cat("```\n")

## ------------------------------------------------------------------------
  indicator.train <- indicator.list[names(indicator.list) %in% colnames(train.mat)]
  training.data <- mix.samples(gene.mat = train.mat,
                               pheno = indicator.train,
                               included.in.X = include.in.X, 
                               nSamples = nSamples, 
                               nPerMixture = nPerMixture, 
                               verbose = F)
  str(training.data)

## ------------------------------------------------------------------------
  indicator.test <- indicator.list[names(indicator.list) %in% colnames(test.mat)]
  test.data <- mix.samples(gene.mat = test.mat,
                           pheno = indicator.test,
                           included.in.X = include.in.X, 
                           nSamples = nSamples, 
                           nPerMixture = nPerMixture, 
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
   DTD.grad.wrapper <- function(tweak, 
                                X = X.matrix, 
                                train.list = training.data){
      Y <- train.list$mixtures 
      C <- train.list$quantities
      grad <- Trace.H.gradient(X = X, Y = Y, C = C, tweak = tweak)
      return(grad)
   }
  # wrapper for evaluation:
  DTD.evCor.wrapper <- function(tweak, 
                                X = X.matrix, 
                                train.list = training.data){
      Y <- train.list$mixtures 
      C <- train.list$quantities
      loss <- evaluate_cor(X = X, Y = Y, C = C, tweak = tweak)/ncol(X)
      return(loss)
  }

## ------------------------------------------------------------------------
  start_tweak <- rep(1, nrow(X.matrix))
  names(start_tweak) <- rownames(X.matrix)
  start.loss <- DTD.evCor.wrapper(start_tweak)
  cat("Start average correlation: ", -start.loss/ncol(X.matrix), "\n")
  catch <- descent_generalized_fista(tweak_vec = start_tweak,
                                     F.GRAD.FUN = DTD.grad.wrapper,
                                     EVAL.FUN = DTD.evCor.wrapper,
                                     maxit = maxit,
                                     save_all_tweaks = T, 
                                     verbose = F)
  str(catch)

## ---- fig.width=5--------------------------------------------------------
  DTD.evCor.wrapper.test <- function(tweak, 
                                X = X.matrix, 
                                test.list = test.data){
      Y = test.list$mixtures 
      C = test.list$quantities
      loss <- evaluate_cor(X = X, Y = Y, C = C, tweak = tweak)/ncol(X)
      return(loss)
  }
  print(ggplot_convergence(fista.output = catch, 
                           EVAL.FUN = DTD.evCor.wrapper.test,
                           main = "DTD Vignette"))

## ---- fig.width=7--------------------------------------------------------
  test.estimations <- est.cs(X = X.matrix, 
                             Y = test.data$mixtures, 
                             gamma.vec = catch$Tweak
                            )
 print(ggplot_true_vs_esti(estimatedC = test.estimations,
                           trueC = test.data$quantities, 
                           norm.columnwise = FALSE)
       )

## ------------------------------------------------------------------------
    print(ggplot_ghistogram(fista.output = catch))

## ---- fig.show="hold", fig.width=3---------------------------------------
  singlePic <- ggplot_gPath(fista.output = catch, 
                            number.pics = 1, 
                            main="All genes")
   print(singlePic$gPath)

## ------------------------------------------------------------------------
  singlePic.subset <- ggplot_gPath(fista.output = catch,
                                   number.pics = 1, 
                                   main = "Gene Subset", 
                                   subset = c("gene121", "gene66"), 
                                   plot_legend = TRUE)
   print(singlePic.subset$gPath)

## ---- fig.width=5--------------------------------------------------------
  multPic <- ggplot_gPath(fista.output = catch, 
                          number.pics = 3,  
                          main = "All Genes, split into 3")
  print(multPic$gPath)

## ------------------------------------------------------------------------
  if(do.cv){
  sequence <- 0.001*2^seq(5, -5, length.out = 10)
  set.seed(2018)
  cv.object <- DTD_cv_lambda(tweak.start = start_tweak, 
                             nfolds = 5, 
                             lambda.seq = sequence, 
                             cv.verbose = TRUE, 
                             train.list = training.data, 
                             F.GRAD.FUN = DTD.grad.wrapper, 
                             EVAL.FUN = DTD.evCor.wrapper, 
                             ST.FUN = soft_thresholding,
                             FACTOR.FUN = nesterov_faktor,
                             line_search_speed = 2,
                             maxit = maxit,
                             save_all_tweaks = FALSE, 
                             NORM.FUN = n2normed, 
                             use_restarts = TRUE,
                             verbose = FALSE
  )
  str(cv.object)
  }

## ---- fig.align='center'-------------------------------------------------
 if(do.cv) print(ggplot_cv(cv.object$cv.obj))

