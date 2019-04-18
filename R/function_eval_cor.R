#' Evaluate correlation
#'
#' The loss-function learning digital tissue deconvolution finds a vector g which minimizes the Loss function L\cr
#' \deqn{L(g) = - \sum cor(true_C, estimatd_C(g))}
#' The evaluate_cor function returns the value of the Loss function.
#'
#' @param X.matrix numeric matrix with cells as columns, and features as rows.
#'  Reference matrix X of the DTD problem. X.matrix can be set to NA (default), if the DTD.model
#'  includes the reference matrix X (default for \code{\link{train_correlatio_model}})
#' @param new.data numeric matrix with samples as columns, and features as rows.
#' In the formula above denoated as Y.
#' @param DTD.model either a numeric vector with length of nrow(X),
#' or a list returned by \code{\link{train_correlatio_model}}, \code{\link{DTD_cv_lambda}},
#' or\code{\link{descent_generalized_fista}}. In the equation above
#'   the DTD.model provides the vector g.
#' @param true.compositions numeric matrix with cells as rows, and mixtures as columns.
#' Each row of C holds the distribution of the cell over all mixtures.
#' @param estimate.c.type string, either "non_negative", or "direct". Indicates how the algorithm finds the solution of
#' \eqn{arg min_C ||diag(g)(Y - XC)||_2}. If estimate.c.type is set to "direct" there is no regularization
#' (see \code{\link{estimate_c}}),
#' if estimate.c.type is set to "non_negative" the estimates "C" must not be negative (non-negative least squares) (see \code{\link{estimate_nn_c}})
#'
#' @return float, value of the Loss function
#'
#' @export
#' @examples
#' library(DTD)
#' random.data <- generate_random_data(n.types = 10,
#'                                     n.samples.per.type = 150,
#'                                     n.features = 250,
#'                                     sample.type = "Cell",
#'                                     feature.type = "gene")
#'
#' # normalize all samples to the same amount of counts:
#' normalized.data <- normalize_to_count(random.data)
#'
#' # extract indicator list.
#' # This list contains the Type of the sample as value, and the sample name as name
#' indicator.list <- gsub("^Cell[0-9]*\\.", "", colnames(random.data))
#' names(indicator.list) <- colnames(random.data)
#'
#' # extract reference matrix X
#' # First, decide which cells should be deconvoluted.
#' # Notice, in the mixtures there can be more cells than in the reference matrix.
#' include.in.X <- paste0("Type", 2:7)
#'
#' percentage.of.all.cells <- 0.2
#' sample.X <- sample_random_X(included.in.X = include.in.X,
#'                             pheno = indicator.list,
#'                             exp.data = normalized.data,
#'                             percentage.of.all.cells = percentage.of.all.cells)
#' X.matrix <- sample.X$X.matrix
#' samples.to.remove <- sample.X$samples.to.remove
#' remaining.mat <- normalized.data[, -which(colnames(normalized.data) %in% samples.to.remove)]
#'
#' indicator.remain <- indicator.list[names(indicator.list) %in% colnames(remaining.mat)]
#' training.data <- mix_samples(exp.data = remaining.mat,
#'                              pheno = indicator.remain,
#'                              included.in.X = include.in.X,
#'                              n.samples = 500,
#'                              n.per.mixture = 100,
#'                              verbose = FALSE)
#'
#'  start.tweak <- rep(1, nrow(X.matrix))
#'  sum.cor <- evaluate_cor(X.matrix = X.matrix,
#'                       new.data = training.data$mixtures,
#'                       true.compositions = training.data$quantities,
#'                       DTD.model = start.tweak,
#'                       estimate.c.type = "direct")
#'
#'  rel.cor <- sum.cor/ncol(X.matrix)
#' cat("Relative correlation: ", -rel.cor, "\n")
evaluate_cor <- function(X.matrix = NA,
                         new.data,
                         true.compositions,
                         DTD.model,
                         estimate.c.type){


  # safety check: estimate.c.type
  ESTIMATE.C.FUN <- test_c_type(test.value = estimate.c.type,
                      output.info = c("evaluate_cor", "estimate.c.type"))
  # end -> estimate.c.type

  # the reference matrix can either be included in the DTD.model, or has to be past
  # via the X.matrix argument:
  if(is.matrix(X.matrix) && !any(is.na(X.matrix))){
    X <- X.matrix
  }else{
    if(is.list(DTD.model) && "reference.X" %in% names(DTD.model)){
      if("reference.X" %in% names(DTD.model)){
        X <- DTD.model$reference.X
      }
      if("best.model" %in% names(DTD.model)){
        if("reference.X" %in% names(DTD.model$best.model)){
          X <- DTD.model$best.model$reference.X
        }
      }
    }
  }
  if(!exists("X")){
    stop("In evaluate_cor: 'X.matrix' must be provided either as the 'X.matrix' argument, or within the DTD.model")
  }

  # as a DTD.model either a list, or only the tweak vector can be used:
  if (is.list(DTD.model)) {
    if ("best.model" %in% names(DTD.model)) {
      if("Tweak" %in% names(DTD.model$best.model)){
        gamma.vec <- DTD.model$best.model$Tweak
      }
    }else {
      if ("Tweak" %in% names(DTD.model)) {
        gamma.vec <- DTD.model$Tweak
      }else {
        stop("In evaluate_cor: DTD.model does not fit")
      }
    }
  }else {
    gamma.vec <- DTD.model
  }

  if(length(gamma.vec) != nrow(X)){
    stop("In evaluate_cor: 'DTD.model' does not fit 'X.matrix'. Check if the provided model includes as many 'Tweak' entries as there are features (=> rows) in 'X.matrix'")
  }

  if(!is.null(names(gamma.vec))){
    if(all(rownames(X) %in% names(gamma.vec))){
      gamma.vec <- gamma.vec[rownames(X)]
    }else{
      stop("In evaluate_cor: There are features within 'X.matrix' where no tweak/g entry in the 'DTD.model' can be found")
    }
  }

  if(!is.numeric(gamma.vec)){
    stop("In evaluate_cor: Used 'Tweak' (from 'DTD.model') is not numeric")
  }

  if(nrow(true.compositions) != ncol(X)){
    stop("In evaluate_cor: 'nrow(true.composition)' does not match 'ncol(X.matrix)'")
  }
  if(ncol(true.compositions) != ncol(new.data)){
    stop("In evaluate_cor: 'ncol(true.composition)' does not match 'ncol(new.data)' ")
  }

  if(all(colnames(true.compositions) %in% colnames(new.data))){
    true.compositions <- true.compositions[, colnames(new.data)]
  }else{
    stop("In evaluate_cor: 'colnames(true.compositions)' do not match 'colnames(new.data)'")
  }

  if(all(rownames(true.compositions) %in% colnames(X))){
    true.compositions <- true.compositions[colnames(X), ]
  }else{
    stop("In evaluate_cor: 'rownames(true.compositions)' do not match 'colnames(X.matrix)'")
  }

  if(all(rownames(new.data) %in% rownames(X))){
    new.data <- new.data[rownames(X), ]
  }else{
    stop("in evaluate_cor: not all'rownames(new.data)' are in 'rownames(X.matrix)'")
  }

  Y <- new.data
  C <- true.compositions
  # estimate C using reference matrix, bulks and the tweak vector:
  esti.cs <- ESTIMATE.C.FUN(X, Y, gamma.vec)

  if(any(dim(esti.cs) != dim(C))){
    stop("evaluate_cor: dimension of estimated C do not fit")
  }

  # initialise loss:
  loss <- 0
  for(l1 in 1:nrow(C)){
    if(sd(esti.cs[l1, ]) != 0){
      # calculate the correlation per Type and add them up
      loss <- loss + stats::cor(C[l1, ], esti.cs[l1, ])
    }
  }

  # Our loss function is set to be a minimization problem, therefore we invert the sign:
  loss <- -loss
  return(loss)
}

