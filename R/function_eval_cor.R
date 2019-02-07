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
#'
#'
#' indicator.train <- indicator.list[names(indicator.list) %in% colnames(train.mat)]
#' training.data <- mix_samples(gene.mat = remaining.mat,
#'                              pheno = indicator.train,
#'                              included.in.X = include.in.X,
#'                              n.samples = 500,
#'                              n.per.mixture = 100,
#'                              verbose = FALSE)
#'
#'  start.tweak <- rep(1, nrow(X.matrix))
#'  sum.cor <- evaluate_cor(X = X.matrix,
#'                       new.data = training.data$mixtures,
#'                       true.compositions = training.data$quantities,
#'                       DTD.model = start.tweak)
#'
#'  rel.cor <- sum.cor/ncol(X.matrix)
#' cat("Relative correlation: ", -rel.cor, "\n")
evaluate_cor <- function(X.matrix = NA,
                         new.data,
                         true.compositions,
                         DTD.model){

  # the reference matrix can either be included in the DTD.model, or has to be past
  # via the X.matrix argument:
  if(any(is.na(X.matrix)) && is.list(DTD.model) && "reference.X" %in% names(DTD.model)){
    X <- DTD.model$reference.X
  }else{
    X <- X.matrix
  }
  # as a DTD.model either a list, or only the tweak vector can be used:
  if(is.list(DTD.model)){
    if("best.model" %in% names(DTD.model)){
      fista.output <- DTD.model$best.model
    }else{
      if("Tweak" %in% names(DTD.model)){
        stop("evaluate_cor: DTD.model does not fit")
      }else{
        fista.output <- DTD.model
      }
    }
  }else{
    tweak <- DTD.model
  }

  Y <- new.data
  C <- true.compositions
  # estimate C using reference matrix, bulks and the tweak vector:
  esti.cs <- estimate_c(X, Y, tweak)

  if(any(dim(esti.cs) != dim(C))){
    stop("evaluate_cor: dimension of estimated C do not fit")
  }

  # initialise loss:
  loss <- 0
  for(l1 in 1:nrow(C)){
    # calculate the correlation per Type and add them up
    loss <- loss + stats::cor(C[l1, ], esti.cs[l1, ])
  }

  # Our loss function is set to be a minimization problem, therefore we invert the sign:
  loss <- -loss
  return(loss)
}

