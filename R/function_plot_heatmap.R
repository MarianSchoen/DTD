#' clustered heatmap of diag(g) * X
#'
#' For a DTD.model the 'ggplot_heatmap' function visualizes
#' \deqn{diag(g) * X} on a subset of features as a clusterd heatmap.\cr
#' Feature subsetting can either be done by a vector of strings, that match the feature names of X.\cr
#' or via 'explained correlation':
#' In order to assess the importance of a feature in the deconvolution process, we can exclude it from a trained model,
#' and observe the change of correlaion on a test set. If the correlation e.g. decreases by 1 %,
#' the gene explains 1 % correlation within the deconvolution model.
#' The 'ggplot_heatmap' function iteratively excludes each feature from the trained model,
#' resulting in a ranking for the genes.
#' The clustering for the heatmp will be calculated using only these features.
#'
#' For an example see section "Explained correlation" in the package vignette `browseVignettes("DTD")`
#'
#' @param DTD.model list as returned by \code{\link{train_deconvolution_model}}, \code{\link{DTD_cv_lambda}},
#' or \code{\link{descent_generalized_fista}}, or a numeric vector, which will be used as 'g'.
#' @param X.matrix numeric matrix, with features/genes as rows, and cell types as column.
#' Each column of X.matrix is a reference expression profile.
#' Notice, if the model includes a 'reference.X' matrix, this parameter can be set to NA.
#' @param test.data numeric matrix with samples as columns, and features as rows.
#' In the deconvolution formula 'test.data' is denoated as Y.
#' If subsetting is done via feature names, 'test.data' will not be used.
#' @param estimate.c.type string, either "non_negative", or "direct". Indicates how the algorithm finds the solution of
#' \eqn{arg min_C ||diag(g)(Y - XC)||_2}. \cr
#' If estimate.c.type is set to "direct" there is no regularization (see \code{\link{estimate_c}}),\cr
#' if estimate.c.type is set to "non_negative" the estimates "C"
#' must not be negative (non-negative least squares) (see (see \code{\link{estimate_nn_c}}))
#' @param feature.subset numeric or a vector of strings. If it is a numeric, "subset" features will be picked from the
#' 'explained correlation' ranking (if 'feature.subset' < 1, this is the fraction of feature, if 'feature.subset' > 1
#' it is the total amount). If it is a vector of strings, these features will be used.
#' @param title string, additionally title (default "")
#' @param log2.expression logical, in the heatmap, should the values be log transformed?
#'
#' @return ggplot object
#' @export
#'
#' @import ggplot2
ggplot_heatmap <- function(DTD.model,
                           X.matrix = NA,
                           test.data=NULL,
                           estimate.c.type = "decide.on.model",
                           title = "",
                           feature.subset = 100,
                           log2.expression = TRUE){
  # safety check: DTD.model
  if (is.list(DTD.model)) {
    if ("best.model" %in% names(DTD.model)) {
      tweak <- DTD.model$best.model$Tweak
    } else {
      if (!"Tweak" %in% names(DTD.model)) {
        stop("In ggplot_heatmap: There is no Tweak entry in the 'DTD.model'")
      } else {
        tweak <- DTD.model$Tweak
      }
      if ("estimate.c.type" %in% names(DTD.model)){
        estimate.c.type <- DTD.model$estimate.c.type
      }
    }
  } else {
    if(is.numeric(DTD.model)){
      tweak <- DTD.model
    }else{
      stop("In ggplot_heatmap: DTD.model is neither a list nor a numeric vector")
    }
  }
  # end -> DTD.model
  # safety check: X.matrix
  if(!is.matrix(X.matrix) || any(is.na(X.matrix)) || !is.numeric(X.matrix)){
    if (is.list(DTD.model) && "reference.X" %in% names(DTD.model)) {
      message("In ggplot_heatmap: provided 'X.matrix' could not be used, therefore used: 'DTD.model$reference.X'")
      X.matrix <- DTD.model$reference.X
    }else{
      stop("In ggplot_heatmap: can't use 'X.matrix'. Neither via X.matrix argument nor included in the DTD.model.")
    }
  }
  # end -> X.matrix
  # safety check: log2.expression
  test <- test_logical(test.value = log2.expression,
                       output.info = c("ggplot_heatmap", "log2.expression"))
  # end -> log2.expression

  # safety check: test.data is moved into "expl.cor" part,
  # because if subset is a list of genes, I don't need it
  # estimate.c.type is moved as well

  # logical, if features has been set
  features.set <- FALSE

  # safety check: subset
  if(length(feature.subset) != 1 && all(is.character(feature.subset))){
    useable.subset <- intersect(feature.subset, rownames(X.matrix))
    if(length(useable.subset) == 0){
      stop("In ggplot_heatmap: 'feature.subset' is provided as vector of character. However, none of them can be found in rownames(X.matrix).")
    }else{
      features <- useable.subset
      features.set <- TRUE
    }
  }else{
    test <- test_numeric(
      feature.subset
      , output.info = c("ggplot_heatmap", "feature.subset")
      , min = 0
      , max = Inf)
    if(feature.subset <= 1){ # fraction is provided
      feature.subset <- round(nrow(X.matrix) * feature.subset)
    }
    feature.subset <- ceiling(min(feature.subset, nrow(X.matrix)))
  }
  # end -> subset
  # safety check: title
  useable.title <- try(as.character(title), silent = TRUE)
  if(any(grepl(x = useable.title, pattern = "Error"))){
    stop("In ggplot_cv: provided 'title' can not be used as.character.")
  }
  # end -> title
  # if there is no test.data, features can not be detected via explained correlation ...
  if(is.null(test.data)){
    # ... therefore, set features to the rownames (complete X)
    if(!features.set){ # If the input feature.subset is a vector of strings, "features" has already been set
      features <- rownames(X.matrix)
      features.set <- TRUE
    }
  }

  if(!features.set){
    # safety check: test.data
    if(is.list(test.data) && length(test.data) == 2){
      if(!all(c("quantities", "mixtures") %in%  names(test.data))){
        stop("In ggplot_heatmap: entries of 'test.data' must be named 'quantities' and 'mixtures'.")
      }else{
        if(!is.matrix(test.data$mixtures)){
          stop("In ggplot_heatmap: 'test.data$mixtures' is not a matrix.")
        }
        if(!is.matrix(test.data$quantities)){
          stop("In ggplot_heatmap: 'test.data$quantities' is not a matrix. Therefore, 'test.data' can not be used.")
        }
      }
    }else{
      stop("In ggplot_heatmap: 'test.data' must be provided as a list with two entries: 'quantities' and 'mixtures'.")
    }
    # end -> test.data
    # safety check: estimate.c.tye
    test <- test_c_type(test.value = estimate.c.type,
                        output.info = c("ggplot_heatmap", "estimate.c.type"))
    # end -> estimate.c.type
    DTD.evCor.wrapper <- function(tweak,
                                  X = X.matrix,
                                  train.list = test.data,
                                  esti.c.type = estimate.c.type) {
      Y <- train.list$mixtures
      C <- train.list$quantities
      loss <- evaluate_cor(
        X.matrix = X,
        new.data = Y,
        true.compositions = C,
        DTD.model = tweak,
        estimate.c.type = esti.c.type
      ) / ncol(X)
      return(loss)
    }

    # correlation achieved with all features:
    maximum.cor <- -DTD.evCor.wrapper(tweak = tweak)

    manipulated.cor <- rep(NA, length(tweak))
    names(manipulated.cor) <- names(tweak)
    for(l.feature in names(tweak)){
      l.tweak <- tweak
      l.tweak[l.feature] <- 0
      # manipulated correlation works as follows:
      # if maximum.cor is e.g. 0.7,
      # the output of 'DTD.evCor.wrapper' is the loss, which is the negative correlation.
      # If it is e.g. -0.68, this results in a manipulated correlation of '+0.02'.
      # meaning, this gene explaines 2% correlation.
      # If the manipulated correlation would be negative, this means excluding this gene is favorable for the model.
      manipulated.cor[l.feature] <- (maximum.cor + DTD.evCor.wrapper(tweak = l.tweak))*100
    }
    sorted.manipulated.cor <- sort(manipulated.cor, decreasing = TRUE)
    features <- names(sorted.manipulated.cor)[1:feature.subset]
 #   upper.axis <- sorted.manipulated.cor[1:feature.subset]
  }
  tmp.rownames <- rownames(X.matrix)
  tmp.colnames <- colnames(X.matrix)

  X.times.g <- diag(tweak) %*% X.matrix
  colnames(X.times.g) <- tmp.colnames
  rownames(X.times.g) <- tmp.rownames
  if(log2.expression){
    tmp.X.matrix <- as.matrix(log2(X.times.g[features, ] + 1))
  }else{
    tmp.X.matrix <- X.times.g[features, ]
  }

  cell.type.cluster <- stats::hclust(stats::dist(x = t(tmp.X.matrix),
                                   method = "euclidean"),
                              method = "average")
  feature.cluster <- stats::hclust(stats::dist(x = tmp.X.matrix,
                                   method = "euclidean"),
                              method = "average")

  melted.X <- reshape2::melt(tmp.X.matrix,
                             value.name = "expression")
  melted.X$Var1 <- factor(
    x = melted.X$Var1,
    levels = rownames(tmp.X.matrix)[feature.cluster$order],
    ordered = TRUE
  )

  melted.X$Var2 <- factor(
    x = melted.X$Var2,
    levels = colnames(tmp.X.matrix)[cell.type.cluster$order],
    ordered = TRUE
  )

  if(log2.expression){
    legend.name <- "log2(expr)"
  }else{
    legend.name <- "expr"
  }

  pic0 <- ggplot(melted.X, aes(x = .data$Var1, y = .data$Var2)) +
    geom_tile(aes(fill = expression)) +
    scale_fill_gradient(
      name = legend.name,
      low = "darkblue",
      high = "yellow"
    ) +
    xlab("features") +
    ylab("Cell types") +
    ggtitle(title) +
    theme(
      axis.text.x = element_text(
        angle = 90
        , hjust = 1)
    )


  return(pic0)
}
