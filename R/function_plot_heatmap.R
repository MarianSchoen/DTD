#' Heatmap of diag(g) * X
#'
#' Plot a heatmap of the reference matrix X times , using a selected subset of features
#' NOT FINISHED, STOP
#' Top features are selected via "explained correlation":
#' In order to assess the importance of a feature in the deconvolution process, we can exclude it from a trained model,
#' and observe the change of correlaion on a test set. If the correlation e.g. decreases by 1 %,
#' the gene explains 1 % correlation within the deconvolution model.
#' The 'ggplot_heatmap' function iteratively excludes each feature from the trained model,
#' resulting in a ranking for the genes.
#' The clustering for the heatmp will be calculated using only these features.
#'
#' For an example see section "Explained correlation" in the package vignette `browseVignettes("DTD")`
#'
#' @param DTD.model list as returned by \code{\link{train_correlatio_model}}, \code{\link{DTD_cv_lambda}},
#' or \code{\link{descent_generalized_fista}}
#' @param X.matrix numeric matrix, with features/genes as rows, and cell types as column.
#' Each column of X.matrix is a reference expression profile
#' @param test.data numeric matrix with samples as columns, and features as rows.
#' In the deconvolution formula 'test.data' is denoated as Y.
#' @param estimate.c.type string, either "non_negative", or "direct". Indicates how the algorithm finds the solution of
#' \eqn{arg min_C ||diag(g)(Y - XC)||_2}. If estimate.c.type is set to "direct" there is no regularization
#' (see \code{\link{estimate_c}}),
#' if estimate.c.type is set to "non_negative" the estimates "C"
#' must not be negative (non-negative least squares) (see (see \code{\link{estimate_nn_c}}))
#' @param feature.subset numeric or a vector of strings. If it is a numeric, "subset" features will be picked from the
#' 'explained correlation' ranking (if 'feature.subset' < 1, this is the fraction of feature, if 'feature.subset' > 1
#' it is the total amount. If it is a vector of strings, these genes will be used.
#' @param main string, used as title in the plot
#'
#' @return ggplot object
#' @export
#'
#' @import ggplot2
#'
ggplot_heatmap <- function(DTD.model,
                           X.matrix = NA,
                           test.data=NULL,
                           estimate.c.type,
                           main = "",
                           feature.subset = 100){
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

  # safety check: test.data is moved into "expl.cor" part,
  # because if subset is a list of genes, I don't need it
  # estimate.c.type is moved as well

  # safety check: subset
  if(length(feature.subset) != 1){
    if(is.character(feature.subset)){
      useable.subset <- intersect(feature.subset, rownames(X.matrix))
      if(length(useable.subset) == 0){
        stop("In ggplot_heatmap: 'feature.subset' is provided as vector of character. However, none of them can be found in rownames(X.matrix).")
      }else{
        features <- useable.subset
      }
    }
  }else{
    test <- test_numeric(feature.subset,
                         output.info = c("ggplot_heatmap", "feature.subset"),
                         min = 0,
                         max = Inf)
    if(feature.subset <= 1){ # fraction is provided
      feature.subset <- round(nrow(X.matrix) * feature.subset)
    }
    feature.subset <- min(feature.subset, nrow(X.matrix))
  }
  # end -> subset
  # safety check: main
  useable.main <- try(as.character(main), silent = TRUE)
  if(any(grepl(x = useable.main, pattern = "Error"))){
    stop("In ggplot_cv: provided 'main' can not be used as.character.")
  }
  # end -> main

  if(is.null(test.data) && !exists("features")){
    features <- rownames(X.matrix)
  }

  if(!exists("features")){
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

  tmp.X.matrix <- X.times.g[features, ]

  cell.type.cluster <- hclust(dist(x = t(tmp.X.matrix),
                                   method = "euclidean"),
                              method = "average")
  feature.cluster <- hclust(dist(x = tmp.X.matrix,
                                   method = "euclidean"),
                              method = "average")

  melted.X <- reshape2::melt(tmp.X.matrix,
                             value.name = "expression")
  melted.X$Var1 <- factor(
    x = melted.X$Var1,
    levels = rownames(tmp.X.matrix)[feature.cluster$order]
  )

  melted.X$Var2 <- factor(
    x = melted.X$Var2,
    levels = colnames(tmp.X.matrix)[cell.type.cluster$order]
  )



  pic0 <- ggplot(melted.X, aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = log(expression))) +
    scale_fill_gradient(
      low = "darkblue",
      high = "yellow"
      ) +
    xlab("features") +
    ylab("Cell types") +
    ggtitle(main) +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1))
  return(pic0)
}