#' Plot the explained correlation per feature
#'
#' In order to assess the importance of a feature in the deconvolution process, we can exclude it from a trained model,
#' and observe the change of correlaion on a test set. If the correlation e.g. decreases by 1 %,
#' the gene explains 1 % correlation within the deconvolution model.
#' The 'ggplot_explained_correlation' function iteratively excludes each feature from the trained model,
#' and visualizes the explained correlation.
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
#' @param print.labels logical, should the feature names be plotted into the picture? (Defaults to FALSE)
#'
#' @return ggplot object
#' @export
#'
#' @import ggplot2
#'
ggplot_explained_correlation <- function(DTD.model,
                                         X.matrix = NA,
                                         test.data,
                                         estimate.c.type,
                                         print.labels = FALSE){

  # safety check: DTD.model
  if (is.list(DTD.model)) {
    if ("best.model" %in% names(DTD.model)) {
      tweak <- DTD.model$best.model$Tweak
    } else {
      if ("Tweak" %in% names(DTD.model)) {
        stop("In ggplot_explained_correlation: There is no Tweak entry in the 'DTD.model'")
      } else {
        tweak <- DTD.model$Tweak
      }
    }
  } else {
    if(is.numeric(DTD.model)){
      tweak <- DTD.model
    }else{
      stop("In ggplot_explained_correlation: DTD.model is neither a list nor a numeric vector")
    }
  }
  # end -> DTD.model
  # safety check: test.data
  if(is.list(test.data) && length(test.data) == 2){
    if(!all(c("quantities", "mixtures") %in%  names(test.data))){
      stop("In ggplot_explained_correlation: entries of 'test.data' must be named 'quantities' and 'mixtures'.")
    }else{
      if(!is.matrix(test.data$mixtures)){
        stop("In ggplot_explained_correlation: 'test.data$mixtures' is not a matrix.")
      }
      if(!is.matrix(test.data$quantities)){
        stop("In ggplot_explained_correlation: 'test.data$quantities' is not a matrix. Therefore, 'test.data' can not be used.")
      }
    }
  }else{
    stop("In ggplot_explained_correlation: 'test.data' must be provided as a list with two entries: 'quantities' and 'mixtures'.")
  }
  # end -> test.data
  # safety check: X.matrix
  if(!is.matrix(X.matrix) || any(is.na(X.matrix)) || !is.numeric(X.matrix)){
    if (is.list(DTD.model) && "reference.X" %in% names(DTD.model)) {
      message("In ggplot_explained_correlation: provided 'X.matrix' could not be used, therefore used: 'DTD.model$reference.X'")
      X.matrix <- DTD.model$reference.X
    }else{
      stop("In ggplot_explained_correlation: can't use 'X.matrix'. Neither via X.matrix argument nor included in the DTD.model.")
    }
  }
  # end -> X.matrix
  # safety check: estimate.c.tye
  test <- test_c_type(test.value = estimate.c.type,
                                output.info = c("ggplot_explained_correlation", "estimate.c.type"))
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

  maximum.cor <- -DTD.evCor.wrapper(tweak = tweak)

  manipulated.cor <- rep(NA, length(tweak))
  names(manipulated.cor) <- names(tweak)
  for(l.feature in names(tweak)){
    l.tweak <- tweak
    l.tweak[l.feature] <- 0
    manipulated.cor[l.feature] <- (maximum.cor + DTD.evCor.wrapper(tweak = l.tweak))*100
  }
  var.times.g <- rowSds(X.matrix) * tweak

  order.expl.cor <- order(manipulated.cor)
  plot.frame <- data.frame("tweak" = tweak[order.expl.cor],
                           "explained.correlation" = manipulated.cor[order.expl.cor])
  plot.frame$name <- factor(names(tweak)[order.expl.cor], levels = names(tweak)[order.expl.cor])


  pic0 <- ggplot(plot.frame, aes(x = name,
                                 y = explained.correlation,
                                 col = tweak)) +

    ylab("Percentage of explained correlation") +
    xlab("") +
    theme(legend.position = "none")


  if(print.labels){
    pic0 <- pic0 +
      geom_label(
        aes(label = name),
        size = 3) +
      theme(axis.text.x = element_blank())
  }else{
    pic0 <- pic0 + geom_point()
  }
  return(pic0)
}
