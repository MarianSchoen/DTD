#' Plot true vs estimated cell composition
#'
#' 'ggplot_true_vs_esti' can be used to evaluate a trained DTD model, in the way that
#' it plots known/true quantities versus estimated quantities per cell type.
#' As input, 'ggplot_true_vs_esti' expects a DTD.model, and test.data.
#' For an example see section "Correlation per cell type" in the package vignette `browseVignettes("DTD")`
#'
#' @param norm.mixturewise logical, in the known AND estimated quantity matrix, should every mixture (=> column) be normalized to sum of 1? (Defaults to FALSE)
#' @param norm.typewise logical, in the known AND estimated quantity matrix, should every type (=> row) be normalized to range from 0 to 1? (Defaults to FALSE)
#' @param title string, additionally title (default "")
#' @param shape.indi vector with length of ncol(test.data$quantities), will be passed to shape argument of geom_point.
#' Idea is to mark samples from different origin. Defaults to NA. If set to NA, the shape will not be changed.
#' @param show.legend logical, should an additional legend be plotted?
#' Notice, this function generates a plot, holding a subfigure for each type of the deconvolution.
#' In every subfigure, the cell type, and the corresponding correlation is shown.
#' In the title of each subfigure the correlation of the included type and the type is shown.
#' This parameter only controls the additional legend. (Defaults to FALSE)
#' @param DTD.model either a numeric vector with length of nrow(X),
#' or a list returned by \code{\link{train_deconvolution_model}}, \code{\link{DTD_cv_lambda}},
#' or\code{\link{descent_generalized_fista}}.
#' @param X.matrix numeric matrix with cells as columns, and features as rows.
#'  Reference matrix X of the DTD problem. X.matrix can be set to NA (default), if the DTD.model
#'  includes the reference matrix X (default for \code{\link{train_deconvolution_model}})
#' @param test.data list of two matrices, named "mixtures" and "quantities".
#' For examples see \code{\link{mix_samples}}, \code{\link{mix_samples_with_jitter}}
#' or the package vignette `browseVignettes("DTD")`.
#' @param estimate.c.type string, either "non_negative", or "direct". Indicates how the algorithm finds the solution of
#' \eqn{arg min_C ||diag(g)(Y - XC)||_2}. \cr
#' If estimate.c.type is set to "direct" there is no regularization (see \code{\link{estimate_c}}),\cr
#' if estimate.c.type is set to "non_negative" the estimates "C"
#' must not be negative (non-negative least squares) (see (see \code{\link{estimate_nn_c}}))
#'
#' @import ggplot2
#' @import reshape2
#'
#' @return ggplot object
#' @export
#'
ggplot_true_vs_esti <- function(DTD.model,
                                X.matrix = NA,
                                test.data,
                                norm.mixturewise = FALSE,
                                norm.typewise = FALSE,
                                estimate.c.type,
                                title = "",
                                shape.indi = NA,
                                show.legend = FALSE) {

  if (is.list(DTD.model)) {
    if ("best.model" %in% names(DTD.model)) {
      tweak <- DTD.model$best.model$Tweak
    } else {
      if (!"Tweak" %in% names(DTD.model)) {
        stop("In ggplot_true_vs_esti: There is no Tweak entry in the 'DTD.model'")
      } else {
        tweak <- DTD.model$Tweak
      }
    }
  } else {
    if(is.numeric(DTD.model)){
      tweak <- DTD.model
    }else{
      stop("In ggplot_true_vs_esti: DTD.model is neither a list nor a numeric vector")
    }
  }


  if(is.list(test.data) && length(test.data) == 2){
    if(!all(c("quantities", "mixtures") %in%  names(test.data))){
      stop("In ggplot_true_vs_esti: entries of 'test.data' must be named 'quantities' and 'mixtures'.")
    }else{
      if(!is.matrix(test.data$mixtures)){
        stop("In ggplot_true_vs_esti: 'test.data$mixtures' is not a matrix.")
      }
      if(!is.matrix(test.data$quantities)){
        stop("In ggplot_true_vs_esti: 'test.data$quantities' is not a matrix. Therefore, 'test.data' can not be used.")
      }
    }
  }else{
    stop("In ggplot_true_vs_esti: 'test.data' must be provided as a list with two entries: 'quantities' and 'mixtures'.")
  }
  if(!is.matrix(X.matrix) || any(is.na(X.matrix)) || !is.numeric(X.matrix)){
    if (is.list(DTD.model) && "reference.X" %in% names(DTD.model)) {
      message("In ggplot_true_vs_esti: provided 'X.matrix' could not be used, therefore used: 'DTD.model$reference.X'")
      X.matrix <- DTD.model$reference.X
    }else{
      stop("In ggplot_true_vs_esti: can't use 'X.matrix'. Neither via X.matrix argument nor included in the DTD.model.")
    }
  }

  # safety check: title
  useable.ylab <- try(as.character(title), silent = TRUE)
  if(any(grepl(x = useable.ylab, pattern = "Error"))){
    stop("In ggplot_true_vs_esti: provided 'title' can not be used as.character.")
  }
  # end -> title

  # safety check: norm.typewise
  test <- test_logical(test.value = norm.typewise,
                       output.info = c("ggplot_true_vs_esti", "norm.typewise"))
  # end -> norm.typewise
  # safety check: show.legend
  test <- test_logical(test.value = show.legend,
                       output.info = c("ggplot_true_vs_esti", "show.legend"))
  # end -> show.legend

  # safety check: norm.mixturewise
  test <- test_logical(test.value = norm.mixturewise,
                       output.info = c("ggplot_true_vs_esti", "norm.mixturewise"))
  # end -> norm.mixturewise

  # safety check: estimate.c.tye
  ESTIMATE_C_FUN <- test_c_type(test.value = estimate.c.type,
                                output.info = c("ggplot_true_vs_esti", "estimate.c.type"))
  # end -> estimate.c.type
  # safety check: title
  title <- test_string(test.value = title,
                       output.info = c("ggplot_true_vs_esti", "title"))
  # end -> title

  estimated.c <- ESTIMATE_C_FUN(X.matrix = X.matrix,
                                new.data = test.data$mixtures,
                                DTD.model = tweak)

  true.c <- test.data$quantities

  if(!all(dim(estimated.c) == dim(true.c))){
    stop("In ggplot_true_vs_esti: dimension of estimated C, and C in 'test.data$quantities' differ" )
  }


  # If shape.indi has not been set, initialize it
  if (any(is.na(shape.indi))) {
    shape.indi <- rep(1, ncol(estimated.c))
  }

  # Norm every mixture (column of the data) to sum of 1
  if (norm.mixturewise) {
    estimated.c <- apply(estimated.c, 2, function(x) {
      return(x / sum(x))
    })
    true.c <- apply(true.c, 2, function(x) {
      return(x / sum(x))
    })
  }
  # Norm every type to range from 0 to 1
  if (norm.typewise) {
    estimated.c <- t(apply(estimated.c, 1, function(x) {
      return((x - min(x)) / (max(x) - min(x)))
    }))
    true.c <- t(apply(true.c, 1, function(x) {
      return((x - min(x)) / (max(x) - min(x)))
    }))
  }

  # check if both rows are sorted in the same way:
  if (any(rownames(estimated.c) != rownames(true.c))) {
    # if not, resort them
    if (all(rownames(estimated.c) %in% rownames(true.c))) {
      true.c <- true.c[rownames(estimated.c), ]
    } else {
      stop("In ggplot_true_vs_esti: There are different rownames in estimated.c and true.c")
    }
  }

  # In the title of the subplots we add the correlation per type, therefore:
  cor.list <- c()
  for (l1 in 1:nrow(estimated.c)) {
    if(stats::sd(true.c[l1, ]) != 0){ # => can't calculate corelation
      cor.list <- c(cor.list, stats::cor(estimated.c[l1, ], true.c[l1, ]))
    }else{
      cor.list <- c(cor.list, 0)
    }
  }
  names(cor.list) <- rownames(estimated.c)

  # Adjust title:
  tit <- paste0("Overall Correlation: ", format(mean(cor.list, na.rm = TRUE), digits = 2))
  if (title != "") {
    tit <- paste0(tit, ";", title)
  }


  # A numeric matrix is hard to plot with ggplot. Therefore, we convert it to long format.
  # Therefore, transpose it:
  esti.frame <- as.data.frame(t(estimated.c))
  # Add the information per mixture (as it is transposed)
  esti.frame$shapeIndi <- factor(shape.indi)
  esti.frame$names <- rownames(esti.frame)
  # And melt it
  estimated <- reshape2::melt(esti.frame, id.vars = c("shapeIndi", "names"))
  # resort the frame, that it can be matched with true:
  estimated <- estimated[order(estimated$names, as.character(estimated$variable)), ]

  # Same for the true.c
  true.frame <- as.data.frame(t(true.c))
  true.frame$shapeIndi <- factor(shape.indi)
  true.frame$names <- rownames(true.frame)
  true <- melt(true.frame, id.vars = c("shapeIndi", "names"))
  # resort the frame, that it can be matched with estimated:
  true <- true[order(true$names, as.character(true$variable)), ]

  # We resorted both frames, and we checked beforehand if all samples match.
  # Therefore the frames should match:
  if (all(true$names == estimated$names) && all(true$variable == estimated$variable)) {
    complete <- estimated
    complete$true <- true$value
  } else {
    stop("Matching Problem between melted true.c and melted estimated.c.")
  }

  # Adjust the levels of the factor "variable", in order to include the correlation
  oldLabels <- levels(complete$variable)
  newLabels <- c()
  for (l1 in oldLabels) {
    newLabels <- c(newLabels, paste0(l1, "\nCor: ", format(cor.list[l1], digits = 2)))
  }

  levels(complete$variable) <- newLabels

  # and plot it
  pic <- ggplot2::ggplot(complete, aes_string(y = "value", x = "true", color = "variable")) +
    ggplot2::ylab("estimated") + xlab("true") +
    ggplot2::ggtitle(tit) +
    ggplot2::facet_grid(. ~ variable, scales = "free") +
    ggplot2::theme(axis.text.x = element_text(angle = 90))

  # If required, remove legend
  if (!show.legend) {
    pic <- pic + theme(legend.position = "none")
  }

  # And add points to the plot, if required with special shape.
  if (length(table(shape.indi)) == 1) {
    pic <- pic + ggplot2::geom_point()
  } else {
    pic <- pic + ggplot2::geom_point(aes_string(shape = "shapeIndi"))
  }
  return(pic)
}
