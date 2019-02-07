#' Plot true vs estimated cell composition
#'
#' 'ggplot_true_vs_esti' can be used to evaluate a trained DTD model, in the way that
#' it plots known/true quantities versus estimated quantities per cell type. Its main inputs are two quantity matrices.
#' For an example see section "Correlation per cell type" in the package vignette `browseVignettes("DTD")`
#'
#' @param estimated.c matrix, with cell types as rows, and mixtures as columns
#' @param true.c matrix, with cell types as rows, and mixtures as columns
#' @param norm.mixturewise boolean, should every sample (=> column) be normalized to sum of 1? (Defaults to FALSE)
#' @param norm.typewise boolean, should every type (=> row) be normalized to range from 0 to 1? (Defaults to FALSE)
#' @param title string, title for plot (Defaults to "")
#' @param shape.indi vector with length of ncol(true.c), used as parameter shape in ggplot.
#' Idea is to mark samples from different origin. If no additional(Defaults to NA)
#' @param show.legend boolean, should an additional legend be plotted? Notice, in this function a figure per type will be generated.
#' In the title of each subfigure the correlation of the included type and the type is shown.
#' This parameter only sets the additional overall legend(Defaults to FALSE)
#'
#' @import ggplot2
#' @import reshape2
#'
#' @return ggplot object
#' @export
#'
ggplot_true_vs_esti <- function(estimated.c,
                                true.c,
                                norm.mixturewise = FALSE,
                                norm.typewise = FALSE,
                                title = "",
                                shape.indi = NA,
                                show.legend = FALSE) {
  # Safety checks:
  if (!all(dim(estimated.c) == dim(true.c))) {
    stop("Dimension of C are uniequal ...\n")
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
      stop("ggplot_true_vs_esti: There are different rownames in estimated.c and true.c")
    }
  }

  # In the title of the subplots we add the correlation per type, therefore:
  cor.list <- c()
  for (l1 in 1:nrow(estimated.c)) {
    cor.list <- c(cor.list, stats::cor(estimated.c[l1, ], true.c[l1, ]))
  }
  names(cor.list) <- rownames(estimated.c)

  # Adjust title:
  tit <- paste0("Overall Correlation: ", format(mean(cor.list), digits = 2))
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
