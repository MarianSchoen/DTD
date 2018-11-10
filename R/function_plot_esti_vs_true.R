#' Plot true vs estimated cell composition
#'
#' 'ggplot_true_vs_esti' can be used to evaluate a trained DTD model,
#' by plotting known/true quantities versus estimated quantities per cell type. Its main inputs are two quantity matrices.
#' For an example see
#' \code{\link{descent_generalized_fista}}
#'
#' @param estimatedC matrix, with cell types as rows, and mixtures as columns.
#' @param trueC matrix, with cell types as rows, and mixtures as columns.
#' @param norm.columnwise boolean, should every column be normalized to sum of 1?
#' @param title string, title for plot
#' @param color.indi vector with length of ncol(trueC), used for colouring in ggplot
#'
#' @return ggplot object
#' @export
#'
ggplot_true_vs_esti <- function(estimatedC, trueC, norm.columnwise = TRUE, title = "", color.indi = NA){
  if(!all(dim(estimatedC) == dim(trueC))){
    stop("Dimension of C are uniequal ...\n")
  }
  if(any(is.na(color.indi))){
    color.indi <- rep(1, ncol(estimatedC))
  }

  if(norm.columnwise){
    estimatedC <- apply(estimatedC, 2, function(x){return(x/sum(x))})
  }

  #check if both rows are sorted in the same way:
  if(any(rownames(estimatedC) != rownames(trueC))){
    # if not, resort them
    if(all(rownames(estimatedC) %in% rownames(trueC))){
      trueC <- trueC[rownames(estimatedC), ]
    }else{
      stop("ggplot_true_vs_esti: There are different rownames in estimatedC and trueC")
    }
  }

  cor.list <- c()
  for(l1 in 1:nrow(estimatedC)){
    cor.list <- c(cor.list, stats::cor(estimatedC[l1, ], trueC[l1, ]))
  }
  names(cor.list) <- rownames(estimatedC)

  tit <- paste0("Overall Correlation: ", format(mean(cor.list), digits=2), "; ", title)


  esti.frame <- as.data.frame(t(estimatedC))
  esti.frame$colorIndi <- factor(color.indi)
  esti.frame$names <- rownames(esti.frame)
  estimated <- melt(esti.frame, id.vars = c("colorIndi", "names"))
  # resort the frame, that it can be matched with true:
  estimated <- estimated[order(estimated$names, as.character(estimated$variable)), ]

  true.frame <- as.data.frame(t(trueC))
  true.frame$colorIndi <- factor(color.indi)
  true.frame$names <- rownames(true.frame)
  true <- melt(true.frame, id.vars = c("colorIndi", "names"))
  # resort the frame, that it can be matched with estimated:
  true <- true[order(true$names, as.character(true$variable)), ]

  if(all(true$names == estimated$names) && all(true$variable == estimated$variable)){
    complete <- estimated
    complete$true <- true$value
  }

  oldLabels <- levels(complete$variable)
  newLabels <- c()
  for(l1 in oldLabels){
    newLabels <- c(newLabels, paste0(l1, "\nCor: ", format(cor.list[l1], digits=2)))
  }

  levels(complete$variable) <- newLabels

  pic <- ggplot(complete, aes_string(y="value", x="true", color = "variable")) +
    ylab("estimated") + xlab("true") +
    ggtitle(tit) +
    facet_grid(.~variable, scales = "free") +
    theme(axis.text.x = element_text(angle=90))

  if(length(table(color.indi)) == 1){
    pic <- pic + geom_point()
  }else{
    pic <- pic + geom_point(aes_string(shape = "colorIndi"))
  }
  return(pic)
}
