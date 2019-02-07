#' Visualizing cross validation output
#'
#' Plot the regularization parameter lambda vs the resulting loss function value.
#' For an example see " Cross validation" in the package vignette `browseVignettes("DTD")`
#'
#' @param DTD.model : list as returned by \code{\link{train_correlatio_model}} or \code{\link{DTD_cv_lambda}}
#' @param main string, used as title in the plot. Defaults to ""
#' @param LAMBDA.TRANS.FUN function, will be applied to the lambda sequence,
#' to get equdistance x ticks (defaults to log2)
#' @param x.lab string, or expression, will be used as the label on the x axis, defaults to "log2(lambda)"
#' @param y.lab string, or expression, will be used as the label on the y axis, defaults to "Loss function"
#'
#' @import ggplot2
#'
#' @return a ggplot plot object
#' @export
ggplot_cv <- function(DTD.model,
                      main = "",
                      LAMBDA.TRANS.FUN = log2,
                      x.lab = expression("log2(" * lambda * ")"),
                      y.lab = "Loss function") {

  # test if DTD.model can be used for plotting:
  if (is.list(DTD.model) && ("cv.obj" %in% names(DTD.model))) {
    cross.val.frame <- DTD.model$cv.obj
  } else {
    if (is.data.frame(DTD.model) && all(c("nonZero", "nFoundModels") %in% rownames(DTD.model))) {
      cross.val.frame <- DTD.model
    }
  }
  if (!exists("cross.val.frame")) {
    stop("ggplot_cv: DTD.model does not fit")
  }

  # safety checks:
  if (!is.data.frame(cross.val.frame) ||
    !any(grepl("nonZero", rownames(cross.val.frame))) ||
    !any(grepl("nFoundModels", rownames(cross.val.frame))) ||
    !is.numeric(LAMBDA.TRANS.FUN(2))
  ) {
    stop("In ggplot_cv: passed \"cross.val.frame\" frame does not fit or LAMBDA.TRANS.FUN does not return numeric.")
  }

  # Extract number of zero coefficient per lambda,
  tmp.nZero <- cross.val.frame["nonZero", ]
  # and number of found models
  tmp.nModel <- cross.val.frame["nFoundModels", ]

  # transpose the data frame, and remove the "nonZero" and "nFoundModels" entries
  tmp.crossValFrame <- as.data.frame(t(cross.val.frame[-which(rownames(cross.val.frame) %in% c("nonZero", "nFoundModels")), ]))

  # add the used lambdas as a numeric column
  tmp.crossValFrame$lambda <- as.numeric(rownames(tmp.crossValFrame))

  # melt the frame to long format
  cvf.melt <- melt(tmp.crossValFrame, value.name = "Loss", id.vars = "lambda")
  # apply lambda transformation function
  cvf.melt$lambda <- LAMBDA.TRANS.FUN(cvf.melt$lambda)

  # check if LAMBDA.TRANS.FUN has been changed, but x.lab not. (making the plot uninterpretable)
  if (LAMBDA.TRANS.FUN(2) != 1 && grepl(x.lab, pattern = "^log2(.)$")) {
    x.lab <- expression("transformed " * lambda)
  }

  # plot the figure ...
  pic <- ggplot2::ggplot(cvf.melt, aes_string(x = "lambda", y = "Loss", group = "lambda")) +
    ggplot2::geom_boxplot(na.rm = TRUE) +
    ggplot2::xlab(x.lab) +
    ggplot2::ylab(y.lab) +
    ggplot2::scale_x_continuous(sec.axis = sec_axis(
      ~.,
      breaks = LAMBDA.TRANS.FUN(tmp.crossValFrame$lambda),
      labels = tmp.nZero,
      name = "# of non-zero-coefficients"
    )) +
    ggplot2::ggtitle(main)

  # ... and return it
  return(pic)
}
