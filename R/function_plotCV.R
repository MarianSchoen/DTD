#' Visualizing cross validation output
#'
#' Plot the regularization parameter lambda vs the resulting loss function value.
#' The lambda with the minimal mean loss in cross validation is marked with a red dot
#' For an example see " Cross validation" in the package vignette `browseVignettes("DTD")`
#'
#' @param DTD.model : list as returned by \code{\link{train_deconvolution_model}} or \code{\link{DTD_cv_lambda}}
#' @param title string, additionally title (default "")
#' @param LAMBDA.TRANS.FUN function, will be applied to the lambda sequence,
#' e.g. to get equidistance x ticks (defaults to log2)
#' @param x.lab string, or expression, will be used as the label on the x axis, defaults to "log2(lambda)"
#' @param y.lab string, or expression, will be used as the label on the y axis, defaults to "Loss function"
#' @param upper.x.axis.info string, either "non-zero" or "geometric-mean".
#' This information will be printed on the upper x axis.
#'
#' @import ggplot2
#'
#' @return a ggplot plot object
#' @export
ggplot_cv <- function(DTD.model,
                      title = "",
                      LAMBDA.TRANS.FUN = log2,
                      upper.x.axis.info = "non-zero",
                      x.lab = expression("log2(" * lambda * ")"),
                      y.lab = "Loss function") {
  # test if DTD.model can be used for plotting:
  if (is.list(DTD.model)) {
    if("cv.obj" %in% names(DTD.model)){
      cross.val.list <- DTD.model$cv.obj
    }else{
      if(length(DTD.model) != 0){
        length.per.lambda <- sapply(DTD.model, length)
        if(all(mean(length.per.lambda) == length.per.lambda)){
          cross.val.list <- DTD.model
        }
      }
    }
  } else {
    stop("In ggplot_cv: DTD.model must be provided as a list")
  }

  if (!exists("cross.val.list")) {
    stop("In ggplot_cv: 'DTD.model' does not fit")
  }

  # test upper.x.axis.info:
  if(!(upper.x.axis.info %in% c("non-zero", "geometric-mean"))){
    stop("In ggplot_cv: 'upper.x.axis.info' must either be 'non-zero' or 'geometric-mean'")
  }else{
    if(upper.x.axis.info == "non-zero"){
      upper.name <- "# of non-zero-coefficients"
      upper.info.FUN <- function(x){sum(x != 0)}
    }
    if(upper.x.axis.info == "geometric-mean"){
      upper.name <- "geometric mean of all coefficients"
      upper.info.FUN <- function(x){
        # definition of geometric mean:
        # in the DTD case, the "x" (=> tweak) must not be negativ, therefore I am not missing any:
        exp(sum(log(x[x > 0]), na.rm=TRUE) / length(x))
      }
    }
  }
  # end -> upper.x.axis.info

  # safety checks:
  if (!is.list(cross.val.list)){
    stop("In ggplot_cv: passed 'cross.val.list' does not fit. It must be a list")
  }

  useable.lambda.trans.fun <- try(LAMBDA.TRANS.FUN(2), silent = TRUE)
  if(any(grepl(x = useable.lambda.trans.fun, pattern = "Error"))){
    stop("In ggplot_cv: provided 'LAMBDA.TRANS.FUN' does not return a numeric, if called with '2'.")
  }

  title <- test_string(test.value = title, output.info = c("ggplot_cv", "title"))
  x.lab <- test_string(test.value = x.lab, output.info = c("ggplot_cv", "x.lab"))
  y.lab <- test_string(test.value = y.lab, output.info = c("ggplot_cv", "y.lab"))
  title <- test_string(test.value = title, output.info = c("ggplot_cv", "title"))

  # for the cross validation plot, we need the following things:
  #   - within a lambda, all test fold values
  #   - the lambdas included in the data frame (due to melting, and plotting)
  #   - the lambda itself, and the mean of the picked lambda
  #   - the "upper.x.axis.info"

  # after the cross validation, find the lambda with best evaluation score
  pick.test.results.function <- function(lambda.list){
    tmp <- lapply(lambda.list, function(each.fold){
      if("cor.test" %in% names(each.fold)){
        return(each.fold$cor.test)
      }else{
        return(Inf)
      }
    })
    test.vec <- unlist(tmp, use.names = FALSE)
    return(test.vec)
  }
  # pick the average mean per lambda:
  test.result.per.lambda <- lapply(cross.val.list,
                                   pick.test.results.function)

  # transform it to a frame => for plotting
  test.results.frame <- data.frame(
    matrix(
      unlist(test.result.per.lambda),
      nrow=length(test.result.per.lambda),
      byrow=TRUE)
    )
  rownames(test.results.frame) <- names(test.result.per.lambda)
  # for each model, pick "upper.x.axis" information (with same function/lapply construct):
  pick.upper.x.fun <- function(lambda.list){
    tmp <- lapply(lambda.list, function(each.fold){
      if("Tweak" %in% names(each.fold)){
        return(upper.info.FUN(each.fold$Tweak))
      }else{
        return(NA)
      }
    })
    upper.info.vec <- mean(
      unlist(
        lapply(
          tmp,
          mean,
          na.rm = TRUE),
        use.names = FALSE),
      na.rm = TRUE
      )
    return(upper.info.vec)
  }
  upper.info.per.lambda <- lapply(cross.val.list,
                                  pick.upper.x.fun)

  upper.info.per.lambda <- round(x = unlist(upper.info.per.lambda, use.names = FALSE), digits = 2)

  # if DTD.model is build by train_correlation_model => pick lambda from the best.model
  if(is.list(DTD.model) && "best.model" %in% names(DTD.model)){
    used.lambda <- LAMBDA.TRANS.FUN(DTD.model$best.model$Lambda)
    used.lambda.median <- stats::median(as.numeric(test.results.frame[as.character(DTD.model$best.model$Lambda), ]))
  }else{
    # find mean per lambda:
    median.per.lambda <- apply(test.results.frame, 1, stats::median, na.rm = TRUE)
    used.lambda.tmp <- as.numeric(names(which.min(median.per.lambda)))
    used.lambda <- LAMBDA.TRANS.FUN(used.lambda.tmp)
    used.lambda.median <- min(median.per.lambda, na.rm = TRUE)
  }
  # add the used lambdas as a numeric column
  test.results.frame$lambda <- as.numeric(rownames(test.results.frame))

  # melt the frame to long format
  cvf.melt <- melt(test.results.frame, value.name = "Loss", id.vars = "lambda")
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
      breaks = LAMBDA.TRANS.FUN(test.results.frame$lambda),
      labels = unlist(upper.info.per.lambda, use.names = FALSE),
      name = upper.name
    )) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_point(data = NULL,
                        aes(x = used.lambda,
                            y = used.lambda.median),
                        color = "red")

  # ... and return it
  return(pic)
}
