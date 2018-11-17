#' Visualizing cross validation output
#'
#' Plot the regularization parameter lambda vs the resulting loss function.
#' For an example see `browseVignettes("DTD")`
#'
#' @param crossValFrame data frame, as return by the DTD_cv_lambda function
#' @param main string, used as title in the plot
#' @param LAMBDA.TRANS.FUN function, will be applied to the lambda sequence, to get equdistance x ticks (defaults to log2)
#' @param x.lab string, or expression, will be used as the label on the x axis
#' @param y.lab string, or expression, will be used as the label on the y axis
#'
#' @import ggplot2
#'
#' @return a ggplot plot object
#' @export
ggplot_cv <- function(crossValFrame, main = "",
                      LAMBDA.TRANS.FUN = log2,
                      x.lab = expression("log2("*lambda*")"),
                      y.lab = "Loss function"){
  # safety check, if the passed frame matches:
  if(!is.data.frame(crossValFrame) ||
     is.null(crossValFrame$lambda) ||
     is.null(crossValFrame$cvm)    ||
     is.null(crossValFrame$cvlo)   ||
     is.null(crossValFrame$cvup)   ||
     is.null(crossValFrame$nfoundModels) ||
     is.null(crossValFrame$nzero)
    ) {
    stop("In plot.cv: passed \"crossValFrame\" frame does not fit")
  }
  crossValFrame$lambda.trans <- LAMBDA.TRANS.FUN(crossValFrame$lambda)

  # check if LAMBDA.TRANS.FUN has been changed, but x.lab not. (making the plot uninterpretable)
  if(LAMBDA.TRANS.FUN(2) != 1 && grepl(x.lab, pattern = "^log2(.)$")){
    x.lab <- expression("transformed "*lambda)
  }


  pic <- ggplot2::ggplot(crossValFrame, aes(x = lambda.trans, y = cvm)) +
         ggplot2::geom_errorbar(aes_string(ymin = "cvlo", ymax="cvup")) +
         ggplot2::geom_point(col = "red") +
         ggplot2::xlab(x.lab) +
         ggplot2::ylab(y.lab) +
         ggplot2::scale_x_continuous(sec.axis = sec_axis(
                                ~.,
                              breaks = crossValFrame$lambda.trans,
                              labels =  paste0(crossValFrame$nfoundModels,  "\n", crossValFrame$nzero),
                              name = "# of Found Models \n # of zero-coefficients")
                              ) +
         ggplot2::ggtitle(main)

  return(pic)
}
