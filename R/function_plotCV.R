#' Visualizing cross validation output
#'
#' @param crossValFrame
#' @param main
#' @param LAMBDA.TRANS.FUN
#'
#' @return
#' @export
#'
#' @examples
ggplot_cv <- function(crossValFrame, main = "", LAMBDA.TRANS.FUN = log2){
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

  pic <- ggplot2::ggplot(crossValFrame, aes(x = lambda.trans, y = cvm)) +
         ggplot2::geom_errorbar(aes_string(ymin = "cvlo", ymax="cvup")) +
         ggplot2::geom_point(col = "red") +
         ggplot2::xlab(expression("log2("*lambda*")")) +
         ggplot2::ylab("Loss function") +
         ggplot2::scale_x_continuous(sec.axis = sec_axis(
                                ~.,
                              breaks = crossValFrame$lambda.trans,
                              labels =  paste0(crossValFrame$nfoundModels,  "\n", crossValFrame$nzero),
                              name = "# of Found Models \n # of zero-coefficients")) +
         ggplot2::ggtitle(main)
}
