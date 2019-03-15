#' Histogram of the g/tweak vector
#'
#' The resulting distribution of g are visualized best in a histogram.
#' As the g vector may get spread over a big range, you can provide a transformation function.
#' This will be applied on the g/tweak vector. If the provided transformation function returns NA,
#' the function will state the number of NAs in the title of the plot.
#' For an example see section "Histogram of g-vector" in the package vignette `browseVignettes("DTD")`
#'
#' @param DTD.model : list, as returned by \code{\link{train_correlatio_model}}, \code{\link{DTD_cv_lambda}},
#' or \code{\link{descent_generalized_fista}}
#' @param n.bins : integer, number of bins in the histogram. Defaults to 50
#' @param TRANSFORM.FUN : function that expects a list of floats, and returns a list of floats
#' Defaults to identity.
#' @param main string, used as title of the plot. Defaults to ""
#' @param x.lab string, used as x label of the plot. Defaults to ""
#'
#' @import ggplot2
#'
#' @return ggplot object
#' @export
#'
ggplot_ghistogram <- function(DTD.model,
                              n.bins = 50,
                              TRANSFORM.FUN = DTD::identity,
                              main = "",
                              x.lab = "g-vec") {

  # test if DTD.model can be used for plotting:
  if (is.list(DTD.model) && ("Tweak" %in% names(DTD.model))) {
    tweak <- DTD.model$Tweak
  }
  if (is.list(DTD.model) && ("best.model" %in% names(DTD.model))) {
    tweak <- DTD.model$best.model$Tweak
  }
  if (is.numeric(DTD.model)) {
      tweak <- DTD.model
    }
  if (!exists("tweak")) {
    stop("In ggplot_ghistogram: 'DTD.model' does not fit")
  }else{
    test <- test_tweak_vec(tweak.vec = tweak,
                           output.info = c("ggplot_ghistogram", "tweak"))
  }

  # safety check: n.bins
  test <- test_integer(test.value = n.bins,
                       output.info = c("ggplot_ghistogram", "n.bins"),
                       min = 1,
                       max = length(tweak))
  # end -> n.bins

  # safety check: main
  useable.main <- try(as.character(main), silent = TRUE)
  if(any(grepl(x = useable.main, pattern = "Error"))){
    stop("In ggplot_ghistogram: provided 'main' can not be used as.character.")
  }
  # end -> main
  # safety check: x.lab
  useable.x.lab <- try(as.character(x.lab), silent = TRUE)
  if(any(grepl(x = useable.x.lab, pattern = "Error"))){
    stop("In ggplot_ghistogram: provided 'x.lab' can not be used as.character.")
  }
  # end -> x.lab


  # Transformation:
  g_vec <- suppressWarnings(TRANSFORM.FUN(tweak))

  # Get number of na:
  nums.NA <- sum(is.na(g_vec)) + sum(is.infinite(g_vec))
  # If there are na, adjust "main"
  if (nums.NA > 0) {
    main <- paste0("Due to TRANSFORM.FUN, there are ", nums.NA, " missing values \n", main)
  }

  # draw the picture:
  pic <- ggplot2::ggplot(data = NULL, aes(x = g_vec)) +
    ggplot2::geom_histogram(bins = n.bins) +
    ggplot2::ggtitle(main) +
    ggplot2::xlab(x.lab)


  return(pic)
}
