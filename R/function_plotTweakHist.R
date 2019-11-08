#' Histogram of the g/tweak vector
#'
#' The resulting distribution of g are visualized best in a histogram.
#' As the g vector may get spread over a big range, you can provide a transformation function.
#' This will be applied on the g/tweak vector. If the provided transformation function returns NA,
#' the function will state the number of NAs in the title of the plot.
#' For an example see section "Histogram of g-vector" in the package vignette `browseVignettes("DTD")`
#'
#' @param DTD.model either a numeric vector ,
#' or a list returned by \code{\link{train_deconvolution_model}}, \code{\link{DTD_cv_lambda}},
#' or\code{\link{descent_generalized_fista}}.
#' @param n.bins : positive integer, number of bins in the histogram. Defaults to 50
#' @param TRANSFORM.FUN : function that expects a list of floats, and returns a list of floats
#' Defaults to identity. (defaults to identity)
#' @param title string, additionally title (default "")
#' @param x.lab string, used as x label of the plot. Defaults to ""
#'
#' @import ggplot2
#'
#' @return ggplot object
#' @export
#'
ggplot_ghistogram <- function(DTD.model,
                              n.bins = NA,
                              TRANSFORM.FUN = log10,
                              title = "",
                              x.lab = "log10(g)") {

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
  if(any(is.na(n.bins))){
    n.bins <- round(0.25 * length(tweak))
  }
  test <- test_integer(test.value = n.bins,
                       output.info = c("ggplot_ghistogram", "n.bins"),
                       min = 1,
                       max = length(tweak))
  # end -> n.bins

  # safety check: title
  title <- test_string(
    test.value = title,
    output.info = c("ggplot_ghistogram", "titel")
  )
  # end -> title

  # safety check: x.lab
  x.lab <- test_string(
    test.value = x.lab,
    output.info = c("ggplot_ghistogram", "titel")
  )
  # end -> x.lab

  # safety check: TRANSFORM.FUN
  if(!is.function(TRANSFORM.FUN)){
    stop("In ggplot_ghistogram: 'TRANSFORM.FUN' is not a function")
  }
  # end -> TRANSFORM.FUN

  # Transformation:
  g_vec <- suppressWarnings(TRANSFORM.FUN(tweak))

  # Get number of na:
  nums.NA <- sum(is.na(g_vec)) + sum(is.infinite(g_vec))
  # If there are na, adjust "title"
  if (nums.NA > 0) {
    title <- paste0("Due to TRANSFORM.FUN, there are ", nums.NA, " missing values \n", title)
  }
  if(nums.NA == length(tweak)){
    stop("In ggplot_ghistogram: after 'TRANSFORM.FUN(tweak), all values are NA or infinite.")
  }

  # draw the picture:
  pic <- ggplot2::ggplot(data = NULL, aes(x = g_vec)) +
    ggplot2::geom_histogram(bins = n.bins) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(x.lab)

  return(pic)
}
