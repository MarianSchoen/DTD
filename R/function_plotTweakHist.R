#' Histogram of the g/tweak vector
#'
#' The resulting distribution of g is visualized best in a histogram.
#' As the g vector may get spread over a big range,
#' you can provide a transformation function.
#' This will be applied on the g/tweak vector. If the provided transformation
#' function returns NA, the function will state the number of NAs in the title
#' of the plot.
#' For an example see section "Histogram of g-vector" in the package vignette
#' `browseVignettes("DTD")`
#'
#' @param DTD.model either a numeric vector with length of nrow(X), or a list
#' returned by \code{\link{train_deconvolution_model}},
#' \code{\link{DTD_cv_lambda_cxx}}, or \code{\link{descent_generalized_fista}}.
#' In the equation above the DTD.model provides the vector g.
#' @param n.bins : positive integer, number of bins in the histogram
#' @param G.TRANSFORM.FUN function, that expects a vector of numerics,
#' and returns a vector of the same length. Will be applied on each intermediate
#''g' vector. Set 'G.TRANSFORM.FUN' to identity if no transformation is required.
#' If you change 'G.TRANSFORM.FUN' don't forget to adjust the 'x.lab' parameter.
#' @param title string, additionally title
#' @param x.lab string, used as x label on the plot
#'
#' @import ggplot2
#'
#' @return ggplot object
#' @export
ggplot_ghistogram <- function(DTD.model,
                              n.bins = NA,
                              G.TRANSFORM.FUN = log10p1,
                              title = "",
                              x.lab = "log10(g+1)") {

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

  # safety check: G.TRANSFORM.FUN
  if(!is.function(G.TRANSFORM.FUN)){
    stop("In ggplot_ghistogram: 'G.TRANSFORM.FUN' is not a function")
  }
  # end -> G.TRANSFORM.FUN

  # Transformation:
  g_vec <- suppressWarnings(G.TRANSFORM.FUN(tweak))

  # Get number of na:
  nums.NA <- sum(is.na(g_vec)) + sum(is.infinite(g_vec))
  # If there are na, adjust "title"
  if (nums.NA > 0) {
    title <- paste0("Due to G.TRANSFORM.FUN, there are ", nums.NA, " missing values \n", title)
  }
  if(nums.NA == length(tweak)){
    stop("In ggplot_ghistogram: after 'G.TRANSFORM.FUN(tweak), all values are NA or infinite.")
  }

  # draw the picture:
  pic <- ggplot2::ggplot(data = NULL, aes(x = g_vec)) +
    ggplot2::geom_histogram(bins = n.bins) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(x.lab)

  return(pic)
}
