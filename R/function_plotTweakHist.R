#' Histogram of the g/tweak vector
#'
#' The resulting distribution of g are visualized best in a histogram.
#' As the g vector may get spread over a big range, you can provide a transformation function.
#' This will be applied on the g/tweak vector. If the provided transformation function returns NA,
#' the function will state the number of NAs in the title of the plot.
#' For an example see `browseVignettes("DTD")`
#'
#' @param fista.output : list with "Tweak" entry. The result of a \code{\link{descent_generalized_fista}} call
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
ggplot_ghistogram <- function(fista.output,
                             n.bins = 50,
                             TRANSFORM.FUN = DTD::identity,
                             main="",
                             x.lab = "g-vec"){

  # Transformation:
  g_vec <- suppressWarnings(TRANSFORM.FUN(fista.output$Tweak))

  # Get number of na:
  nums.NA <- sum(is.na(g_vec)) + sum(is.infinite(g_vec))
  # If there are na, adjust "main"
  if(nums.NA > 0){
    main <- paste0("Due to TRANSFORM.FUN, there are ", nums.NA, " missing values \n", main)
  }

  # draw the picture:
  pic <- ggplot2::ggplot(data = NULL, aes(x=g_vec)) +
            ggplot2::geom_histogram(bins=n.bins) +
            ggplot2::ggtitle(main) +
            ggplot2::xlab(x.lab)


  return(pic)
}
