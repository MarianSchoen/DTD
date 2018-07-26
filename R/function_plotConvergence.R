#' Plot loss curve
#'
#' The "ggplot_correlation" function uses ggplot2 and reshape to visualize the decrease of the loss-function.
#'
#' As input parameter it needs the output of \code{\link{descent_generalized_fista}}. If fista has been evoked
#' with save_all_tweaks = T, and a test set is available, the loss-function can be evaluated for all tweak vectors,
#' and both training and test loss curves will be plotted (=> detect overfitting).\cr
#' The following example needs data generation, and evoking the FISTA algorithm as shown in
#' \code{\link{descent_generalized_fista}}.
#'
#' @param fista.output
#' @param test.set
#'
#' @return
#' @export
#'
#' @import ggplot2
#' @import reshape2
#'
#'
#' @examples
#' ggplot_correlation(fista.output = catch, test.set = NA)
ggplot_correlation <- function(fista.output, test.set){

  return("upsidupsi")
}
