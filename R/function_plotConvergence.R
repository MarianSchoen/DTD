#' Plot loss curve
#'
#' The "ggplot_convergence" function uses ggplot2 and reshape to visualize the decrease of the loss-function after
#' a model has been trained.
#'
#' As input parameter it needs the output of \code{\link{descent_generalized_fista}}. If fista has been evoked
#' with save_all_tweaks = T, and a test set is available, the loss-function can be evaluated for all tweak vectors,
#' and both training and test loss curves will be plotted (=> detect overfitting).\cr
#' For an example see
#' \code{\link{descent_generalized_fista}}
#'
#' @param fista.output list, with "Convergence" entry as return by the descent_generalized_fista function
#' @param test.set list with "mixtures" matrix, and "quantity" matrix as returned by mix.samples or mix.samples.jitter function
#' @param X.matrix numeric matrix, reference matrix of the DTD problem
#' @param main string, additionally title
#' @param EVAL.FUN function, that takes a single input (the tweak vector of the fista.output)
#'
#' @return ggplot object
#' @export
#'
#' @import ggplot2
#' @import reshape2
#'
ggplot_convergence <- function(fista.output,
                               test.set = NA,
                               EVAL.FUN,
                               main = ""){
  # convergence can be plotted for training AND test set, if:
  #   - a test set is provided
  #   - the X.matrix is provided
  #   - fista.output has "History" entry
  # otherwise only trainings convergence will be plotted
  if(! (is.na(test.set) || is.null(fista.output$History))){
    # calculate convergence for every saved tweak_vec:
    cor.in.test <- c()
    for(l.iteration in 1:length(fista.output$Convergence)){
    cor.in.test <- c(cor.in.test,
                     EVAL.FUN(tweak = fista.output$History[, l.iteration])
                     )
    }
    # build a data.frame holding training, test and iter
    convergence <- data.frame("training" = fista.output$Convergence,
                              "test" = cor.in.test,
                              "iter" = 1:length(fista.output$Convergence))


  }else{
    # if the convergence can not be calculated on the test set,
    # the data.frame only holds training and iter:
    convergence <- data.frame("training" = fista.output$Convergence,
                              "iter" = 1:length(fista.output$Convergence))
  }
  # melt convergence:
  convergence.melt <- melt(convergence, id.var = "iter")
  # set title:
  tit <- paste0("Loss-function curve during FISTA optimization \n", main)

  # build the ggplot object
  pic <- ggplot(convergence.melt, aes_string(x = "iter", y = "value", col = "variable")) +
                geom_point() + xlab("Iteration") + ylab("Loss-Function") +
                ggtitle(tit)

  return(pic)
}
