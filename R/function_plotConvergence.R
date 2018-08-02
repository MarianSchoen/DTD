#' Plot loss curve
#'
#' The "ggplot_correlation" function uses ggplot2 and reshape to visualize the decrease of the loss-function after
#' a model has been trained.
#'
#' As input parameter it needs the output of \code{\link{descent_generalized_fista}}. If fista has been evoked
#' with save_all_tweaks = T, and a test set is available, the loss-function can be evaluated for all tweak vectors,
#' and both training and test loss curves will be plotted (=> detect overfitting).\cr
#' The following example needs data generation, and evoking the FISTA algorithm as shown in
#' \code{\link{descent_generalized_fista}}
#'
#' @param fista.output list, with "Convergence" entry as return by the descent_generalized_fista function
#' @param test.set list with "mixtures" matrix, and "quantity" matrix as returned by mix.samples or mix.samples.jitter function
#' @param X.matrix numeric matrix, reference matrix of the DTD problem
#'
#' @return
#' @export
#'
#' @import ggplot2
#' @import reshape2
#'
ggplot_correlation <- function(fista.output, test.set = NA, X.matrix = NA){
  # correlation can be plotted for training and test set, if:
  #   - a test set is provided
  #   - the X.matrix is provided
  #   - fista.output has "History" entry
  # otherwise only trainings correlation will be plotted
  if(! (is.na(test.set) || is.na(X.matrix) || is.null(fista.output$History))){
    cor.in.test <- c()
    for(l.iteration in 1:length(fista.output$Convergence)){
    cor.in.test <- c(cor.in.test, evaluate_cor(X = X.matrix,
                                                Y = test.set$mixtures,
                                                C = test.set$quantities,
                                                tweak = fista.output$History[, l.iteration])
                      )
    }
    convergence <- data.frame("trainig" = fista.output$Convergence,
                              "test" = cor.in.test,
                              "iter" = 1:length(fista.output$Convergence))

    convergence.melt <- melt(convergence, id.var = "iter")
  }else{
    convergence <- data.frame("trainig" = fista.output$Convergence,
                              "iter" = 1:length(fista.output$Convergence))

    convergence.melt <- melt(convergence, id.var = "iter")
  }

  pic <- ggplot(convergence.melt, aes(x=iter, y = value, col = variable)) +
                geom_point() + xlab("Iteration") + ylab("Loss-Function") +
                ggtitle("Loss-function curve during FISTA optimization")

  return(pic)
}
