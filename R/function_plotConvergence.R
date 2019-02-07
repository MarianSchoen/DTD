#' Plot loss curve
#'
#' The "ggplot_convergence" function uses ggplot2 and reshape to visualize the decrease of the loss-function after
#' a model has been trained.
#'
#' As input parameter it needs the output of \code{\link{train_correlatio_model}}, \code{\link{DTD_cv_lambda}},
#' or\code{\link{descent_generalized_fista}}.
#' If fista has been evoked with save_all_tweaks = T, and a evaluation function is available,
#' the loss-function can be evaluated for all intermediate tweak vectors,
#' and both training and test loss curves will be plotted (=> detect overfitting).\cr
#' For an example see section "Visualization of learn curve" in the package vignette `browseVignettes("DTD")`
#'
#' @param DTD.model list, as returned by the train_correlation_model function.
#' @param EVAL.FUN function, that takes a single input (the tweak vector of the fista.output),
#' and returns the loss function on the test set. Default is NA
#' @param main string, additionally title (default "")
#'
#' @return ggplot object
#' @export
#'
#' @import ggplot2
#' @import reshape2
#'
ggplot_convergence <- function(DTD.model,
                               EVAL.FUN = NA,
                               main = "") {
  # convergence can be plotted for training AND test set, if:
  #   - the complete model of 'DTD.model' has the "History" entry
  #   - if EVAL.FUN returns a numeric value for every entry in "History"
  # otherwise only trainings convergence will be plotted

  # as a DTD.model either a list, or only the tweak vector can be used:
  if(is.list(DTD.model)){
    if("best.model" %in% names(DTD.model)){
      fista.output <- DTD.model$best.model
    }else{
      if(all(c("Tweak", "Convergence") %in% names(DTD.model))){
        stop("ggplot_convergence: DTD.model does not fit")
      }else{
        fista.output <- DTD.model
      }
    }
  }else{
    stop("ggplot_convergence: DTD.model is not a list")
  }


  # Test if EVAL.FUN has been set:
  test.eval.fun <- !suppressWarnings(is.na(EVAL.FUN))
  # If EVAl.FUN has been set, therefore "test.eval.fun" is TRUE, test if it returns a numeric value:
  useable.eval.fun <- ifelse(test.eval.fun, is.numeric(EVAL.FUN(fista.output$Tweak)), FALSE)

  # If the eval.fun is "useable", and there is a History, then test can be evaluated:
  if (useable.eval.fun && !is.null(fista.output$History)) {
    # calculate convergence for every saved tweak_vec:
    cor.in.test <- c()
    for (l.iteration in 1:length(fista.output$Convergence)) {
      cor.in.test <- c(
        cor.in.test,
        EVAL.FUN(tweak = fista.output$History[, l.iteration])
      )
    }
    # build a data.frame holding training, test and iter
    convergence <- data.frame(
      "training" = fista.output$Convergence,
      "test" = cor.in.test,
      "iter" = 1:length(fista.output$Convergence)
    )
  } else {
    # if the convergence can not be calculated on the test set,
    # the data.frame only holds training and iter:
    convergence <- data.frame(
      "training" = fista.output$Convergence,
      "iter" = 1:length(fista.output$Convergence)
    )
  }
  # melt convergence:
  convergence.melt <- reshape2::melt(convergence, id.var = "iter")
  # set title:
  tit <- paste0("Loss-function curve during FISTA optimization \n", main)

  # build the ggplot object
  pic <- ggplot2::ggplot(convergence.melt, aes_string(x = "iter", y = "value", col = "variable")) +
    ggplot2::geom_point() +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab("Loss-Function") +
    ggplot2::ggtitle(tit)

  return(pic)
}
