#' Plot loss curve
#'
#' The "ggplot_convergence" function uses ggplot2 and reshape2 to visualize the
#' decrease of the loss-function after a model has been trained.
#'
#' As input parameter it needs the output of
#' \code{\link{train_deconvolution_model}}, an DTD cross validation object,
#' or the output of an FISTA optimization run.
#' If the `DTD.model` includes a 'History' entry, and a `test.data` is available
#' the loss function can be evaluated for each intermediate steps of the
#' optimization. Then, the resulting picture includes two convergence paths,
#' one for the training data, and one for the test.data.\cr
#' For an example see section "Visualization of learn curve" in the package
#' vignette `browseVignettes("DTD")`
#'
#' @param DTD.model either a numeric vector with length of nrow(X), or a list
#' returned by \code{\link{train_deconvolution_model}},
#' \code{\link{DTD_cv_lambda_cxx}}, or \code{\link{descent_generalized_fista}}.
#' @param test.data list, with two entries, a numeric matrix each,
#' named 'mixtures' and 'quantities' For examples see \code{\link{mix_samples}},
#' \code{\link{mix_samples_with_jitter}} or the package vignette
#' `browseVignettes("DTD")`.
#' @param estimate.c.type string, either "non_negative", or "direct".
#' Indicates how the algorithm finds the solution of
#' \eqn{arg min_C ||diag(g)(Y - XC)||_2}.
#' \itemize{
#'    \item If 'estimate.c.type' is set to "direct",
#'  there is no regularization (see \code{\link{estimate_c}}),
#'    \item if 'estimate.c.type' is set to "non_negative",
#'  the estimates "C" must not be negative (non-negative least squares)
#' (see (see \code{\link{estimate_nn_c}}))
#' }
#' @param title string, additionally title
#' @param X.matrix numeric matrix, with features/genes as rows,
#' and cell types as column. Each column of X.matrix is a reference
#' expression profile. A trained DTD model includes X.matrix, it has been
#' trained on. Therefore, X.matrix should only be set, if the 'DTD.model'
#' is not a DTD model.
#' @return ggplot object
#' @export
#'
#' @import ggplot2
#' @import reshape2
#'
ggplot_convergence <- function(
  DTD.model,
  X.matrix = NA,
  test.data = NULL,
  estimate.c.type = "decide.on.model",
  title = "") {
  
  # convergence can be plotted for training AND test set, if:
  #   - the complete model of 'DTD.model' has the "History" entry
  #   - if test.data and reference matrix are provided
  #   - estimate.c.type is a "interpretable" type
  # otherwise only trainings convergence will be plotted

  # as a DTD.model either a list, or only the tweak vector can be used:
  if (is.list(DTD.model)) {
    if ("best.model" %in% names(DTD.model)) {
      fista.output <- DTD.model$best.model
    } else {
      if (!all(c("Tweak", "Convergence") %in% names(DTD.model))) {
        stop("In ggplot_convergence: DTD.model does not fit")
      } else {
        fista.output <- DTD.model
      }
    }
    if ("estimate.c.type" %in% names(DTD.model)){
      estimate.c.type <- DTD.model$estimate.c.type
    }
  } else {
    stop("In ggplot_convergence: DTD.model is not a list")
  }

  # check if there is an History entry
  if(!"History" %in% names(fista.output)){
    if(!is.null(test.data)){
      message("In ggplot_convergence: There is no 'History' entry in the DTD model, therefore convergence can not be shown on 'test.data'.")
    }
    test.data <- NULL
  }

  # check if test.data can be used:
  if(!is.null(test.data)){
    if(is.list(test.data) && length(test.data) == 2){
      if(!all(c("quantities", "mixtures") %in%  names(test.data))){
        message("In ggplot_convergence: entries of 'test.data' must be named 'quantities' and 'mixtures'. Therefore, 'test.data' can not be used. Plotting only training convergence")
        test.data <- NULL
      }else{
        if(!is.matrix(test.data$mixtures)){
          message("In ggplot_convergence: 'test.data$mixtures' is not a matrix. Therefore, 'test.data' can not be used. Plotting only training convergence.")
          test.data <- NULL
        }
        if(!is.matrix(test.data$quantities)){
          message("In ggplot_convergence: 'test.data$quantities' is not a matrix. Therefore, 'test.data' can not be used. Plotting only training convergence")
          test.data <- NULL
        }
      }
    }else{
      message("In ggplot_convergence: 'test.data' must be provided as a list with two entries: 'quantities' and 'mixtures'. Therefore, 'test.data' can not be used. Plotting only training convergence")
      test.data <- NULL
    }
  }

  # the reference matrix can either be included in the DTD.model, or has to be past
  # via the X.matrix argument:
  if(!is.numeric(X.matrix) || any(is.na(X.matrix))){
    if (is.list(DTD.model) && "reference.X" %in% names(DTD.model)) {
        X.matrix <- DTD.model$reference.X
      }else{
        message("In ggplot_convergence: provide X matrix either via X.matrix argument or included in the DTD.model. Plotting only training convergence.")
        # kind of 'around the tree', but with this X will not be used
        test.data <- NULL
      }
  }

  # safety check: title
  title <- test_string(
    test.value = title
    , output.info = c("ggplot_convergence", "title"))
  # end -> title


  # safety check for estimate.c.type:
  test <- test_c_type(
    test.value = estimate.c.type
    , output.info = c("ggplot_convergence", "estimate.c.type"))
  # end -> estimate.c.type

  DTD.evCor.wrapper <- function(
    tweak,
    X = X.matrix,
    data.list = test.data,
    esti.c.type = estimate.c.type
    ) {
    Y <- data.list$mixtures
    C <- data.list$quantities
    loss <- evaluate_cor(
      X.matrix = X,
      new.data = Y,
      true.compositions = C,
      DTD.model = tweak,
      estimate.c.type = esti.c.type
    ) / ncol(X)
    return(loss)
  }

  # check if dimensions fit (I)
  if(!is.null(test.data)){
    if( (! is.null(rownames(fista.output$History))) && 
        all(rownames(fista.output$History) %in% rownames(X.matrix))  ){
      X.matrix <- X.matrix[rownames(fista.output$History), , drop = FALSE]
    }else{
      # print(rownames(fista.output$History))
      # print(rownames(X.matrix))
      message("In ggplot_convergence: rownames('X.matrix') does not fit 'DTD.model' (rownames of History entry). Therefore convergence can not be shown on 'test.data'.")
      test.data <- NULL
    }
  }

  # If there is a History, then test can be evaluated:
  if( !is.null(fista.output$History) &&
      !(is.null(test.data) || is.na(test.data))
      ) {
    # calculate convergence for every saved tweak_vec:
    cor.in.test <- c()

    for (l.iteration in 1:length(fista.output$Convergence)) {
      cor.in.test <- c(
        cor.in.test,
        DTD.evCor.wrapper(
          tweak = fista.output$History[, l.iteration, drop = FALSE]
          )
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
  tit <- paste0("Loss-function curve during FISTA optimization \n", title)

  # build the ggplot object
  pic <- ggplot2::ggplot(convergence.melt, aes_string(x = "iter", y = "value", col = "variable")) +
    ggplot2::geom_point() +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab("Loss-Function") +
    ggplot2::ggtitle(tit)

  return(pic)
}
