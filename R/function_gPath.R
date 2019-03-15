#' Visualize the path of each \eqn{g_i} over all iterations
#'
#' With this function the path of each \eqn{g_i} over all iterations can be plotted.
#' Notice, that if there are many genes in your model, it may be hard to distinguish between each path.
#' As a solution the parameter "number.pics" can be set to a higher integer.
#' Then, the visualization will be split into more pictures.
#' In each picture all \eqn{g_i} get collected that end up in the same quantile range.
#' E.g if you split into 3 pictures, the first picture includes all genes that result into the
#' quantile range from 0% Qu to 33% Qu of all g.
#' There are parameters (G.TRANSFORM.FUN and ITER.TRANSFORM.FUN) to transform the g vector, and iteration number.
#' These might help to make the plot more understandable, e.g. if the distribution of the g vector is dominated by
#' same outliers, applying a log transformation might help.
#' In most of the scenarios the major changes in the g vector occur in the early iterations. Focus on this part via a log transformation.
#' For an example see section "g-Path" in the package vignette `browseVignettes("DTD")`
#'
#' @param DTD.model : list as returned by \code{\link{train_correlatio_model}}, \code{\link{DTD_cv_lambda}},
#' or \code{\link{descent_generalized_fista}}
#' @param number.pics : integer, into how many pictures should the resutlt be split. (defaults to 3)
#' @param G.TRANSFORM.FUN : function, that expects a vector of numerics, and returns a vector of the same length.
#' Will be applied on fista.output$Tweak. Set G.TRANSFORM.FUN to identity if no transformation is required.
#' If you change G.TRANSFORM.FUN don't forget to adjust the y.lab parameter. (default is identity)
#' @param ITER.TRANSFORM.FUN : function, that expects a vector of numerics, and returns a vector of the same length.
#' Will be applied on the iteration/x.axis of the plot. Set ITER.TRANSFORM.FUN to identity if no transformation is required.
#' If you change ITER.TRANSFORM.FUN don't forget to adjust the x.lab parameter (Default: log10)
#' @param y.lab string, used as y label on the plot (Default: "g)
#' @param x.lab string, used as x label on the plot (default is "log10(iteration)")
#' @param plot.legend logical, should the legend be plotted? Notice that the legend will
#' be plotted in a additional figure, and can be visualized via grid::grid.draw (Default: FALSE)
#' @param subset vector of strings, or NA that match the rownames of fista.output$History.
#' Only these entries will be visualized. If set to NA, all entries will be used. (Default: NAA)
#' @param main string, used as title of the plot (Default: "")
#'
#' @import ggplot2
#' @import reshape2
#' @return list, with "gPath" entry. "gPath" will be a ggplot object.
#' Depending on "plot.legend" the list has a second entry named "legend". "legend" will be a grid object.
#' @export
#'
ggplot_gpath <- function(DTD.model,
                         number.pics = 3,
                         G.TRANSFORM.FUN = DTD::identity,
                         ITER.TRANSFORM.FUN = log10,
                         y.lab = "g",
                         x.lab = "log10(iteration)",
                         subset = NA,
                         main = "",
                         plot.legend = FALSE) {

  # safety check: number.pics
  test <- test_integer(number.pics,
                       output.info = c("ggplot_gpath", "number.pics"),
                       min = 1,
                       max = Inf)
  # end -> number.pics

  # safety check: y.lab
  useable.ylab <- try(as.character(y.lab), silent = TRUE)
  if(any(grepl(x = useable.ylab, pattern = "Error"))){
    stop("In ggplot_gpath: provided 'y.lab' can not be used as.character.")
  }
  # end -> y.lab

  # safety check: x.lab
  useable.ylab <- try(as.character(x.lab), silent = TRUE)
  if(any(grepl(x = useable.ylab, pattern = "Error"))){
    stop("In ggplot_gpath: provided 'x.lab' can not be used as.character.")
  }
  # end -> x.lab

  # safety check: main
  useable.ylab <- try(as.character(main), silent = TRUE)
  if(any(grepl(x = useable.ylab, pattern = "Error"))){
    stop("In ggplot_gpath: provided 'main' can not be used as.character.")
  }
  # end -> main



  # safety check: plot.legend
  test <- test_logical(test.value = plot.legend,
                       output.info = c("ggplot_gpath", "plot.legend"))
  # end -> plot.legend

  # for gPath, the following elements are needed:
  # - 'History' of learning

  # Either it is provided as a list ...
  if(is.list(DTD.model)){
    if("best.model" %in% names(DTD.model)){
      fista.output <- DTD.model$best.model
    }else{
      if(all(c("Tweak", "History") %in% names(DTD.model))){
        stop("In ggplot_gpath: 'DTD.model' can not be used (provide a DTD.model with 'History' entry)")
      }else{
        fista.output <- DTD.model
      }
    }
    f.history <- fista.output$History
    tweak <- fista.output$Tweak
  }else{
    # ... or only the History matrix is provided:
    if(!is.matrix(DTD.model)){
      stop("In ggplot_gpath: 'DTD.model' can not be used (provide a DTD.model with 'History' entry)")
    }else{
      f.history <- DTD.model
      tweak <- f.history[, ncol(f.history)]
    }
  }
  # safety check: G.TRANSFORM.FUN
  useable.g.trans.fun <- try(G.TRANSFORM.FUN(tweak), silent = TRUE)
  if(any(grepl(x = useable.g.trans.fun, pattern = "Error")) || any(!is.numeric(useable.g.trans.fun))){
    stop("In ggplot_gpath: 'G.TRANSFORM.FUN' does not return numeric vector.")
  }
  # end -> main

  # safety check: ITER.TRANSFORM.FUN
  useable.iter.trans.fun <- try(ITER.TRANSFORM.FUN(tweak), silent = TRUE)
  if(any(grepl(x = useable.iter.trans.fun, pattern = "Error")) || any(!is.numeric(useable.iter.trans.fun))){
    stop("In ggplot_gpath: 'ITER.TRANSFORM.FUN' does not return numeric vector.")
  }
  # end -> main

  # if:
  # - subset is not na,
  # - any subset is within rownames
  if (!all(is.na(subset)) && any(subset %in% rownames(f.history))) {
    subset <- subset[subset %in% rownames(f.history)]
    f.history <- f.history[subset, , drop = FALSE]
    tweak <- tweak[subset]
  } else {
    if (!all(is.na(subset))) {
      message("In ggplot_gpath: subset could not be used, therefore complete tweak, and history will be used\n")
    }
  }

  # We start by calculating in which quantile range each gene falls:
  # Therefore, we calculate how many quantile ranges are necessary (depending on the number.pics parameter)
  pic.sep <- as.numeric(format(x = seq(0, 1, length.out = (number.pics + 1))[2:(number.pics + 1)], digits = 2))
  # Next we calculate the value of the quantiles in our g vector ...
  quantile.values <- sapply(X = pic.sep, FUN = stats::quantile, x = abs(tweak))
  # ... and name them without the "%" sign
  names(quantile.values) <- gsub("%", "", names(quantile.values))

  # Now we know the quantile values, next we have to test in which quantile ranges our value fall.
  # Therefore the following function helps.
  # It takes a numeric value x, and returns the position of the first quantile value that is below it
  quantile.apply.function <- function(x, values = quantile.values) {
    winner <- which(values >= abs(x))[1]
    return(winner)
  }
  # apply the function ...
  quantile.per.gene <- sapply(tweak, quantile.apply.function)
  # ... set names
  names(quantile.per.gene) <- rownames(f.history)

  # For easy visualization with the ggplot2 package we need the f.history matrix in long format:
  f.h.melt <- reshape2::melt(f.history,
    varnames = c("geneName", "iteration"),
    value.name = "g"
  )

  # add the "q"unatile "p"er "g"ene information (with the names, not the positions!)
  f.h.melt$qpg <- factor(as.numeric(names(quantile.values[quantile.per.gene[f.h.melt$geneName]])))
  # Reset the levels to more interpretable
  levels(f.h.melt$qpg) <- paste0("below ", levels(f.h.melt$qpg), "% Quantile")

  # Transform g and iter
  f.h.melt$g <- G.TRANSFORM.FUN(f.h.melt$g)
  f.h.melt$iteration <- ITER.TRANSFORM.FUN(f.h.melt$iteration)


  # Plot the picture (notice that this plot is with the legend!)
  pics <- ggplot2::ggplot(
    f.h.melt,
    aes_string(x = "iteration", y = "g", group = "geneName", colour = "geneName")
  ) +
    ggplot2::geom_line() +
    ggplot2::ylab(y.lab) +
    ggplot2::xlab(x.lab) +
    ggplot2::ggtitle(main) +
    ggplot2::facet_grid(. ~ qpg)

  ret <- list()
  # Store the picture WITHOUT the legend
  ret[["gPath"]] <- pics + theme(legend.position = "none")

  # Only if required, extract the legend from "pics" and provide as a entry in ret
  if (plot.legend) {
    tmp <- ggplot2::ggplot_gtable(ggplot_build(pics))
    tmp.leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[tmp.leg]]
    ret[["legend"]] <- legend
  }
  # return ret (including the picture without legend, and the legend if required)
  return(ret)
}
