% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_eval_cor.R
\name{evaluate_cor}
\alias{evaluate_cor}
\title{Evaluate correlation}
\usage{
evaluate_cor(
  X.matrix = NA,
  new.data,
  true.compositions,
  DTD.model,
  estimate.c.type
)
}
\arguments{
\item{X.matrix}{numeric matrix, with features/genes as rows,
and cell types as column. Each column of X.matrix is a reference
expression profile. A trained DTD model includes X.matrix, it has been
trained on. Therefore, X.matrix should only be set, if the 'DTD.model'
is not a DTD model.}

\item{new.data}{numeric matrix with samples as columns,
and features/genes as rows.}

\item{true.compositions}{numeric matrix with cells as rows,
and mixtures as columns. In the formula above named true_C.
Each row of C holds the distribution of the cell over all mixtures.}

\item{DTD.model}{either a numeric vector with length of nrow(X), or a list
returned by \code{\link{train_deconvolution_model}},
\code{\link{DTD_cv_lambda_cxx}}, or \code{\link{descent_generalized_fista}}.
In the equation above the DTD.model provides the vector g.}

\item{estimate.c.type}{string, either "non_negative", or "direct".
Indicates how the algorithm finds the solution of
\eqn{arg min_C ||diag(g)(Y - XC)||_2}.
\itemize{
   \item If 'estimate.c.type' is set to "direct",
 there is no regularization (see \code{\link{estimate_c}}),
   \item if 'estimate.c.type' is set to "non_negative",
 the estimates "C" must not be negative (non-negative least squares)
(see (see \code{\link{estimate_nn_c}}))
}}
}
\value{
float, value of the Loss function
}
\description{
The loss-function learning digital tissue deconvolution finds a vector g
which minimizes the Loss function L\cr
\deqn{L(g) = - \sum cor(true_C, estimatd_C(g))}
The evaluate_cor function returns the value of the Loss function.
}
\examples{
library(DTD)
random.data <- generate_random_data(
    n.types = 10,
    n.samples.per.type = 150,
    n.features = 250,
    sample.type = "Cell",
    feature.type = "gene"
    )

# normalize all samples to the same amount of counts:
normalized.data <- normalize_to_count(random.data)

# extract indicator list.
# This list contains the Type of the sample as value, and the sample name as name
indicator.list <- gsub("^Cell[0-9]*\\\\.", "", colnames(random.data))
names(indicator.list) <- colnames(random.data)

# extract reference matrix X
# First, decide which cells should be deconvoluted.
# Notice, in the mixtures there can be more cells than in the reference matrix.
include.in.X <- paste0("Type", 2:7)

percentage.of.all.cells <- 0.2
sample.X <- sample_random_X(
    included.in.X = include.in.X,
    pheno = indicator.list,
    expr.data = normalized.data,
    percentage.of.all.cells = percentage.of.all.cells
    )
X.matrix <- sample.X$X.matrix
samples.to.remove <- sample.X$samples.to.remove
remaining.mat <- normalized.data[, -which(colnames(normalized.data) \%in\% samples.to.remove)]

indicator.remain <- indicator.list[names(indicator.list) \%in\% colnames(remaining.mat)]
training.data <- mix_samples(
    expr.data = remaining.mat,
    pheno = indicator.remain,
    included.in.X = include.in.X,
    n.samples = 500,
    n.per.mixture = 100,
    verbose = FALSE
    )

 start.tweak <- rep(1, nrow(X.matrix))
 sum.cor <- evaluate_cor(
   X.matrix = X.matrix
   , new.data = training.data$mixtures
   , true.compositions = training.data$quantities
   , DTD.model = start.tweak
   , estimate.c.type = "direct"
   )

 rel.cor <- sum.cor/ncol(X.matrix)
cat("Relative correlation: ", -rel.cor, "\n")
}
