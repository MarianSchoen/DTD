% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_norms.R
\name{n1normed}
\alias{n1normed}
\title{Norm a vector such that ||x||_1 is equal to the length of x}
\usage{
n1normed(x)
}
\arguments{
\item{x}{numeric vector}
}
\value{
numeric vector with same length as x
}
\description{
Norm a vector such that ||x||_1 is equal to the length of x
}
\examples{
set.seed(1)
same.values <- rnorm(n = 10)
norm(as.matrix(same.values), type = "O")
normed.values <- n1normed(same.values)
norm(as.matrix(normed.values), type = "O")
}