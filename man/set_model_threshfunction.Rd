% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface_cxx.R
\name{set_model_threshfunction}
\alias{set_model_threshfunction}
\title{set_model_threshfunction}
\usage{
set_model_threshfunction(model, threshfnname)
}
\arguments{
\item{model}{input model, as constructed by, e.g., empty_model()}

\item{threshfnname}{the desired value for the threshold function. May currently be only "softmax".}
}
\value{
the updated model with the thresh function set.
}
\description{
set_model_threshfunction
}
