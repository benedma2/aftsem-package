% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aftsem_print.R
\name{print.aftsem}
\alias{print.aftsem}
\title{Print method for aftsem xs}
\usage{
\method{print}{aftsem}(x, ...)
}
\arguments{
\item{x}{An x of class "aftsem", typically the result of a call to `aftsem_fit`.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
The function is called for its side effect, which is printing the summary to the console. It invisibly returns NULL.
}
\description{
Prints a summary of an aftsem model fit x.
}
\details{
The `print.aftsem` method provides a user-friendly summary of the model fit, including
the method used for parameter estimation, convergence status, estimated parameters, number of iterations, and
the percentage of censored observations.
}
\seealso{
\code{\link{aftsem_fit}} for model fitting.
}
