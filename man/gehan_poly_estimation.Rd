% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aftsem_gehan_poly.R
\name{gehan_poly_estimation}
\alias{gehan_poly_estimation}
\title{Estimation of Regression Parameters from Smoothed Gehan Function}
\usage{
gehan_poly_estimation(y, Z, delta, binit, epsilon, optimx.alg, use.grad)
}
\arguments{
\item{y}{A numeric vector of the response variable, survival times.}

\item{Z}{A matrix of covariates.}

\item{delta}{A censoring indicator vector where 1 indicates an uncensored observation and 0 indicates a censored observation.}

\item{binit}{Initial values for the beta coefficients.}

\item{epsilon}{Smoothing parameter.}

\item{optimx.alg}{Optimalization algorithm that will be used (see optimx package documentation for more information)}

\item{use.grad}{Indicator wheter numerical or excact gradient will be used, default is FALSE == numerical}
}
\value{
A list containing:
- `BETA`: The estimated beta coefficients.
- `RESID`: The residuals from the model fit.
- `ITERS`: The number of iterations performed during optimization.
}
\description{
Estimates regression parameters by optimizing a smoothed version of Gehan's statistic.
}
\details{
The `gehan_poly_estimation` function performs estimation of regression parameters 
by minimizing the smoothed Gehan's loss function.
}
