#' Semi-parametric AFT Model Fitting
#'
#' @description
#' Fits a semi-parametric accelerated failure time (AFT) model to the provided data using various methods.
#'
#' @param Z A matrix of covariates.
#' @param y A vector of the response variable, typically survival times.
#' @param delta A censoring indicator vector where 1 indicates an uncensored observation and 0 indicates a censored observation.
#' @param betafirst The initial estimate of the beta coefficients.
#' @param method The method of estimation to use, one of "buckley", "gehan", "jin", or "gehan-poly".
#' @param control A list of control parameters including `eps` for convergence criterion and `maxiter` for the maximum number of iterations.
#' @param intercept Logical; if TRUE, include an intercept in the model.
#' @param resample The number of resamples to use for Monte Carlo estimation of variance; relevant for certain methods only.
#' @param nobs The number of observations in the data.
#'
#' @details
#' The `aftsem_fit` function provides a way to fit a semi-parametric AFT model
#' to survival data with potential RIGHT censoring. Depending on the chosen method, 
#' different estimation techniques are used, such as Buckley-James or Gehan's method.
#' If resampling is required for the method, the function will generate resamples from an exponential distribution.
#'
#' @return
#' Returns a list object of class "aftsem" containing the following components:
#' - `converged`: Logical indicating if the fitting procedure converged.
#' - `beta`: The estimated beta coefficients.
#' - `iters`: The number of iterations performed.
#' - `resid`: The residuals from the model fit. NOT THE MARTINGALE RESIDUALS
#' - `sampling.used`: Logical indicating if sampling was used.
#' - `intercept`: The estimated intercept, included if `intercept = TRUE`.
#' - `beta_star`: The beta coefficients estimated for each resample, included if resampling was used.
#' - `fe`: Number of calls of function in minimalization proccess (only available for gehan-poly and gehan-heller method)
#' - `covariance` Covariance matrix (only available for gehan-heller method)
aftsem_fit <- function (Z,
                        y,
                        delta,
                        betafirst,
                        method,
                        control,
                        intercept,
                        resample,
                        nobs)
{
  
  
  fit.result <- NULL
  beta_star <- NULL
  R <- NULL
  sampling.used <- FALSE
  
  
  if(method %in% c("gehan","jin") && resample >= 2) sampling.used <- TRUE
  
  if(sampling.used)
  {
    # set.seed(42) you can uncomment this line if you need to debug
    # generate matrix for variance estimation, we use exponentional distribution
    R <- matrix(rexp(nobs*resample), ncol=resample)
  }
  
  # dummy
  rs <- matrix(rep(1,nobs), ncol=1)
  
  #
  # Feel free to add more methods, you may want to change the fit.result and fit lists
  #
  
  switch (method,
    "buckley" =
      {
        fit.result <- estimate_buckley(y,Z,delta,betafirst,control$eps,control$maxiter)
      },
    "gehan" = 
      {
        fit.result <- gehan_estimation(y,Z,delta,rs,control$quantile.method)
        if (sampling.used) beta_star <- gehan_estimation(y,Z,delta,R,control$quantile.method)
      },
    "jin" =
      {
        if(sampling.used) beta_star <- gehan_estimation(y,Z,delta,R,control$quantile.method)
        if(sampling.used == FALSE)
        {
          beta_star <- rs #dummy
          R <- rs #dummy
        }
        fit.result <- estimate_jin(y,Z,delta,betafirst,control$eps,control$maxiter,R,beta_star,sampling.used)
        if(sampling.used) beta_star <- fit.result$DISTRIBUTION
      },
    "gehan-poly" = 
      {
        fit.result <- gehan_poly_estimation(y,Z,delta,betafirst,control$gehan_eps, control$optimx.alg, control$use.grad)
      },
    "gehan-heller" = 
      {
        fit.result <- gehan_heller_estimation(y,Z,delta,betafirst,control$optimx.alg, control$variance.estimation, control$use.grad)
      },
      {
        stop("Unknown method specified")
      }
  )
  
  
  
  fit <- list(converged = fit.result$CONVERGED, beta = fit.result$BETA, iters = fit.result$ITERATIONS, resid = fit.result$RESID,sampling.used = sampling.used, that = fit.result$THAT)
  class(fit) <- "aftsem"
  fit$heller.estimated = control$variance.estimation
  if(intercept)
  {
    fit$intercept = fit.result$INTERCEPT
  }
  if(sampling.used)
  {
    fit$beta_star <- beta_star
  }
  
  if(method %in% c("gehan-heller","gehan-poly"))
  {
    fit$fe <- fit.result$FN_CALLS
  }
  if(method == "gehan-heller" && control$variance.estimation == TRUE)
  {
    fit$covariance = fit.result$COVARIANCE
  }

  return(fit)
}
