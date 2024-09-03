#' Estimation of Regression Parameters from Smoothed Gehan Function
#'
#' @description
#' Estimates regression parameters by optimizing a smoothed version of Gehan's statistic.
#'
#' @param y A numeric vector of the response variable, survival times.
#' @param Z A matrix of covariates.
#' @param delta A censoring indicator vector where 1 indicates an uncensored observation and 0 indicates a censored observation.
#' @param binit Initial values for the beta coefficients.
#' @param epsilon Smoothing parameter.
#' @param optimx.alg Optimalization algorithm that will be used (see optimx package documentation for more information)
#' @param use.grad Indicator wheter numerical or excact gradient will be used, default is FALSE == numerical
#'
#' @details
#' The `gehan_poly_estimation` function performs estimation of regression parameters 
#' by minimizing the smoothed Gehan's loss function.
#'
#' @return
#' A list containing:
#' - `BETA`: The estimated beta coefficients.
#' - `RESID`: The residuals from the model fit.
#' - `ITERS`: The number of iterations performed during optimization.
#'
gehan_poly_estimation <- function(y,
                                  Z,
                                  delta,
                                  binit,
                                  epsilon,
                                  optimx.alg,
                                  use.grad)
{
  n <- nrow(Z)
  p <- ncol(Z)
  
  res <- NULL
  binit <- as.vector(binit) #because gehan estimation
  
  #
  # Sensitive on initial beta estimate, better with numerical estimation
  #
  
  f_epsilon_grad <- function(beta)
  {
    cppverg <- compute_f_epsilon_grad(beta,y,Z,delta,epsilon)
    return(cppverg)
  }
  
  
  f_epsilon <- function(beta)
  {
    cppver <- compute_f_epsilon(beta,y,Z,delta,epsilon)
    return(cppver)
  }
  
  
  optim.res <- NULL
  
  if(use.grad == FALSE)
  {
    optim.res <- optimx(par = binit,
                        fn = f_epsilon,
                        #gr = f_epsilon_grad,
                        lower = -Inf,
                        upper = Inf,
                        method = optimx.alg)
  }
  
  if(use.grad)
  {
    optim.res <- optimx(par = binit,
                        fn = f_epsilon,
                        gr = f_epsilon_grad,
                        lower = -Inf,
                        upper = Inf,
                        method = optimx.alg)
  }
  

  beta <- optim.res[1:p]
  beta <- unlist(beta)
  beta <- as.vector(beta)
  
  w <- rep(1, n)
  w <- as.vector(w)
  e <- y - Z %*% beta
  km <- km_e(e, delta, w)
  y_star <- delta * y + (1 - delta) * (km + Z %*% beta)
  resid <- y_star - Z%*%beta
  
  beta <- as.matrix(beta)
  
  
  code <- optim.res$convcode
  bool.code <- FALSE
  
  if(code == 0)
  {
    bool.code <- TRUE
  }
  
  res = list(BETA = beta, RESID = resid, ITERS = optim.res$niter, FN_CALLS = optim.res$fevals, CONVERGED = bool.code, THAT = y_star)
  
  return(res)
}