#' Gehan-Heller Estimation of regression parameters
#'
#' @param y Numeric vector of survival times or times to event/censoring.
#' @param Z Numeric matrix of covariates with observations in rows and covariates in columns.
#' @param delta Numeric vector indicating censoring, with 1 for an event and 0 for censored observations.
#' @param binit Numeric vector or matrix for initial estimates of regression coefficients.
#' @param optimx.alg Optimalization algorithm that will be used (see optimx package documentation for more information)
#' @param variance.estimation If covariance matrix will be estimated
#' @param use.grad Indicator wheter numerical or excact gradient will be used, default is FALSE == numerical
#' 
#' Covariance estimation is programmed but not tested!
#' 
#' @note The recommend use is with numerical aproximation of gradient. The true gradiet can be sensitive for initial beta values (binit).
#' For Covariance estimation please set the variance.estimation in control list to TRUE.
#' 
#' @return A list containing the estimated regression coefficients (`BETA`), residuals (`RESID`),
#'         and the number of iterations taken by the optimization routine (`ITERS`).
#'

gehan_heller_estimation <- function(y,
                                    Z,
                                    delta,
                                    binit,
                                    optimx.alg,
                                    variance.estimation,
                                    use.grad
                                    )
{
  n <- nrow(Z)
  p <- ncol(Z)
  
  binit <- as.vector(binit) #because gehan estimation
  res <- NULL
  
  inx <- which(delta==1)
  

  e_s <- y[inx] - as.matrix(Z[inx,]) %*% binit
  a <- sd(e_s) * n^(-(1/5))
  
 #
 # Sensitive on initial beta estimate, better with numerical estimation
 #
  
  heller_grad <- function(beta)
  {
   cppverg <- compute_heller_grad(beta,y,Z,delta,a)
   return(cppverg)
  }
  
  
  heller <- function(beta)
  {
    cppver <- compute_heller(beta,y,Z,delta,a)
    return(cppver)
  }
  
  optim.res <- NULL
  
  if(use.grad == FALSE)
  {
    optim.res <- optimx(par = binit,
                        fn = heller,
                        #gr = heller_grad,
                        lower = -Inf,
                        upper = Inf,
                        method = optimx.alg,
                        )
  }
  
  if(use.grad)
  {
    optim.res <- optimx(par = binit,
                        fn = heller,
                        gr = heller_grad,
                        lower = -Inf,
                        upper = Inf,
                        method = optimx.alg,
                        )
  }
  
  
  beta <- optim.res[1:p]
  beta <- unlist(beta)
  beta <- as.vector(beta)
  
  weights <- rep(1, n)
  e <- y - Z %*% beta
  km <- km_e(e, delta, weights)
  y_star <- delta * y + (1 - delta) * (km + Z %*% beta)
  resid <- y_star - Z%*%beta
  
  if(variance.estimation)
  {
    covariance = compute_covariance(delta,Z,e,a)
    covariance <- as.matrix(covariance)
  }
  
  beta <- as.matrix(beta)
  
  code <- optim.res$convcode
  bool.code <- FALSE
  
  if(code == 0)
  {
    bool.code <- TRUE
  }
  
  res = list(BETA = beta, RESID = resid, ITERS = optim.res$niter, FN_CALLS = optim.res$fevals, CONVERGED = bool.code, THAT = y_star)
  
  if(variance.estimation)
  {
    res$COVARIANCE = covariance
  }
  
  
  return(res)
}