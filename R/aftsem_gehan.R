#' Gehan's Estimation for Survival Data
#'
#' This function performs Gehan's estimation of regression parameters proposed by Jin
#'
#' @param y A numeric vector of survival times.
#' @param Z A matrix of covariates
#' @param delta A numeric vector indicating censoring status
#' @param rsmat A resampling matrix 
#' @param init A logical value indicating whether to return the initial fit object
#'   (default is `FALSE`). If `FALSE`, only the coefficients are returned.
#' @param m Method for quantreg optimalization
#'
#' @return If `init = FALSE` and `change == 1`, returns a list with elements `INTERCEPT`,
#'   `RESID`, `ITERS`, `CONVERGED`, `BETA`. Otherwise, returns a matrix of resampled
#'   Gehan estimates.
#'
#' @author Zherzen Jin
#'
#'
#' @note This function is a slightly different version from the original by Zherzen Jin,
#'   part of the now not available `lss` program.
gehan_estimation <- function (y,
                              Z,
                              delta,
                              rsmat,
                              m,
                              init = FALSE)
{
  # number should be larger than beta^T sum_k { sum_l { delta_k * (Z_l - Z_k) } }
  ll <- length(y)
  M <- 1000*(ll)^4
  
  dimnum<-dim(Z)
  
  res <- NULL
  
  rowsZ<-dimnum[1]
  colsZ<-dimnum[2]
  change <- ncol(rsmat)
  
  # handles distribution and estimate
  res.gehan <- matrix(0,nrow = colsZ, ncol = change)
  
  #
  # sum_i { sum_j { delta_i * (y_i - y_j) } }
  #
  y_i <- rep(y,rep(rowsZ,rowsZ)) # each y n*times
  delta_i <- rep(delta,rep(rowsZ,rowsZ)) # each delta n * times
  y_j <- rep(y,rowsZ) # vector of y n*times
  y_all <- delta_i*(y_i-y_j)
  
  #
  # beta^T * sum_k { sum_l { delta_k * (Z_l - Z_k) } }
  #
  z_k <- matrix(rep(as.vector(Z),rep(rowsZ,rowsZ*colsZ)),nrow=rowsZ*rowsZ)
  z_i <- t(matrix(rep(as.vector(t(Z)),rowsZ),nrow=colsZ))
  z_all <- z_k - z_i   
  
  for(i in 1:change)
  {
    resamples <- rep(rsmat[,i],rep(rowsZ,rowsZ))*rep(rsmat[,i],rowsZ)
    z_fin <- z_all*resamples*delta_i
    X <- apply(z_fin,2,sum)
    X <- rbind(z_fin,-X)
    Y <- c(y_all*resamples,M)
    
    fit <- rq(Y ~ X - 1, tau = 0.5, method = m) # BR is Koenker & D`orey algorithm, -1 because we dont want intercept
    res.gehan[,i] <- fit$coef
    
  }
  
  
  if(change == 1 && init == FALSE)
  {
    beta <- as.vector(res.gehan[,1])
    weights <- rep(1, ll)
    e <- y - Z %*% beta
    km <- km_e(e, delta, weights)
    y_star <- delta * y + (1 - delta) * (km + Z %*% beta)
    resid <- y_star - Z%*%beta
    
    res <- list(INTERCEPT = NULL, RESID = resid, ITERS = "used rq", CONVERGED = TRUE, BETA = res.gehan, THAT = y_star)
  }
  else
  {
    res <- res.gehan
  }
  
  return(res)
}