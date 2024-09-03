#' Summary function for aftsem package
#'
#' Provides a summary of an aftsem model fit, including the model call,
#' residuals, initial and final coefficient estimates, method, convergence status,
#' number of iterations, number of observations, percent of censored observations,
#' and if available, the estimated covariance matrix of the coefficients, standard
#' deviations, z-values, and p-values for a Wald test.
#'
#' @param object An object of aftsemfit
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class `summaryaftsem` that contains summary information of the
#' fitted aftsem model.
#'
#' @export
#'
summary.aftsem <- function(object,...)
{
  aft <- object
  summobject <- NULL # result of summary
  
  if(is.null(aft$beta))
  {
    return(object) # same as print if algorithm didnt converged or dont have distribution estimation
  }
  
  if(aft$method == "gehan-poly" || aft$method == "buckley")
  {
    return(object)
  }
  
  if(aft$method == "gehan" || aft$method == "jin")
  {
    if(aft$sampling.used == FALSE)
    {
      return(object)
    }
  }
  
  if(aft$method == "gehan-heller" && aft$heller.estimated == FALSE)
  {
    return(object)
  }
  
  summobject$call <- aft$call
  summobject$resid <- aft$resid
  summobject$binit <- aft$betafirst
  summobject$method <- aft$method
  summobject$converged <- aft$converged
  summobject$beta <- aft$beta
  summobject$iters <- aft$iters
  summobject$nobs <- aft$nobs
  summobject$intercept <- aft$intercept
  one_percent = aft$nobs/100
  summobject$percent_censored <- round(aft$censored / one_percent, digits = 2)
  
  
  if(aft$sampling.used)
  {
    b_overline <- apply(aft$beta_star, 1, mean)
    diff_matrix <- aft$beta_star - b_overline
    summobject$cov <- diff_matrix %*% t(diff_matrix)/(aft$resample - 1) # estimation of covariance matrix by sample variance
    
    #####################################################################################
    
    summobject$sd <- sqrt(diag(summobject$cov)) # getting standard deviation
    summobject$zvalue <- summobject$beta/summobject$sd # Wald test
    summobject$pvalue <- (1 - pnorm(abs(summobject$zvalue))) * 2
    summobject$coef_matrix <- matrix(c(aft$beta,summobject$sd,summobject$zvalue,summobject$pvalue),ncol = 4, dimnames = list(aft$cnames,c("|Estimate|", "|Std. Error|", "|Z value|", "|Pr(>|Z|)|")))
  }
  
  
  # for now leave it like this, but later it would need some refactoring
  
  if(aft$method == "gehan-heller" && aft$heller.estimated)
  {
    summobject$cov = aft$covariance # already computed from hellers estimation
    
    #####################################################################################
    
    summobject$sd <- sqrt(diag(summobject$cov)) # getting standard deviation
    summobject$zvalue <- summobject$beta/summobject$sd # Wald test
    summobject$pvalue <- (1 - pnorm(abs(summobject$zvalue))) * 2
    summobject$coef_matrix <- matrix(c(aft$beta,summobject$sd,summobject$zvalue,summobject$pvalue),ncol = 4, dimnames = list(aft$cnames,c("|Estimate|", "|Std. Error|", "|Z value|", "|Pr(>|Z|)|")))
  }
  
  
  class(summobject) <- "summaryaftsem"
  summobject
}

#' Print method for objects of class `summaryaftsem`
#'
#' @param x An object of class `summaryaftsem`
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The function prints object `summaryaftsem`
#' @export

print.summaryaftsem <- function(x,...)
{
  cat("AFT Semiparametric Model Summary\n")
  cat("================================\n\n")
  
  cat(sprintf("Model Call: %s\n", deparse(x$call)))
  cat("\n")
  
  cat("Used parameter estimation method: \n")
  print(x$method)
  cat("\n")
  
  cat("Convergence Status: ")
  cat(ifelse(x$converged, "Converged\n", "Not Converged\n"))
  cat("\n")
  
  printCoefmat(x$coef_matrix, digits = 3, signif.stars = TRUE, signif.legend = TRUE, P.values = TRUE, has.Pvalue = TRUE, na.print = '-')
  
  if(!is.null(x$intercept))
  {
    cat("Estimated Intercept:\n")
    print(x$intercept)
  }
  
  if(!is.null(x$iters))
  {
    cat(sprintf("Number of Iterations: %d\n", x$iters))
  }
  cat(sprintf("Number of Observations: %d\n", x$nobs))
  cat(sprintf("Percent of Censored Observations: %f\n", x$percent_censored))
  cat("\n")
  
  
}
