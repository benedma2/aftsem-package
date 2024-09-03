#' Print method for aftsem xs
#'
#' @description
#' Prints a summary of an aftsem model fit x.
#'
#' @param x An x of class "aftsem", typically the result of a call to `aftsem_fit`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#' The `print.aftsem` method provides a user-friendly summary of the model fit, including
#' the method used for parameter estimation, convergence status, estimated parameters, number of iterations, and
#' the percentage of censored observations.
#'
#' @return
#' The function is called for its side effect, which is printing the summary to the console. It invisibly returns NULL.
#'
#'
#' @seealso
#' \code{\link{aftsem_fit}} for model fitting.
#'
#' @export

print.aftsem <- function(x, ...)
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
  
  cat("Estimated Parameters:\n")
  if (!is.null(x$beta))
  {
    print(x$beta)
  }
  else
  {
    cat("Not available\n")
  }
  cat("\n")
  
if(!is.null(x$intercept))
{
  cat("Estimated Intercept:\n")
  print(x$intercept)
}
  
  
  cat(sprintf("Number of Iterations: %d\n", x$iters))
  cat(sprintf("Number of Observations: %d\n", x$nobs))
  
  
  one_percent = x$nobs/100
  percent_censored <- round(x$censored / one_percent, digits = 2)
  
  cat(sprintf("Percent of Censored Observations: %f\n", percent_censored))
  cat("\n")
  
  cat("Initial Beta Estimate:\n")
  if (!is.null(x$betafirst))
  {
    print(x$betafirst)
  }
  else
  {
    cat("Not available\n")
  }
  cat("\n")
  
  if(x$method == "jin" || x$method == "buckley")
  {
    cat(sprintf("Epsilon (Convergence Tolerance): %.5f\n", x$epsilon))
    cat(sprintf("Maximum Iterations Allowed: %d\n", x$max_iterations))
    cat(sprintf("Algorithm stopped after: %d iterations\n",x$iters))
  }
  cat("\n")
}