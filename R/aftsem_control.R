#' Control list for package
#'
#' @param eps Convergence criterion
#' @param maxiter Maximum iterations for algorithms
#' @param gehan_eps Epsilon value for polynomial Gehan optimalization
#' @param optimx.alg Algorithm that will be used in optimx minimalization (see optimx documentation for more details)
#' @param variance.estimation If hellers sd will be estimated
#' @param quantile.method Method used for quantile regression minimalization
#' @param use.grad If excact gradient will be used instead of the numerical one, default is numerical == FALSE
#' 
#' @note When alternating the control list, one must write other variables also. Example: When user want to estimate the Hellers covariance matrix he would need
#' to change the control list -> aftsem(....., control = list(variance.estimation = TRUE, use.grad = FALSE, optimx.alg = "BFGS))
#'
#' @return list of parameters above
#'
aftsem.control <- function(eps = 10^-5,
                           maxiter = 15,
                           gehan_eps = 10^-6,
                           optimx.alg = "BFGS",
                           variance.estimation = FALSE,
                           quantile.method = "br",
                           use.grad = FALSE)
{
  list(eps = eps, maxiter = maxiter, gehan_eps = gehan_eps, optimx.alg = optimx.alg, variance.estimation = variance.estimation, quantile.method = quantile.method, use.grad = use.grad)
}