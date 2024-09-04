#' Accelerated Failure Time Semiparametric Model
#'
#' @param formula A formula expression, of the form \code{response ~ predictors}. Response must be a Surv object
#' @param data An optional data.frame in which to interpret the variables in the \code{formula}.
#' @param control Control parameters for the AFT model.
#' @param method A character string specifying the method to be used (buckley,jin,gehan,gehan-heller,gehan-poly).
#' @param binit Initial values for the regression parameters.
#' @param ties A method to handle ties in the failure times. If ties = NULL only warning will be printed. If ties = jitter, the data will be augumented
#' @param na.action A method to deal with missing values (na.fail)
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#' @param resample Number of resamples for variance estimation for gehan and jin methods.
#' @param ... Additional arguments.
#' @return A list representing the fit
#'  - `call`: Call of the function
#'  - `cnames`: Column names
#'  - `method`: Method of estimation
#'  - `nobs`: Number of observations
#'  - `censored`: Number of censored observations
#'  - `betafirst`: Initial beta
#'  - `epsilon`: Epsilon in convergence criterion
#'  - `max_iterations`: Max iterations for buckley and jin method
#'  - `resample`: Resample number
#'  - `objects from aftsem.fit`: All the object from fit function
#' @export
#'
#' @examples
#' # Generating example data
#' library(survival)
#' set.seed(123) # for reproducibility
#' n <- 100 # number of observations
#' Z <- matrix(rnorm(n*2), ncol = 2) # two covariates
#' beta <- c(0.5, -0.25) # true coefficients
#' times <- exp(Z %*% beta + rnorm(n)) # simulated survival times
#' censoring <- runif(n,0,30)
#' observed_times <- times
#' delta <- 1 * (times<=censoring)
#'
#' # Fit the model
#' 
#' 
#' fit <- aftsem(Surv(log(observed_times), delta) ~ Z[,1] + Z[,2],
#'               method = "buckley",
#'               binit = "auto",
#'               ties = "NULL",
#'               na.action = na.omit,
#'               subset = NULL
#' )
#' 
#' # Print the summary
#' summary(fit)
#' 
#'
aftsem <- function(formula,
                   data,
                   control = aftsem.control(),
                   method = "buckley",
                   binit = "auto",
                   ties = NULL,
                   na.action = na.omit,
                   subset = NULL,
                   resample = 0,
                   ...)
{
  
  ##############################################################################
  
  
  Call <- match.call() # save call for return list
  mf <- match.call(expand.dots = FALSE) # working with this call, expand.dots includes parameters in ...
  mn <- match(c("formula", "data", "subset", "na.action"), names(mf), 0) # vector of positions
  mf <- mf[c(1, mn)] # take only columns from mn
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, sys.parent())
  Terms <- attr(mf, "terms")
  zv <- as.character(attr(Terms, "variables")) # variables
  yv <- attr(Terms, "response") # response variable
  if (yv > 0) 
  {
    zv <- zv[-yv] # remove response variable from xv
  }
  
  y <- model.extract(mf, "response")
  
  if (!inherits(y, "Surv")) 
  {
    stop("The response variable must be of type 'Surv' with right censoring.")
  }
  
  if (is.null(ties)) 
  {
    time_var <- y[, 1]
    if (any(duplicated(time_var))) 
    {
      warning("There are ties in the failure times.")
    }
  }
  
  if (is.null(ties) == FALSE)
  {
    if(ties == "jitter")
    {
      perturbed_times <- y[, 1] + seq_along(y[, 1]) * .Machine$double.eps^0.5
      y[, 1] <- perturbed_times
    }
  }
  
  
  ##############################################################################
  
  Z <- model.matrix(Terms, mf)
  intercept <- FALSE
  
  if (all(Z[, 1] == 1)) 
  {
    Z <- Z[, -1] # remove intercept
    intercept <- TRUE
  }
  
  Z <- as.matrix(Z)
  
  #if (any(y[, 1] < 0)) 
  #{
  #  stop("Survival times must be positive!")
  #}
  
  nobs <- nrow(Z)
  nvar <- ncol(Z)
  betafirst <- NULL
  
  if (is.numeric(binit)) 
  {
    if (length(binit) != nvar) 
    {
      stop("Initial beta values have wrong dimension")
    }
    betafirst <- binit
  } 
  else 
  {
    if (!(binit %in% c("gehan", "lm","auto"))) 
    {
      stop("Aftsem supports only least-squares or gehan estimation as initial value.")
    } 
    else 
    {
      if(binit != "auto")
      {
        warning("We recommend to put binit value as 'auto'.")
      }
      
      if(binit == "auto" && method == "buckley") binit <- "lm"
      if(binit == "auto" && method == "jin") binit <- "gehan"
      if(binit == "auto" && method == "gehan-poly") binit <- "lm"
      if(binit == "lm" && method == "gehan") binit <- "gehan"
      if(binit == "auto" && method == "gehan-heller") binit <- "lm"
      
      
      fit <- NULL
      if(binit == "lm")
      {
        ifelse(intercept, fit <- lm(y[, 1] ~ Z), fit <- lm(y[, 1] ~ Z - 1)) 
        betafirst <- fit$coef
        if (intercept) 
        {
          betafirst <- betafirst[-1]
        }
      }
      if (binit == "gehan")
      {
        if (method == "gehan") betafirst <- rep(0,nvar)
        else
        {
          rs <- matrix(rep(1,nobs), ncol=1)
          gehanres <- gehan_estimation(y[, 1], Z, y[, 2], rs,control$quantile.method, init = TRUE)
          betafirst <- gehanres
        }
      }
    }
  }
  
  if (resample < 0) resample <- 0
  if (resample > 0 && resample < 2) resample <- 3
  
  stime <- y[,1]
  delta <- y[,2]
  
  # rank the observations
  if (method %in% c("gehan","gehan-heller","gehan-poly"))
  {
   ord <- order(stime)
   stime <- stime[ord]
   delta <- delta[ord]
   Z <- Z[ord,]
   Z <- as.matrix(Z) # if nvar is 1, the ranking of matrix Z changes its class to int. 
  }
  
  
  fit <- aftsem_fit(Z, stime, delta, betafirst, method, control, intercept, resample, nobs)
  
  fit$call <- Call
  fit$cnames <- dimnames(Z)[[2]]
  fit$method <- method
  fit$nobs <- nobs
  fit$censored <- nobs - sum(delta)
  fit$betafirst <- betafirst
  fit$epsilon <- control$eps
  fit$max_iterations <- control$maxiter
  fit$resample <- resample
  
  return(fit)
}
