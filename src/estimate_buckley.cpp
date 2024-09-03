// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
Rcpp::List estimate_buckley(const arma::vec &  y, // survival times
                            const arma::mat &  Z, // covariate matrix
                            const arma::vec &  delta, // censoring status
                            const arma::vec &  beta0, // initial beta estimation
                            const double epsilon, // convergence criterion epsilon
                            const int max_iters // maximum iterations 
                            )
{
  // 1) asign starting values
  // 2) compute the expectation
  // 3) compute \beta^{i+1}
  // 4) check if converges, otherwise repeat step 2
  
  int nobs = Z.n_rows;
  bool converged = false;
  
  // precompute matrixes and vectors for later use
  arma::rowvec z_sum = arma::sum(Z);
  arma::vec beta_curr = beta0;
  arma::vec beta_prev = beta_curr;
  
  arma::rowvec z_mean = z_sum / static_cast<double>(nobs);
  arma::mat Zmean = arma::repmat(z_mean, nobs, 1);
  
  arma::mat Zstar = Z - Zmean;
  arma::mat Zsolve = arma::pinv(Zstar.t()*Zstar);
  arma::vec weights = arma::ones<arma::vec>(nobs); // for possible variance estimation
  
  int iter = 0; // iteration tracker
  double alpha = 0; //intercept
  
  arma::vec e (nobs);
  arma::vec km (nobs);
  arma::vec y_star (nobs);
  
  while(iter < max_iters)
  {
    iter++;
    beta_prev = beta_curr;
    
    e = y - Z*beta_curr; // residuals
    km = km_e(e,delta,weights); // integral
    y_star = delta % y + (1 - delta) % (km + Z*beta_curr);
    beta_curr = Zsolve * Zstar.t()*y_star;  
    
    if (arma::norm(beta_curr - beta_prev) < epsilon * (arma::norm(beta_curr) + 0.000001)) // convergence criterion was met
    {
      alpha = buckley_intercept(y,y_star,z_mean,delta,beta_curr); // cannot be realible estimated, it is not advised to use
      converged = true;
      break;
    }
  }
  
  Rcpp::List out;
  if(!converged)
    out["BETA"] = beta_curr;
  else
    out["BETA"] = beta_curr;
  
  out["THAT"] = y_star;
  out["RESID"] = y_star - Z*beta_curr;
  out["CONVERGED"] = converged;
  out["ITERATIONS"] = iter;
  out["INTERCEPT"] = alpha;
  return out;
}