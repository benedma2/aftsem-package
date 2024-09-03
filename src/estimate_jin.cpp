// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
Rcpp::List estimate_jin(const arma::vec &  y, // survival times
                        const arma::mat &  Z, // covariate matrix
                        const arma::vec &  delta, // censoring vector
                        const arma::vec &  beta0, // initial beta estimates
                        const double epsilon, // convergence criterion epsilon
                        const int max_iters, // maximum iterations allowed
                        const arma::mat & Resample_mat, // resample mat, containing random variables R
                        arma::mat & Beta_star, // initial beta estimates for resampling from gehan estimation
                        bool sampling_used) // bool indicator if resampling will be used
{
  
  int nobs = Z.n_rows;
  int nvar = Z.n_cols;
  
  //
  int resample_size = Resample_mat.n_cols;
  //
  
  bool converged = false;
  
  // precompute later used matrixes and vectors, this saves a lot of computer memory  and its more efficient
  arma::rowvec z_sum = arma::sum(Z);
  arma::vec beta_curr = beta0;
  arma::vec beta_prev = beta_curr;
  
  arma::rowvec z_mean = z_sum / static_cast<double>(nobs);
  arma::mat Zmean = arma::repmat(z_mean, nobs, 1);
  
  arma::mat Zstar = Z - Zmean;
  
  // get solving matrix
  // notice that is used ::pinv, this is Moore-Penrose pseudo-inverse, computation is done using SVD decomposition, see armadillo documentation for more info https://arma.sourceforge.net/docs.html#pinv
  arma::mat Zsolve = arma::pinv(Zstar.t()*Zstar);
  
  arma::vec dummy_res = arma::ones<arma::vec>(nobs); // for normal beta, for Beta_star we will use resample mat
  
  /****************************************************************************/
  // precomputes solve matrixes for resampling, again this is more efficient
  arma::cube Zsolve_star(nvar, nvar, resample_size);
  if(sampling_used)
  {
    for(int i = 0; i < resample_size; ++i)
    {
      Zsolve_star.slice(i) = arma::inv(Zstar.t() * (Zstar.each_col() % Resample_mat.col(i))); //Zstar.each_col() % Resample_mat.col(i)
    }
  }
  /****************************************************************************/
  
  int iter = 0; // iteration tracker
  
  arma::vec e (nobs);
  arma::vec km (nobs);
  arma::vec y_star (nobs);
  
  /****************************************************************************/
  arma::mat tmp (nobs,nvar);
  arma::vec tmpy (nobs);
  /****************************************************************************/
  
  double y_mean = 0;
  
  while(iter < max_iters)
  {
    iter++;
    beta_prev = beta_curr;
    
    e = y - Z*beta_curr; // residuals
    km = km_e(e,delta,dummy_res); // integral
    y_star = delta % y + (1 - delta) % (km + Z*beta_curr);
    y_mean = arma::mean(y_star);
    beta_curr = Zsolve * Zstar.t()*(y_star-y_mean);
    
    /****************************************************************************/
    // Resampling
    if(sampling_used)
    {
      // basically we perform the excact same operation as before, but now with perturbed fumction, also we need to repeat this process for all columns in beta star
      for(int i = 0; i < resample_size; ++i)
      {
        e = y - Z*Beta_star.col(i);
        km = km_e(e,delta,Resample_mat.col(i));
        y_star = delta % y + (1-delta) % (km + Z*Beta_star.col(i));
        y_mean = arma::mean(y_star);
        tmp = (Zstar.each_col() % Resample_mat.col(i));
        tmpy = (y_star-y_mean);
        Beta_star.col(i) = Zsolve_star.slice(i) * (tmp.t() * tmpy);
      }
    }  
    /****************************************************************************/
    
    if (arma::norm(beta_curr - beta_prev) < epsilon * (arma::norm(beta_curr) + 0.000001)) // convergence criterion was met
    {
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
  
  /****************************************************************************/
  if(sampling_used)
    out["DISTRIBUTION"] = Beta_star;
  /****************************************************************************/
  
  return out;
}
