// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]




// [[Rcpp::export]]
/*
 * Hellers estimating function with standard normal distribution function
 */
double compute_heller(const arma::vec &beta, // regression coefficients
                      const arma::vec & y, // survival times
                      const arma::mat & Z, // matrix of covariates
                      const arma::vec & delta, // censoring vecor
                      double a) // bandwith
  
{
  
  int n = y.size();
  arma::vec e = y - Z*beta;
  double loss = 0;
  double diff_e_ij, diff_e_ji, term_ij, term_ji;
  
  for (int i = 0; i < (n-1); ++i)
  {
    for (int j = i+1; j < n; ++j)
    {
      diff_e_ij = (e[j] - e[i])/a;
      diff_e_ji = -diff_e_ij;
      
      term_ij = delta[i] * (e[j] - e[i]) *  arma::normcdf(diff_e_ij) + a *  arma::normpdf(diff_e_ij);
      term_ji = delta[j] * (e[i] - e[j]) *  arma::normcdf(diff_e_ji) + a *  arma::normpdf(diff_e_ji);
      
      loss += term_ij + term_ji;
    }
  }
  
  return loss/n;
}

/******************************************************************************/

// [[Rcpp::export]]
/*
 * Hellers gradient
 */
arma::vec compute_heller_grad(const arma::vec &beta, // regression coefficients
                              const arma::vec & y, // survival times
                              const arma::mat & Z, // covariate matrix
                              const arma::vec & delta, // censoring vector
                              double a) // bandwith
{
  int n = y.size();
  arma::vec e = y - Z*beta;
  arma::vec gr = arma::zeros<arma::vec>(Z.n_cols);
  
  for(int i = 0; i < n; ++i)
  {
    for(int j = 0; j < n; ++j)
    {
      arma::rowvec diff = Z.row(i) - Z.row(j);
      gr += delta[i] * diff.t() * (1-arma::normcdf((e[i] - e[j])/a));
    }
  }
  
  return(gr/n);
}

/******************************************************************************/

/*
 * A matrix in covariance matrix estimation
 */
arma::mat computeA(const arma::vec & delta,// censoring vector
                   const arma::mat & Z, // covariate matrix
                   const arma::vec & e, // vector of residuals
                   const double & a) // bandwith
{
  int n = Z.n_rows;
  int c = Z.n_cols;
  arma::mat A = arma::zeros<arma::mat>(c, c);
  
  for(int i = 0; i < n; ++i)
  {
    for(int j = 0; j < n; ++j)
    {
      arma::rowvec diff = Z.row(i) - Z.row(j);
      A += delta[i] * (diff.t() * diff) * arma::normpdf((e[i] - e[j])/a)/a;
    }
  }
  
  A = A/n;
  return A;
}

/*
 * B matrix in covariance matrix estimation
 */

arma::mat computeB(const arma::vec & delta, // censoring vector
                   const arma::mat & Z, // covariate matrix
                   const arma::vec & e, // vector of residuals
                   const double & a) // bandwith
{
  int n = Z.n_rows;
  int c = Z.n_cols;
  
  double multiplic_first = 0;
  double multiplic_second = 0;
  
  arma::mat B = arma::zeros<arma::mat>(c,c);
  
  for(int i = 0; i < n; ++i)
  {
    for(int j = 0; j < n; ++j)
    {
      for(int k = 0; k < n; ++k)
      {
        if(k == j)
            continue; // Hellers condition
        
        arma::rowvec first = Z.row(i) - Z.row(j);
        arma::rowvec second = Z.row(i) - Z.row(k);
        
        multiplic_first = ((delta[i] * (1 - arma::normcdf((e[i] - e[j])/a))) - (delta[j] * (1- arma::normcdf((e[j] - e[i])/a))));
        multiplic_second = ((delta[i] * (1 - arma::normcdf((e[i] - e[k])/a))) - (delta[k] * (1- arma::normcdf((e[k] - e[i])/a))));
        
        B += (first.t() * second) * multiplic_first * multiplic_second;
      }
    }
  }
  
  return B/std::pow(n,2);
}

/*
 * Estimation of covariance matrix --> A^-1 * B * A^-1
 */

// [[Rcpp::export]]
arma::mat compute_covariance(const arma::vec & delta, // censoring vector
                             const arma::mat & Z, // covariate matrix
                             const arma::vec & e, // vector of residuals
                             const double & a) // bandwith
{
  arma::mat A = computeA(delta,Z,e,a);
  //A.print();
  arma::mat B = computeB(delta,Z,e,a);
  //B.print();
  // pinv is pseudo inversion, returns Moore-Penrose pseudo-inverse
  arma::mat res = arma::pinv(A)*B*arma::pinv(A);
  
  return res;
}
