// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


/*
 * Definition of c_epsilon function
 */

double c_epsilon(double x, // difference of residuals
                 double epsilon) // chosen epsilon for comparing
{
  if (x < -epsilon)
  {
    return -x;
  } 
  else if (x <= epsilon && x >= -epsilon)
  {
    return -1 / (16 * std::pow(epsilon, 3)) * std::pow((x + epsilon), 4) + 
      1 / (4 * std::pow(epsilon, 2)) * std::pow((x + epsilon), 3);
  }
  else
  {
    return 0;
  }
}


/*
 * Polynomial smoothed Gehan estimator
 */

// [[Rcpp::export]]
double compute_f_epsilon(const arma::vec &beta, // regression coefficients
                         const arma::vec & y, // survival times
                         const arma::mat & Z, // matrix of covariates
                         const arma::vec & delta, // censoring vector
                         double epsilon) // chosen epsilon for comparing

{
  int n = y.size();
  arma::vec e = y - Z*beta;
  double loss = 0;
  
  
  
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      loss += delta[i] * c_epsilon(e[i]-e[j],epsilon);
    }
  }
  
  return loss/n;
}

/******************************************************************************/

/*
 * definition of k_epsilon
 */
double k_epsilon(double x, // difference of reziduals
                 double epsilon) // chosen epsilon for comparing
{
  
    if (x <= -epsilon)
    {
      return 1;
    } 
    if ((x > -epsilon) && (x <= epsilon))
    {
      return (-1 /(4 * std::pow(epsilon, 3)) * std::pow((epsilon - x), 3) + 
        3/(4 * std::pow(epsilon, 2)) * std::pow((epsilon - x), 2));
    }
    return 0;
}

/*
 * Implementation of gradient of k_epsilon function
 */

// [[Rcpp::export]]
arma::vec compute_f_epsilon_grad(const arma::vec &beta, // regression coefficients
                                 const arma::vec & y, // survival times
                                 const arma::mat & Z, // covariate matrix
                                 const arma::vec & delta, // censoring vector
                                 double epsilon) // chosen epsilon for comparing
{
  int n = y.size();
  arma::vec e = y - Z*beta;
  arma::vec gr = arma::zeros<arma::vec>(Z.n_cols);
  
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      arma::rowvec diff = Z.row(i) - Z.row(j);
      double ke = k_epsilon(e[i] - e[j], epsilon);
      gr += delta[i] * diff.t() * ke;
    }
  }
  
  return gr/n;
}