#include "utils.h"

#define eps std::numeric_limits<double>::epsilon() // for comparing



/*
 * Computes the buckley_james estimation of intercept
 * However it is recomended not to use it, because Wei in his paper https://doi.org/10.1002/sim.4780111409 states, that it cant be consistently estimated
 */

double buckley_intercept(const arma::vec & y, // survival time
                         const arma::vec & y_star, // delta % y + (1 - delta) % (km + Z*beta_curr);
                         const arma::rowvec & z_mean, // vector of means of all columns from matrix Z
                         const arma::vec & delta, // censoring status
                         const arma::vec & beta) // estimated regression coefficients
{
  double res = 0;
  double help = arma::mean((delta%y + (1-delta)%y_star));
  res = help - (z_mean*beta).eval().at(0,0);
  return res;
}



/*
 * Calculates the fraction with integral of Kaplan-Meier estimation of distribution function of residuals
 * Function is inspired by function lss by Jin https://www.sciencedirect.com/science/article/pii/S0169260706003002
 * Jin function is sadly not available anymore, so I tried to write my own version and optimize it using Armadillo functions
 */

// [[Rcpp::export]]
arma::vec km_e (const arma::vec & e, // residuals
                const arma::vec & delta, // censoring status
                const arma::vec & weights) // either vector of ones or resample weights used in variance estimation
{
 
  // first, we need to rank the observations, because Kaplan-Meier estimator is a step function, that only takes new value when observed event happen 
  arma::uvec ord = arma::sort_index(e);
  arma::vec e_ord = e(ord);
  arma::vec delta_ord = delta(ord);
  
  int nobs = e_ord.n_elem;
  
  // dont forget to order the weights too!
  arma::vec w_ord = weights(ord);
  
  // the KM estimator looks like this
  // S_hat = CUMPROD(1-d/n)
  // so we need to compute subjects at risk at each time, also we need to compute the number of deaths at each time
  // i had decided to compare with std::numeric_limits<double>::epsilon() for better numerical stability
  arma::vec tt = join_vert(arma::conv_to<arma::vec>::from(diff(e_ord) > eps), arma::vec({1}));
  arma::vec tt_1 = join_vert(arma::vec({1}), arma::conv_to<arma::vec>::from(diff(e_ord) > eps));
  arma::vec lin_space = arma::regspace(1, nobs); // (1,2,3,4.....,nobs)
  arma::vec repeats = diff(join_vert(lin_space.elem(find(tt_1 > eps)), arma::vec({static_cast<double>(nobs + 1)})));
  
  arma::vec r = reverse(cumsum(reverse(w_ord)));
  arma::vec d = cumsum(w_ord % delta_ord);
  d = d.elem(find(tt > eps));
  
  // at time=0 no event could have happen (right censoring only)
  d = join_vert(arma::vec({d(0)}), diff(d));
  
  // tt_1 because at time 0 there had to be subjects at risk
  arma::vec s_hat = cumprod(1 - d / r.elem(find(tt_1 > eps))); // KM
  
  int total = arma::accu(repeats); // Total length of the result vector
  arma::vec s_hat_all(total); // final vector of all observations
  int k = 0;
  for (unsigned int i = 0; i < s_hat.n_elem; ++i)
  {
    s_hat_all.subvec(k, k + repeats(i) - 1).fill(s_hat(i));
    k += repeats(i);
  }
  s_hat = s_hat_all;
  arma::vec e_dif = join_vert(diff(e_ord), arma::vec({0}));
  arma::vec e_hat = reverse(cumsum(reverse(e_dif % s_hat)));

  
  e_hat.elem(find(s_hat < eps)).fill(0);
  s_hat.elem(find(s_hat < eps)).fill(1);
  e_hat = e_hat / s_hat + e_ord;
  arma::vec km = e_hat;
  
  // return in original order
  km.elem(ord) = e_hat; 
  return km;
}



