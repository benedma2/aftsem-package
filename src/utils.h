#ifndef TILS_H
#define TILS_H

#include "RcppArmadillo.h"

arma::vec km_e (const arma::vec & e, const arma::vec & delta, const arma::vec & weights);
double buckley_intercept(const arma::vec & y, const arma::vec & y_star, const arma::rowvec & z_mean, const arma::vec & delta, const arma::vec & beta);

#endif