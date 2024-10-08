# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

compute_f_epsilon <- function(beta, y, Z, delta, epsilon) {
    .Call(`_aftsem_compute_f_epsilon`, beta, y, Z, delta, epsilon)
}

compute_f_epsilon_grad <- function(beta, y, Z, delta, epsilon) {
    .Call(`_aftsem_compute_f_epsilon_grad`, beta, y, Z, delta, epsilon)
}

compute_heller <- function(beta, y, Z, delta, a) {
    .Call(`_aftsem_compute_heller`, beta, y, Z, delta, a)
}

compute_heller_grad <- function(beta, y, Z, delta, a) {
    .Call(`_aftsem_compute_heller_grad`, beta, y, Z, delta, a)
}

compute_covariance <- function(delta, Z, e, a) {
    .Call(`_aftsem_compute_covariance`, delta, Z, e, a)
}

estimate_buckley <- function(y, Z, delta, beta0, epsilon, max_iters) {
    .Call(`_aftsem_estimate_buckley`, y, Z, delta, beta0, epsilon, max_iters)
}

estimate_jin <- function(y, Z, delta, beta0, epsilon, max_iters, Resample_mat, Beta_star, sampling_used) {
    .Call(`_aftsem_estimate_jin`, y, Z, delta, beta0, epsilon, max_iters, Resample_mat, Beta_star, sampling_used)
}

km_e <- function(e, delta, weights) {
    .Call(`_aftsem_km_e`, e, delta, weights)
}

