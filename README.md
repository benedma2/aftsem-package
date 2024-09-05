# `aftsem` 📦

The `aftsem` package implements several algorithms for estimating coefficients in semiparametric accelerated failure time (AFT) models. These models are widely used in survival analysis to handle censored data, where the goal is to estimate regression coefficients that relate covariates to failure times.

## Status

The package is under active development.

## Key Features ✨

The package offers two main approaches for estimating the coefficients:

1. **Least Squares Methods** 🧮
   - Based on the classical Buckley-James algorithm. 
   - Two basic algorithms are implemented here:
     - [BUCKLEY, Jonathan; JAMES, Ian. Linear Regression with Censored Data. Biometrika. 
     1979, roč. 66, č. 3, s. 429–436 . issn 00063444. Available on: http:
	//www.jstor.org/stable/2335161.] 📖
     - [JIN, Zhezhen; LIN, DY; YING, Zhiliang. On least-squares regression with censored data.
Biometrika. 2006, roč. 93, č. 1, s. 147–161.] 📖

2. **Rank-Based Methods** 📊
   - Inspired by the methods of Jin and Tsiatis. 
   - Three rank-based algorithms are implemented:
     - [ JIN, Zhezhen; LIN, D. Y.; WEI, L. J.; YING, Zhiliang. Rank-Based Inference for the Ac-
celerated Failure Time Model. Biometrika . 2003, roč. 90, č. 2, s. 341–353 [cit. 2024-
02-13]. issn 00063444. Available on: http://www.jstor.org/stable/30042044.] 📖
     - [HELLER, Glenn. Smoothed rank regression with censored data. Journal of the American
Statistical Association. 2007, roč. 102, č. 478, s. 552–559.] 📖
     - [CHUNG, Matthias; LONG, Qi; JOHNSON, Brent A. A tutorial on rank-based coefficient
estimation for censored data in small-and large-scale problems. Statistics and computing.
2013, roč. 23, s. 601–614.] 📖

## Installation 🛠️

You can install the development version of the `aftsem` package from CRAN:

``` r
> install.packages("aftsem")
> library(aftsem)
```
