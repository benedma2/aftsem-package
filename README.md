# `aftsem` 📦

The aftsem package provides multiple algorithms for estimating coefficients in semiparametric accelerated failure time (AFT) models, commonly used in survival analysis to address censored data. These models estimate regression coefficients that relate covariates to failure times. The package supports modeling right-censored data and enables comparisons of various estimation algorithms in terms of precision and speed. Additionally, it is designed for easy extension, allowing users to implement their own methods as needed.

## Status

The package is under active development.

## Key Features ✨

The package offers two main approaches for estimating the coefficients:

1. **Least Squares Methods** 🧮
   - Based on the classical Buckley-James algorithm. 
   - Two basic algorithms are implemented here:
     - [BUCKLEY, Jonathan; JAMES, Ian. Linear Regression with Censored Data. Biometrika. 
     1979, doi:10.2307/2335161 Available on: http:
	//www.jstor.org/stable/2335161.] 📖
     - [JIN, Zhezhen; LIN, DY; YING, Zhiliang. On least-squares regression with censored data. doi:10.1093/biomet/93.1.147] 📖

2. **Rank-Based Methods** 📊
   - Inspired by the methods of Jin and Tsiatis. 
   - Three rank-based algorithms are implemented:
     - [ JIN, Zhezhen; LIN, D. Y.; WEI, L. J.; YING, Zhiliang. Rank-Based Inference for the Ac-
celerated Failure Time Model. doi:10.1093/biomet/90.2.341 Available on: http://www.jstor.org/stable/30042044.] 📖
     - [HELLER, Glenn. Smoothed rank regression with censored data. Journal of the American
Statistical Association. doi:10.1198/016214506000001257 ] 📖
     - [CHUNG, Matthias; LONG, Qi; JOHNSON, Brent A. A tutorial on rank-based coefficient
estimation for censored data in small-and large-scale problems. Statistics and computing. doi:10.1007/s11222-012-9333-9] 📖

## Installation 🛠️

You can install the `aftsem` package from CRAN:

``` r
> install.packages("aftsem")
> library(aftsem)
```

