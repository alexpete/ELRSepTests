
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ELRSepTests

<!-- badges: start -->

<!-- badges: end -->

`ELRSepTests` provides tools for empirical likelihood ratio tests for
assessing the degree of separability in the covariance structure of
two-way data. Examples include matrix-, multivariate functional-, or
surface-valued data

The package implements methods based on (constrained) empirical
likelihood evaluation and optimization, combining R and C++ (via Rcpp)
for efficient computation.

## Installation

You can install the development version of ELRSepTests from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("alexpete/ELRSepTests")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ELRSepTests)

set.seed(123)

# Small matrix-valued sample: n observations of a 3 x 3 matrix
n <- 20
p1 <- 3
p2 <- 3

# Everything is independent, so covariance is separable (and therefore 
# also weakly and partially separable)
X <- array(rnorm(n * p1 * p2), dim = c(n, p1, p2))

# Run the main testing function
res <- ELRSepTests(
  X = X,
  JTest = 2,
  LTest = 2,
  nullHyp = c('ParSep', 'WkSep', 'Sep'),
  B = 100
)


res$tStats
#>     ParSep      WkSep        Sep 
#> 4.56986975 6.87955536 0.05054832
res$bootRes$bootPval
#>    ParSep     WkSep       Sep 
#> 0.6039604 0.3762376 0.9207921
```
