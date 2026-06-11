
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

## Dependencies

The package requires the packages `MASS`, `Rcpp`, `RcppArmadillo`,
`melt`, and `covsep` to be installed. All but `covsep` are available on
CRAN for R version 4.6.0, and `covsep` can be installed using the [CRAN
archive](https://cran.r-project.org/web/packages/covsep/index.html).

## Installation

You can install the development version of ELRSepTests from
[GitHub](https://github.com/) with:

``` r
# Option 1: using remotes (widely supported)
install.packages("remotes")
remotes::install_github("alexpete/ELRSepTests")

# Option 2: using pak (faster, recommended)
install.packages("pak")
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
  JTest = 2L,
  LTest = 2L,
  nullHyp = c('ParSep', 'WkSep', 'Sep'),
  B = 100
)


res$tStats
#>            (2, 2)
#> ParSep 4.56986975
#> WkSep  6.87955536
#> Sep    0.05054868
res$bootPval
#>           (2, 2)
#> ParSep 0.7029703
#> WkSep  0.3861386
#> Sep    0.9306931
```
