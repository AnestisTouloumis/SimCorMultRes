
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SimCorMultRes: Simulates Correlated Multinomial Responses

[![Github
version](https://img.shields.io/badge/GitHub%20-1.9.0-orange.svg)](%22commits/master%22)
[![R-CMD-check](https://github.com/AnestisTouloumis/SimCorMultRes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AnestisTouloumis/SimCorMultRes/actions/workflows/R-CMD-check.yaml)
[![Project Status: Active The project has reached a stable, usable state
and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/AnestisTouloumis/SimCorMultRes/branch/master/graph/badge.svg)](https://codecov.io/gh/AnestisTouloumis/SimCorMultRes?branch=master)

[![CRAN
Version](https://www.r-pkg.org/badges/version/SimCorMultRes?color=blue)](https://cran.r-project.org/package=SimCorMultRes)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/SimCorMultRes?color=blue)](https://cranlogs.r-pkg.org/badges/grand-total/SimCorMultRes)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/SimCorMultRes)](https://cran.r-project.org/package=SimCorMultRes)

## Installation

You can install the release version of `SimCorMultRes`:

``` r
install.packages("SimCorMultRes")
```

The source code for the release version of `SimCorMultRes` is available
on CRAN at:

- <https://CRAN.R-project.org/package=SimCorMultRes>

Or you can install the development version of `SimCorMultRes`:

``` r
# install.packages('devtools')
devtools::install_github("AnestisTouloumis/SimCorMultRes")
```

The source code for the development version of `SimCorMultRes` is
available on github at:

- <https://github.com/AnestisTouloumis/SimCorMultRes>

To use `SimCorMultRes`, you should load the package as follows:

``` r
library("SimCorMultRes")
```

## Usage and functions

This package provides five core functions to simulate correlated binary
(`rbin`), nominal (`rmult.bcl`) and ordinal (`rmult.acl`, `rmult.clm`
and `rmult.crm`) responses, which are drawn as realizations of a latent
regression model for continuous random vectors as proposed by Touloumis
(2016):

- `rbin` to simulate correlated binary responses under a marginal model
  with logit, probit, cloglog and cauchit link function,
- `rmult.bcl` to simulate correlated nominal multinomial responses under
  a marginal baseline-category logit model,
- `rmult.acl` to simulate correlated ordinal responses under a marginal
  adjacent-category logit model,
- `rmult.clm` to simulate correlated ordinal responses under a marginal
  cumulative link model,
- `rmult.crm` to simulate correlated ordinal responses under a marginal
  continuation-ratio link model.

All five functions, assume that you provide either the correlation
matrix of the multivariate normal distribution in NORTA (via
`cor.matrix`) or the values of the latent responses (via the `rlatent`).
Based on a simulation study (see Section 3.5 of the vignette and dataset
`simulation`), it is indicated that the correlation matrix of the
multivariate normal distribution used in the NORTA method (via
`cor.matrix`) can be considered a reliable approximation of the actual
correlation matrix of the latent responses generated by the NORTA
method. This appears to be the case irrespective of the marginal
distributions of the latent responses for all the threshold approaches
implemented in `SimCorMultRes`.

There are also two utility functions:

- `rnorta` for simulating continuous or discrete random vectors with
  prescribed marginal distributions using the NORTA method,
- `rsmvnorm` for simulating continuous random vectors from a
  multivariate normal distribution.

## Example

The following R code illustrates how to use the core function `rbin`:

``` r
## See Example 3.5 in the Vignette.
set.seed(123)
## define number of random vectors
sample_size <- 100
## define size of each random vector
cluster_size <- 4
## define intercept of the binary probit regression model
beta_intercepts <- 0
## define coefficients of the explanatory variables
beta_coefficients <- 0.2
## provide explanatory variables
x <- rep(rnorm(sample_size), each = cluster_size)
## define correlation matrix for the multivariate normal distribution in NORTA
latent_correlation_matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
## use rbin function to create the desired dataset
simulated_binary_responses <- rbin(clsize = cluster_size, intercepts = beta_intercepts,
    betas = beta_coefficients, xformula = ~x, cor.matrix = latent_correlation_matrix,
    link = "probit")
library("gee")
binary_gee_model <- gee(y ~ x, family = binomial("probit"), id = id, data = simulated_binary_responses$simdata)
#> Beginning Cgee S-function, @(#) geeformula.q 4.13 98/01/27
#> running glm to get initial regression estimate
#> (Intercept)           x 
#>   0.1315121   0.2826005
summary(binary_gee_model)$coefficients
#>              Estimate Naive S.E.  Naive z Robust S.E. Robust z
#> (Intercept) 0.1315121 0.06399465 2.055048   0.1106696 1.188331
#> x           0.2826006 0.07191931 3.929412   0.1270285 2.224703
```

Additional examples can be found in Touloumis (2016) and in the vignette
of `SimCorMultRes`. To access these two documents, run the following
command:

``` r
browseVignettes("SimCorMultRes")
```

## How to cite

    To cite 'SimCorMultRes' in publications, please use:

      Touloumis A (2016). "Simulating Correlated Binary and Multinomial
      Responses under Marginal Model Specification: The SimCorMultRes
      Package." _The R Journal_, *8*(2), 79-91. R package version 1.9.0,
      <https://journal.r-project.org/archive/2016/RJ-2016-034/index.html>.

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {Simulating Correlated Binary and Multinomial Responses under
             Marginal Model Specification: The SimCorMultRes Package},
        author = {Anestis Touloumis},
        year = {2016},
        journal = {The R Journal},
        volume = {8},
        number = {2},
        note = {R package version 1.9.0},
        pages = {79-91},
        url = {https://journal.r-project.org/archive/2016/RJ-2016-034/index.html},
      }

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Touloumis2016" class="csl-entry">

Touloumis, A. (2016) [Simulating Correlated Binary and Multinomial
Responses under Marginal Model Specification: The SimCorMultRes
Package](https://journal.r-project.org/archive/2016/RJ-2016-034/index.html).
*The R Journal*, **8**, 79–91.

</div>

</div>
