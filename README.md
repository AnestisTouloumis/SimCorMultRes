<!-- README.md is generated from README.Rmd. Please edit that file -->
SimCorMultRes: Simulates Correlated Multinomial Responses
=========================================================

[![Travis-CI Build Status](https://travis-ci.org/AnestisTouloumis/SimCorMultRes.svg?branch=master)](https://travis-ci.org/AnestisTouloumis/SimCorMultRes) [![develVersion](https://img.shields.io/badge/devel%20version-1.5.0.9-brightgreen.svg?style=flat)](https://github.com/AnestisTouloumis/SimCorMultRes) [![Project Status: Active The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Last-changedate](https://img.shields.io/badge/last%20change-2017--07--12-brightgreen.svg)](/commits/master)

[![CRAN Version](http://www.r-pkg.org/badges/version/SimCorMultRes?color=blue)](https://cran.r-project.org/package=SimCorMultRes) [![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/SimCorMultRes?color=blue)](http://cranlogs.r-pkg.org/badges/grand-total/SimCorMultRes) [![CRAN Downloads](http://cranlogs.r-pkg.org/badges/SimCorMultRes)](http://cran.rstudio.com/web/packages/SimCorMultRes/index.html)

Author
------

Anestis Touloumis: <https://brighton.academia.edu/AnestisTouloumis>

School of Computing, Engineering and Mathematics, University of Brighton.

Installation
------------

You can install the release version of **SimCorMultRes**:

``` r
install.packages("SimCorMultRes")
```

The source code for the release version of **SimCorMultRes** is available on CRAN at:

-   <https://CRAN.R-project.org/package=SimCorMultRes>

Or you can install the development version of **SimCorMultRes**:

``` r
# install.packages('devtools')
devtools::install_github("AnestisTouloumis/SimCorMultRes")
```

The source code for the development version of **SimCorMultRes** is available on github at:

-   <https://github.com/AnestisTouloumis/SimCorMultRes>

To use **SimCorMultRes**, you should load the package as follows:

``` r
library(SimCorMultRes)
```

Usage and functions
-------------------

This package provides functions to simulate correlated binary, ordinal and nominal responses, which are drawn as realizations of a latent regression model for continuous random vectors as proposed by Touloumis (2016).

There are four core functions:

-   `rbin` to simulate correlated binary responses,
-   `rmult.bcl` to simulate correlated nominal multinomial responses,
-   `rmult.clm` to simulate correlated ordinal responses under a marginal cumulative link model,
-   `rmult.clm` to simulate correlated ordinal responses under a marginal continuation-ratio link model.

There are also two utility functions:

-   `rnorta` for simulating continuous or discrete random vectors with prescribed marginal distributions using the NORTA method,
-   `rsmvnorm` for simulating continuous random vectors from a multivariate normal distribution.

Example
-------

The following R code illustrates how to use the core function `rbin`:

``` r
## See Example 3.4 in the Vignette.
set.seed(123)
N <- 5000
clsize <- 4
intercepts <- 0
betas <- 0.2
cor.matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
x <- rep(rnorm(N), each = clsize)
CorBinRes <- rbin(clsize = clsize, intercepts = intercepts, betas = betas, xformula = ~x, 
    cor.matrix = cor.matrix, link = "probit")
library(gee)
binGEEmod <- gee(y ~ x, family = binomial("probit"), id = id, data = CorBinRes$simdata)
#> Beginning Cgee S-function, @(#) geeformula.q 4.13 98/01/27
#> running glm to get initial regression estimate
#> (Intercept)           x 
#> 0.002636705 0.204827031
summary(binGEEmod)$coefficients
#>                Estimate  Naive S.E.    Naive z Robust S.E.   Robust z
#> (Intercept) 0.002636705 0.008929290  0.2952872  0.01572132  0.1677153
#> x           0.204827031 0.009114596 22.4724192  0.01610695 12.7166857
```

Additional examples can be found in Touloumis (2016) and in the vignette of **SimCorMultRes**.

``` r
browseVignettes("SimCorMultRes")
```

How to cite
-----------


    To cite R package SimCorMultRes in publications, please use:

      Touloumis, A. (2016). Simulating Correlated Binary and
      Multinomial Responses under Marginal Model Specification: The
      SimCorMultRes Package. The R Journal 8:2, 79-91.

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {Simulating Correlated Binary and Multinomial Responses under 
             Marginal Model Specification: The SimCorMultRes Package},
        author = {Anestis Touloumis},
        year = {2016},
        journal = {The R Journal},
        volume = {8},
        number = {2},
        pages = {79-91},
        url = {https://journal.r-project.org/archive/2016/RJ-2016-034/index.html},
      }

References
==========

Touloumis, Anestis. 2016. “Simulating Correlated Binary and Multinomial Responses Under Marginal Model Specification: The SimCorMultRes Package.” *The R Journal* 8 (2): 79–91. <https://journal.r-project.org/archive/2016/RJ-2016-034/index.html>.
