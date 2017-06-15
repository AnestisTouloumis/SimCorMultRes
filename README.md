<!-- README.md is generated from README.Rmd. Please edit that file -->
SimCorMultRes: GEE Solver for Correlated Nominal or Ordinal Multinomial Responses
=================================================================================

[![Travis-CI Build Status](https://travis-ci.org/AnestisTouloumis/SimCorMultRes.svg?branch=master)](https://travis-ci.org/AnestisTouloumis/SimCorMultRes) [![develVersion](https://img.shields.io/badge/devel%20version-1.4.3-brightgreen.svg?style=flat)](https://github.com/AnestisTouloumis/SimCorMultRes) [![Project Status: Active The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Last-changedate](https://img.shields.io/badge/last%20change-2017--06--15-brightgreen.svg)](/commits/master)

[![CRAN Version](http://www.r-pkg.org/badges/version/SimCorMultRes?color=blue)](https://cran.r-project.org/package=SimCorMultRes) [![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/SimCorMultRes?color=blue)](http://cranlogs.r-pkg.org/badges/grand-total/SimCorMultRes) [![CRAN Downloads](http://cranlogs.r-pkg.org/badges/SimCorMultRes)](http://cran.rstudio.com/web/packages/SimCorMultRes/index.html)

Author
------

Anestis Touloumis <https://brighton.academia.edu/AnestisTouloumis>

School of Computing, Engineering and Mathematics, University of Brighton.

Installation
------------

You can install the release version of the **SimCorMultRes** package from CRAN:

``` r
install.packages("SimCorMultRes")
```

or the development version from github:

``` r
install.packages("devtools")  # if you have not installed 'devtools' package
devtools::install_github("AnestisTouloumis/SimCorMultRes")
```

The source code for the **SimCorMultRes** package is available on github at

-   <https://github.com/AnestisTouloumis/SimCorMultRes>.

To use **SimCorMultRes**, you should load the package as follows:

``` r
library(SimCorMultRes)
```

Usage and functions
-------------------

The **SimCorMultRes** package provides functions to simulate correlated binary or multinomial responses. These responses are drawn as realizations of a latent regression model for continuous random vectors with the correlation structure expressed in terms of the latent correlation (Touloumis 2016).

The package contains four core functions:

-   `rbin` to simulate correlated binary responses,
-   `rmult.bcl` to simulate correlated nominal multinomial responses,
-   `rmult.clm` to simulate correlated ordinal responses under a marginal cumulative link model,
-   `rmult.clm` to simulate correlated ordinal responses under a marginal continuation-ratio link model.

The **SimCorMultRes** package also offers two utility functions:

-   `rnorta` for simulating random vectors using the NORTA method,
-   `rsmvnorm` for simulating continuous random vectors from a multivariate normal distribution.

Example
-------

The following R code illustrates how to use the core functions of **SimCorMultRes**.

``` r
x <- 5
```

Getting help
------------

A detailed description of **SimCorMultRes** can be found in Touloumis (2016). The vignette of **SimCorMultRes** contains additional examples and it be accessed by executing the following R commands:

``` r
browseVignettes("SimCorMultRes")
```

How to cite
-----------

``` r
citation("SimCorMultRes")
#> 
#> To cite R package SimCorMultRes in publications, please use:
#> 
#>   Touloumis, A. (2016). Simulating Correlated Binary and
#>   Multinomial Responses under Marginal Model Specification: The
#>   SimCorMultRes Package. The R Journal 8:2, 79-91.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Simulating Correlated Binary and Multinomial Responses under 
#>          Marginal Model Specification: The SimCorMultRes Package},
#>     author = {Anestis Touloumis},
#>     year = {2016},
#>     journal = {The R Journal},
#>     volume = {8},
#>     number = {2},
#>     pages = {79-91},
#>     url = {https://journal.r-project.org/archive/2016/RJ-2016-034/index.html},
#>   }
```

References
==========

Touloumis, Anestis. 2016. “Simulating Correlated Binary and Multinomial Responses Under Marginal Model Specification: The SimCorMultRes Package.” *The R Journal* 8 (2): 79–91. <https://journal.r-project.org/archive/2016/RJ-2016-034/index.html>.
