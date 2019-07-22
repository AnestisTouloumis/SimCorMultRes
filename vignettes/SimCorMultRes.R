## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  tidy = TRUE,
  collapse = TRUE,
  comment = "#>"
  )

## ---- tidy=TRUE----------------------------------------------------------
# parameter vector
betas <-  c(1, 3, 2, 1.25, 3.25, 1.75, 0.75, 2.75, 2.25, 0, 0, 0)
# sample size
sample_size <- 500 
# number of nominal response categories
categories_no <- 4  
# cluster size
cluster_size <- 3 
set.seed(1) 
# time-stationary covariate x_{i1}
x1 <- rep(rnorm(sample_size), each = cluster_size)
# time-varying covariate x_{it2}
x2 <- rnorm(sample_size * cluster_size) 
# create covariates dataframe 
xdata <- data.frame(x1, x2) 
set.seed(321)
library(SimCorMultRes)
# latent correlation matrix for the NORTA method
equicorrelation_matrix <- toeplitz(c(1, rep(0.95,cluster_size - 1)))
identity_matrix <- diag(categories_no)
latent_correlation_matrix <- kronecker(equicorrelation_matrix, identity_matrix) 
# simulation of clustered nominal responses
simulated_nominal_responses <- rmult.bcl(clsize = cluster_size,
                                         ncategories = categories_no, 
                                         betas = betas, xformula = ~ x1 + x2,
                                         xdata = xdata,
                                         cor.matrix = latent_correlation_matrix)
suppressPackageStartupMessages(library("multgee"))
# fitting a GEE model
nominal_gee_model <- nomLORgee(y ~ x1 + x2,
                               data = simulated_nominal_responses$simdata, id = id,
                               repeated = time, LORstr="time.exch")
# checking regression coefficients
round(coef(nominal_gee_model), 2)

## ------------------------------------------------------------------------
set.seed(12345)
# sample size
sample_size <- 500
# cluster size
cluster_size <- 4
# category-specific intercepts
beta_intercepts <- c(-1.5, -0.5, 0.5, 1.5)
# time-varying regression parameters associated with covariates
beta_coefficients <- matrix(c(1, 2, 3, 4), 4, 1)
# time-stationary covariate 
x <- rep(rnorm(sample_size), each = cluster_size)
# latent correlation matrix for the NORTA method
latent_correlation_matrix <- toeplitz(c(1, 0.85, 0.5, 0.15))
# simulation of ordinal responses
simulated_ordinal_responses <- rmult.clm(clsize = cluster_size,
                                         intercepts = beta_intercepts,
                                         betas = beta_coefficients,
                                         xformula = ~ x,
                                         cor.matrix = latent_correlation_matrix,
                                         link = "probit")
# first eight rows of the simulated dataframe
head(simulated_ordinal_responses$simdata, n = 8)

## ------------------------------------------------------------------------
set.seed(1)
# sample size
sample_size <- 500
# cluster size
cluster_size <- 4
# category-specific intercepts
beta_intercepts <- c(-1.5, -0.5, 0.5, 1.5)
# regression parameters associated with covariates
beta_coefficients <- 1
# time-varying covariate
x <- rnorm(sample_size * cluster_size)
# number of ordinal response categories
categories_no <- 5
# correlation matrix for the NORTA method
latent_correlation_matrix <- diag(1, (categories_no - 1) * cluster_size) +
  kronecker(toeplitz(c(0, rep(0.24, categories_no - 2))), matrix(1, cluster_size, cluster_size))
# simulation of ordinal responses
simulated_ordinal_responses <- rmult.crm(clsize = cluster_size,
                                         intercepts = beta_intercepts,
                                         betas = beta_coefficients,
                                         xformula = ~ x,
                                         cor.matrix = latent_correlation_matrix,
                                         link = "probit")
# first six clusters with ordinal responses
head(simulated_ordinal_responses$Ysim)

## ---- tidy=TRUE----------------------------------------------------------
# intercepts
beta_intercepts <- c(3, 2, 1)
# parameter vector
beta_coefficients <- c(1, 1)
# sample size
sample_size <- 500 
# cluster size
cluster_size <- 3 
set.seed(321) 
# time-stationary covariate x_{i1}
x1 <- rep(rnorm(sample_size), each = cluster_size)
# time-varying covariate x_{it2}
x2 <- rnorm(sample_size * cluster_size) 
# create covariates dataframe 
xdata <- data.frame(x1, x2) 
# correlation matrix for the NORTA method
equicorrelation_matrix <- toeplitz(c(1, rep(0.95, cluster_size - 1)))
identity_matrix <- diag(4)
latent_correlation_matrix <- kronecker(equicorrelation_matrix, identity_matrix) 
# simulation of clustered ordinal responses
simulated_ordinal_responses <- rmult.acl(clsize = cluster_size,
                                         intercepts = beta_intercepts,
                                         betas = beta_coefficients,
                                         xformula = ~ x1 + x2, xdata = xdata,
                                         cor.matrix = latent_correlation_matrix)
suppressPackageStartupMessages(library("multgee"))
# fitting a GEE model
ordinal_gee_model <- ordLORgee(y ~ x1 + x2,
                               data = simulated_ordinal_responses$simdata,
                               id = id, repeated = time, LORstr = "time.exch",
                               link = "acl")
# checking regression coefficients
round(coef(ordinal_gee_model), 2)

## ------------------------------------------------------------------------
set.seed(123)
# sample size
sample_size <- 100
# cluster size
cluster_size <- 4
# intercept
beta_intercepts <- 0
# regression parameter associated with the covariate
beta_coefficients <- 0.2
# correlation matrix for the NORTA method
latent_correlation_matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
# time-stationary covariate
x <- rep(rnorm(sample_size), each = cluster_size)
# simulation of clustered binary responses
simulated_binary_responses <- rbin(clsize = cluster_size, 
                                   intercepts = beta_intercepts, 
                                   betas = beta_coefficients, xformula = ~ x,
                                   cor.matrix = latent_correlation_matrix,
                                   link = "probit")
library(gee)
# fitting a GEE model
binary_gee_model <- gee(y ~ x, family = binomial("probit"), id = id,
                        data = simulated_binary_responses$simdata)
# checking the estimated coefficients
summary(binary_gee_model)$coefficients

## ------------------------------------------------------------------------
set.seed(8)
# simulation of epsilon variables
library(evd)
simulated_latent_variables1 <- rmvevd(sample_size, dep = sqrt(1 - 0.9),
                                      model = "log", d = cluster_size)
simulated_latent_variables2 <- rmvevd(sample_size, dep = sqrt(1 - 0.9),
                                      model = "log", d = cluster_size)
simulated_latent_variables <- simulated_latent_variables1 -
  simulated_latent_variables2
# simulation of clustered binary responses
simulated_binary_responses <- rbin(clsize = cluster_size,
                                   intercepts = beta_intercepts,
                                   betas = beta_coefficients, xformula = ~ x,
                                   rlatent = simulated_latent_variables)
# fitting a GEE model
binary_gee_model <- gee(y ~ x, family = binomial("logit"), id = id,
                        data = simulated_binary_responses$simdata)
# checking the estimated coefficients
summary(binary_gee_model)$coefficients

## ------------------------------------------------------------------------
set.seed(123)
# sample size
sample_size <- 5000
# cluster size
cluster_size <- 4
# intercept
beta_intercepts <- qnorm(0.8)
# pseudo-covariate
x <- rep(0, each = cluster_size * sample_size)
# regression parameter associated with the covariate
beta_coefficients <- 0
# correlation matrix for the NORTA method
latent_correlation_matrix <- diag(cluster_size)
# simulation of clustered binary responses
simulated_binary_responses <- rbin(clsize = cluster_size,
                                   intercepts = beta_intercepts, 
                                   betas = beta_coefficients,
                                   xformula = ~x,
                                   cor.matrix = latent_correlation_matrix,
                                   link = "probit")
library(gee)
# simulated marginal probabilities
colMeans(simulated_binary_responses$Ysim)

## ---- tidy=TRUE----------------------------------------------------------
# sample size
sample_size <- 5000
# cluster size
cluster_size <- 3 
# pseudo-covariate 
x <- rep(0, each = cluster_size * sample_size)
# parameter vector
betas <-  c(log(0.1/0.4), 0, log(0.2/0.4), 0, log(0.3/0.4), 0, 0, 0)
# number of nominal response categories
categories_no <- 4  
set.seed(1) 
# correlation matrix for the NORTA method
latent_correlation_matrix <- kronecker(diag(cluster_size), diag(categories_no)) 
# simulation of clustered nominal responses
simulated_nominal_responses <- rmult.bcl(clsize = cluster_size,
                                         ncategories = categories_no,
                                         betas = betas, xformula = ~ x,
                                         cor.matrix = latent_correlation_matrix)
# simulated marginal probabilities
apply(simulated_nominal_responses$Ysim, 2, table) / sample_size

## ---- comment=""---------------------------------------------------------
citation("SimCorMultRes")

