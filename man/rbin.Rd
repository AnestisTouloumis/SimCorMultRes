% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbin.R
\name{rbin}
\alias{rbin}
\title{Simulating Correlated Binary Responses Conditional on a Marginal Model
Specification}
\usage{
rbin(clsize = clsize, intercepts = intercepts, betas = betas,
  xformula = formula(xdata), xdata = parent.frame(), link = "logit",
  cor.matrix = cor.matrix, rlatent = NULL)
}
\arguments{
\item{clsize}{integer indicating the common cluster size.}

\item{intercepts}{numerical (or numeric vector of length \code{clsize})
containing the intercept(s) of the marginal model.}

\item{betas}{numerical vector or matrix containing the value of the marginal
regression parameter vector associated with the covariates (i.e., excluding
\code{intercepts}).}

\item{xformula}{formula expression as in other marginal regression models
but without including a response variable.}

\item{xdata}{optional data frame containing the variables provided in
\code{xformula}.}

\item{link}{character string indicating the link function in the marginal
model. Options include \code{'probit'}, \code{'logit'}, \code{'cloglog'},
\code{'cauchit'} or \code{'identity'}. Required when \code{rlatent = NULL}.}

\item{cor.matrix}{matrix indicating the correlation matrix of the
multivariate normal distribution when the NORTA method is employed
(\code{rlatent = NULL}).}

\item{rlatent}{matrix with \code{clsize} columns containing realizations of
the latent random vectors when the NORTA method is not preferred. See
details for more info.}
}
\value{
Returns a list that has components: \item{Ysim}{the simulated binary
responses. Element (\eqn{i},\eqn{t}) represents the realization of
\eqn{Y_{it}}.} \item{simdata}{a data frame that includes the simulated
response variables (y), the covariates specified by \code{xformula},
subjects' identities (id) and the corresponding measurement occasions
(time).} \item{rlatent}{the latent random variables denoted by
\eqn{e^{B}_{it}} in \cite{Touloumis (2016)}.}
}
\description{
Simulates correlated binary responses assuming a regression model for the
marginal probabilities.
}
\details{
The formulae are easier to read from either the Vignette or the Reference
Manual (both available
\href{https://CRAN.R-project.org/package=SimCorMultRes}{here}).

The assumed marginal model is \deqn{Pr(Y_{it} = 1 |x_{it})=F(\beta_{t0}
+\beta^{'}_{t} x_{it})} where \eqn{F} is the cumulative distribution
function determined by \code{link}. For subject \eqn{i}, \eqn{Y_{it}} is the
\eqn{t}-th binary response and \eqn{x_{it}} is the associated covariates
vector. Finally, \eqn{\beta_{t0}} and \eqn{\beta_{t}} are the intercept and
regression parameter vector at the \eqn{t}-th measurement occasion.

The binary response \eqn{Y_{it}} is obtained by extending the approach of
\cite{Emrich and Piedmonte (1991)} as suggested in \cite{Touloumis (2016)}.

When \eqn{\beta_{t0}=\beta_{0}} for all \eqn{t}, then \code{intercepts}
should be provided as a single number. Otherwise, \code{intercepts} must be
provided as a numeric vector such that the \eqn{t}-th element corresponds to
the intercept at measurement occasion \eqn{t}.

\code{betas} should be provided as a numeric vector only when
\eqn{\beta_{t}=\beta} for all \eqn{t}. Otherwise, \code{betas} must be
provided as a numeric matrix with \code{clsize} rows such that the
\eqn{t}-th row contains the value of \eqn{\beta_{t}}. In either case,
\code{betas} should reflect the order of the terms implied by
\code{xformula}.

The appropriate use of \code{xformula} is \code{xformula = ~ covariates},
where \code{covariates} indicate the linear predictor as in other marginal
regression models.

The optional argument \code{xdata} should be provided in ``long'' format.

The NORTA method is the default option for simulating the latent random
vectors denoted by \eqn{e^{B}_{it}} in \cite{Touloumis (2016)}. To import
simulated values for the latent random vectors without utilizing the NORTA
method, the user can employ the \code{rlatent} argument. In this case,
element (\eqn{i,t}) of \code{rlatent} represents the realization of
\eqn{e^{B}_{it}}.
}
\examples{
## See Example 3.5 in the Vignette.
set.seed(123)
sample_size <- 5000
cluster_size <- 4
beta_intercepts <- 0
beta_coefficients <- 0.2
latent_correlation_matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
x <- rep(rnorm(sample_size), each = cluster_size)
simulated_binary_dataset <- rbin(clsize = cluster_size,
  intercepts = beta_intercepts, betas = beta_coefficients,
  xformula = ~x, cor.matrix = latent_correlation_matrix, link = "probit")
library(gee)
binary_gee_model <- gee(y ~ x, family = binomial("probit"), id = id,
  data = simulated_binary_dataset$simdata)
summary(binary_gee_model)$coefficients

## See Example 3.6 in the Vignette.
set.seed(8)
library(evd)
simulated_latent_variables1 <- rmvevd(sample_size, dep = sqrt(1 - 0.9),
  model = "log", d = cluster_size)
  simulated_latent_variables2 <- rmvevd(sample_size, dep = sqrt(1 - 0.9),
  model = "log", d = cluster_size)
simulated_latent_variables <- simulated_latent_variables1 -
  simulated_latent_variables2
simulated_binary_dataset <- rbin(clsize = cluster_size,
  intercepts = beta_intercepts, betas = beta_coefficients,
  xformula = ~x, rlatent = simulated_latent_variables)
binary_gee_model <- gee(y ~ x, family = binomial("logit"), id = id,
  data = simulated_binary_dataset$simdata)
summary(binary_gee_model)$coefficients
}
\references{
Cario, M. C. and Nelson, B. L. (1997) \emph{Modeling and
generating random vectors with arbitrary marginal distributions and
correlation matrix}. Technical Report, Department of Industrial Engineering
and Management Sciences, Northwestern University, Evanston, Illinois.

Emrich, L. J. and Piedmonte, M. R. (1991) A method for generating
high-dimensional multivariate binary variates. \emph{The American
Statistician} \bold{45}, 302--304.

Li, S. T. and Hammond, J. L. (1975) Generation of pseudorandom numbers with
specified univariate distributions and correlation coefficients. \emph{IEEE
Transactions on Systems, Man and Cybernetics} \bold{5}, 557--561.

Touloumis, A. (2016) Simulating Correlated Binary and Multinomial Responses
under Marginal Model Specification: The SimCorMultRes Package. \emph{The R
Journal} \bold{8}, 79--91.
}
\seealso{
\code{\link{rmult.bcl}} for simulating correlated nominal
responses, \code{\link{rmult.clm}}, \code{\link{rmult.crm}} and
\code{\link{rmult.acl}} for simulating correlated ordinal responses.
}
\author{
Anestis Touloumis
}
