 \name{rmult.bcl}
\alias{rmult.bcl}
\title{
Simulating Correlated Nominal Responses Conditional on a Marginal Baseline-Category Logit Model Specification
}
\description{
Simulates correlated nominal responses assuming a baseline-category logit model for the marginal probabilities.
}
\usage{
rmult.bcl(clsize, ncategories, betas, xformula, xdata, cor.matrix, rlatent = NULL)
}
\arguments{
  \item{clsize}{
integer indicating the common cluster size.
}
  \item{ncategories}{
integer indicating the number of nominal response categories. 
}
 \item{betas}{
numerical vector or matrix containing the value of the marginal regression parameter vector.
}
  \item{xformula}{
formula expression as in other marginal regression models but without including a response variable.  
}
  \item{xdata}{
optional data frame containing the variables provided in \code{xformula}.
} 
\item{cor.matrix}{
matrix indicating the correlation matrix of the multivariate normal distribution when the NORTA method is employed (\code{rlatent = NULL}).
}
 \item{rlatent}{
matrix with \code{(clsize * ncategories)} columns containing realizations of the latent random vectors when the NORTA method is not preferred. See details for more info.
}
}
\details{
The formulae are easier to read from either the Vignette or the Reference Manual (both available \href{https://CRAN.R-project.org/package=SimCorMultRes}{here}).

The assumed marginal baseline category logit model is 
\deqn{log \frac{Pr(Y_{it}=j |x_{it})}{Pr(Y_{it}=J |x_{it})}=(\beta_{tj0}-\beta_{tJ0}) + (\beta^{'}_{tj}-\beta^{'}_{tJ}) x_{it}=\beta^{*}_{tj0}+ \beta^{*'}_{tj} x_{it}}
For subject \eqn{i}, \eqn{Y_{it}} is the \eqn{t}-th nominal response and \eqn{x_{it}} is the associated covariates vector. Also \eqn{\beta_{tj0}} is the \eqn{j}-th category-specific intercept at the \eqn{t}-th measurement occasion and \eqn{\beta_{tj}} is the \eqn{j}-th category-specific regression parameter vector at the \eqn{t}-th measurement occasion.

The nominal response \eqn{Y_{it}} is obtained by extending the principle of maximum random utility (\cite{McFadden, 1974}) as suggested in \cite{Touloumis (2016)}.

\code{betas} should be provided as a numeric vector only when \eqn{\beta_{tj0}=\beta_{j0}} and \eqn{\beta_{tj}=\beta_j} for all \eqn{t}. Otherwise, \code{betas} must be provided as a numeric matrix with \code{clsize} rows such that the \eqn{t}-th row contains the value of (\eqn{\beta_{t10},\beta_{t1},\beta_{t20},\beta_{t2},...,\beta_{tJ0}, \beta_{tJ}}). In either case, \code{betas} should reflect the order of the terms implied by \code{xformula}.

The appropriate use of \code{xformula} is \code{xformula = ~ covariates}, where \code{covariates} indicate the linear predictor as in other marginal regression models.

The optional argument \code{xdata} should be provided in ``long'' format.

The NORTA method is the default option for simulating the latent random vectors denoted by \eqn{e^{NO}_{itj}} in \cite{Touloumis (2016)}. In this case, the algorithm forces \code{cor.matrix} to respect the assumption of choice independence. To import simulated values for the latent random vectors without utilizing the NORTA method, the user can employ the \code{rlatent} argument. In this case, row \eqn{i} corresponds to subject \eqn{i} and columns \eqn{(t-1)*\code{ncategories}+1,...,t*\code{ncategories}} should contain the realization of \eqn{e^{NO}_{it1},...,e^{NO}_{itJ}}, respectively, for \eqn{t=1,\ldots,\code{clsize}}. 
}
\value{
Returns a list that has components:
\item{Ysim}{the simulated nominal responses. Element (\eqn{i},\eqn{t}) represents the realization of \eqn{Y_{it}}.}
\item{simdata}{a data frame that includes the simulated response variables (y), the covariates specified by \code{xformula}, subjects' identities (id) and the corresponding measurement occasions (time).}
\item{rlatent}{the latent random variables denoted by \eqn{e^{NO}_{it}} in \cite{Touloumis (2016)}.}
}
\references{
Cario, M. C. and Nelson, B. L. (1997)  \emph{Modeling and generating random vectors with arbitrary marginal distributions and correlation matrix}. Technical Report, Department of Industrial Engineering and Management Sciences, Northwestern University, Evanston, Illinois.

Li, S. T. and Hammond, J. L. (1975)  Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transactions on Systems, Man and Cybernetics} \bold{5}, 557--561.   

McFadden, D. (1974) \emph{Conditional logit analysis of qualitative choice behavior}. New York: Academic Press, 105--142.

Touloumis, A. (2016) Simulating Correlated Binary and Multinomial Responses under Marginal Model Specification: The SimCorMultRes Package. \emph{The R Journal} \bold{8}, 79--91.

Touloumis, A., Agresti, A. and Kateri, M. (2013) GEE for multinomial responses using a local odds ratios parameterization. \emph{Biometrics} \bold{69}, 633--640.
}
\author{
Anestis Touloumis
}
\seealso{
\code{\link{rbin}} for simulating correlated binary responses, \code{\link{rmult.clm}} and \code{\link{rmult.crm}} for simulating correlated ordinal responses.
}
\examples{
## See Example 3.1 in the Vignette.
betas <- c(1, 3, 2, 1.25, 3.25, 1.75, 0.75, 2.75, 2.25, 0, 0, 0)
N <- 500
ncategories <- 4
clsize <- 3
set.seed(1)
x1 <- rep(rnorm(N), each = clsize)
x2 <- rnorm(N * clsize)
xdata <- data.frame(x1, x2)
library(SimCorMultRes)
cor.matrix <- kronecker(toeplitz(c(1, rep(0.95, clsize - 1))), diag(ncategories))
CorNorRes <- rmult.bcl(clsize = clsize, ncategories = ncategories, betas = betas, 
    xformula = ~x1 + x2, xdata = xdata, cor.matrix = cor.matrix)
suppressPackageStartupMessages(library("multgee"))
fit <- nomLORgee(y ~ x1 + x2, data = CorNorRes$simdata, id = id, repeated = time, 
    LORstr = "time.exch")
round(coef(fit), 2)
}
