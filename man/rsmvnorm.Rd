\name{rsmvnorm}
\alias{rsmvnorm}
\title{
Simulating Continuous Random Vectors from a Multivariate Normal Distribution
}
\description{
Utility function to simulate continuous random vectors from a multivariate normal distribution such that all marginal distributions are univariate standard normal.
}
\usage{
rsmvnorm(R, cor.matrix)
}
\arguments{
  \item{R}{
integer indicating the sample size.
}
  \item{cor.matrix}{
matrix indicating the correlation matrix of the multivariate normal distribution. 
}
}
\details{
Checks are made to ensure that \code{cor.matrix} is a positive definite correlation matrix. The positive definiteness of \code{cor.matrix} is assessed via eigenvalues.
}
\value{
Returns \code{R} random vectors of size \code{ncol(cor.matrix)}.  
}
\author{
Anestis Touloumis
}
\examples{
## Simulating 10000 bivariate random vectors with correlation parameter
## equal to 0.4.
set.seed(1)
R <- 10000
cor.matrix <- toeplitz(c(1, 0.4))
SimBivariateNormal <- rsmvnorm(R = R, cor.matrix = cor.matrix)
colMeans(SimBivariateNormal)
apply(SimBivariateNormal, 2, sd)
cor(SimBivariateNormal)
}
