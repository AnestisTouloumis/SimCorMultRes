#' Simulating Continuous Random Vectors from a Multivariate Normal Distribution
#'
#' Utility function to simulate continuous random vectors from a multivariate
#' normal distribution such that all marginal distributions are univariate
#' standard normal.
#'
#' Checks are made to ensure that \code{cor.matrix} is a positive definite
#' correlation matrix. The positive definiteness of \code{cor.matrix} is
#' assessed via eigenvalues.
#'
#' @param R integer indicating the sample size.
#' @param cor.matrix matrix indicating the correlation matrix of the
#' multivariate normal distribution.
#' @return Returns \code{R} random vectors of size \code{ncol(cor.matrix)}.
#' @author Anestis Touloumis
#' @examples
#' ## Simulating 10000 bivariate random vectors with correlation parameter
#' ## equal to 0.4.
#' set.seed(1)
#' R <- 10000
#' cor.matrix <- toeplitz(c(1, 0.4))
#' SimBivariateNormal <- rsmvnorm(R = R, cor.matrix = cor.matrix)
#' colMeans(SimBivariateNormal)
#' apply(SimBivariateNormal, 2, sd)
#' cor(SimBivariateNormal)
#' @export
rsmvnorm <- function(R = R, cor.matrix = cor.matrix) {
    if (all.equal(R, as.integer(R)) != TRUE | R < 1)
        stop("'R' must be a positive integer")
    if (!is.numeric(cor.matrix))
        stop("'cor.matrix' must be numeric")
    cor.matrix <- as.matrix(cor.matrix)
    if (!isSymmetric(cor.matrix))
        stop("'cor.matrix' must be a symmetric matrix")
    if (any(diag(cor.matrix) != 1))
        stop("the diagonal elements of 'cor.matrix' must be equal to one")
    if (any(cor.matrix > 1) | any(cor.matrix < -1))
        stop("all the elements of 'cor.matrix' must be on [-1,1]")
    if (any(eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <=
        0))
        stop("'cor.matrix' must be a positive definite matrix")
    p <- ncol(cor.matrix)
    ans <- matrix(rnorm(R * p), R, p) %*% chol(cor.matrix)
    ans
}
