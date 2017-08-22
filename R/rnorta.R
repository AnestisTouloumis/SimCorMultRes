#' Simulating Random Vectors using the NORTA Method
#' 
#' Utility function to simulate random vectors with predefined marginal
#' distributions via the NORTA method.
#' 
#' Checks are made to ensure that \code{cor.matrix} is a positive definite
#' correlation matrix. The positive definiteness of \code{cor.matrix} is
#' assessed via eigenvalues.
#' 
#' The \eqn{t}-th character string in \code{distr} indicates the quantile
#' function of the \eqn{t}-th marginal distribution. See
#' \code{\link{Distributions}} for the most common distributions. Quantile
#' functions supported by other R packages are allowed provided that these
#' packages have been uploaded first. However, note that no checks are made to
#' ensure that the character strings in \code{distr} correspond to valid names
#' of quantile functions.
#' 
#' If \code{qparameters = NULL} then the default parameter values for the
#' quantile functions specified by \code{distr} are used. Otherwise,
#' \code{qparameters} should be provided as a list of \code{ncol(cor.matrix)}
#' lists such that the \eqn{t}-th list contains the desired parameter values of
#' the \eqn{t}-th quantile function.
#' 
#' @param R integer indicating the sample size.
#' @param cor.matrix matrix indicating the correlation matrix of the
#' multivariate normal distribution employed in the NORTA method.
#' @param distr character string vector of length \code{ncol(cor.matrix)}
#' naming the quantile functions of the desired marginal distributions.
#' @param qparameters list of \code{ncol(cor.matrix)} lists indicating the
#' parameter values of the quantile functions specified by \code{distr}.
#' @return Returns \code{R} random vectors of size \code{ncol(cor.matrix)} with
#' marginal distributions specified by \code{distr} (and \code{qparameters}).
#' @author Anestis Touloumis
#' @references Cario, M. C. and Nelson, B. L. (1997) \emph{Modeling and
#' generating random vectors with arbitrary marginal distributions and
#' correlation matrix}. Technical Report, Department of Industrial Engineering
#' and Management Sciences, Northwestern University, Evanston, Illinois.
#' 
#' Li, S. T. and Hammond, J. L. (1975) Generation of pseudorandom numbers with
#' specified univariate distributions and correlation coefficients. \emph{IEEE
#' Transactions on Systems, Man and Cybernetics} \bold{5}, 557--561.
#' 
#' Touloumis, A. (2016) Simulating Correlated Binary and Multinomial Responses
#' under Marginal Model Specification: The SimCorMultRes Package. \emph{The R
#' Journal} \bold{8}, 79--91.
#' @examples
#' ## An example with standard logistic as marginal distribution.
#' set.seed(1)
#' R <- 1000
#' LatentCorrelation <- toeplitz(c(1, rep(0.8, 2)))
#' LatentCorrelation
#' CommonMarginals <- rep("qlogis", 3)
#' SimLogistic <- rnorta(R = R, cor.matrix = LatentCorrelation, distr = CommonMarginals)
#' 
#' ## The following lines exemplify the NORTA method.
#' set.seed(1)
#' SimNormal <- rsmvnorm(R = R, cor.matrix = LatentCorrelation)
#' all(SimLogistic == qlogis(pnorm(SimNormal)))
#' 
#' ## Change the marginal distributions to standard normal, standard
#' ## logistic and standard extreme value distribution.
#' set.seed(1)
#' DiffMarginals <- c("qnorm", "qlogis", "qgumbel")
#' SimDiffMars <- rnorta(R = R, cor.matrix = LatentCorrelation, distr = DiffMarginals)
#' cor(SimDiffMars)
#' colMeans(SimDiffMars)
#' apply(SimDiffMars, 2, sd)
#' 
#' ## Same as above but using parameter values other than the default ones.
#' set.seed(1)
#' qpars <- list(c(mean = 1, sd = 9), c(location = 2, scale = 1), c(loc = 3, 
#'     scale = 1))
#' SimDiffMars2 <- rnorta(R = R, cor.matrix = LatentCorrelation, distr = DiffMarginals, 
#'     qparameters = qpars)
#' cor(SimDiffMars2)
#' colMeans(SimDiffMars2)
#' apply(SimDiffMars2, 2, sd)
#' @export
rnorta <- function(R = R, cor.matrix = cor.matrix, distr = distr, 
                   qparameters = NULL) {
    if (all.equal(R, as.integer(R)) != TRUE | R < 1) 
        stop("'R' must be a positive integer")
    QuantileFunctions <- as.character(distr)
    ans <- rsmvnorm(R = R, cor.matrix = cor.matrix)
    if (length(QuantileFunctions) != ncol(cor.matrix)) 
        stop("'distr' must be a ", ncol(cor.matrix),
             "-variate vector of strings naming a valid quantile function")
    if (!is.null(qparameters)) {
        qparameters <- as.list(qparameters)
        if (length(qparameters) != ncol(cor.matrix)) 
            stop("'qparameters' must be provided as a list of length ", 
                ncol(cor.matrix))
    }
    ans <- pnorm(ans)
    for (i in seq_len(ncol(cor.matrix))) {
        QuantileFunction <- get(QuantileFunctions[i], mode = "function")
        if (!is.function(QuantileFunction)) 
            stop("Character string ", i, 
                 " in `distr' does not correspond to a valid function")
        if (!is.null(qparameters)) 
            formals(QuantileFunction)[pmatch(names(qparameters[[i]]),
                                      methods::formalArgs(QuantileFunction))] <- 
            qparameters[[i]]
        ans[, i] <- QuantileFunction(ans[, i])
    }
    ans
}
