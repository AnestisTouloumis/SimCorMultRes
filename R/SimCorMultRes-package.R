#' Simulating Correlated Multinomial Responses
#'
#' Functions to simulate correlated multinomial responses (three or more
#' nominal or ordinal response categories) and correlated binary responses
#' subject to a marginal model specification.
#'
#' The simulated correlated binary or multinomial responses are drawn as
#' realizations of a latent regression model for continuous random vectors with
#' the correlation structure expressed in terms of the latent correlation.
#'
#' For an ordinal response scale, the multinomial variables are simulated
#' conditional on a marginal cumulative link model
#' (\code{\link{rmult.clm}}), a marginal continuation-ratio model
#' (\code{\link{rmult.crm}}) or a marginal adjacent-category logit model
#' (\code{\link{rmult.acl}}).
#'
#' For a nominal response scale, the multinomial responses are simulated
#' conditional on a marginal baseline-category logit model
#' (\code{\link{rmult.bcl}}).
#'
#' Correlated binary responses are simulated using the function
#' \code{\link{rbin}}.
#'
#' The threshold approaches that give rise to the implemented marginal models
#' are fully described in \cite{Touloumis (2016)} and in the Vignette.
#'
#' The formulae are easier to read from either the Vignette or the Reference
#' Manual (both available
#' \href{https://CRAN.R-project.org/package=SimCorMultRes}{here}).
#'
#' @name SimCorMultRes-package
#' @aliases SimCorMultRes-package SimCorMultRes
#' @author Anestis Touloumis
#'
#' Maintainer: Anestis Touloumis \email{A.Touloumis@@brighton.ac.uk}
#' @references Cario, M. C. and Nelson, B. L. (1997) \emph{Modeling and
#' generating random vectors with arbitrary marginal distributions and
#' correlation matrix}. Technical Report, Department of Industrial Engineering
#' and Management Sciences, Northwestern University, Evanston, Illinois.
#'
#' Emrich, L. J. and Piedmonte, M. R. (1991) A method for generating
#' high-dimensional multivariate binary variates. \emph{The American
#' Statistician} \bold{45}, 302--304.
#'
#' Li, S. T. and Hammond, J. L. (1975) Generation of pseudorandom numbers with
#' specified univariate distributions and correlation coefficients. \emph{IEEE
#' Transactions on Systems, Man and Cybernetics} \bold{5}, 557--561.
#'
#' McCullagh, P. (1980) Regression models for ordinal data. \emph{Journal of
#' the Royal Statistical Society B} \bold{42}, 109--142.
#'
#' McFadden, D. (1974) \emph{Conditional logit analysis of qualitative choice
#' behavior}. New York: Academic Press, 105--142.
#'
#' Touloumis, A. (2016) Simulating Correlated Binary and Multinomial Responses
#' under Marginal Model Specification: The SimCorMultRes Package. \emph{The R
#' Journal} \bold{8}, 79--91.
#'
#' Touloumis, A., Agresti, A. and Kateri, M. (2013) GEE for multinomial
#' responses using a local odds ratios parameterization. \emph{Biometrics}
#' \bold{69}, 633--640.
#'
#' Tutz, G. (1991) Sequential models in categorical regression.
#' \emph{Computational Statistics & Data Analysis} \bold{11}, 275--295.
#' @importFrom evd qgumbel
#' @importFrom methods formalArgs
#' @importFrom stats as.formula formula get_all_vars model.frame model.matrix
#' na.omit pnorm qcauchy qlogis qunif rnorm terms toeplitz update
#' @keywords internal
"_PACKAGE"
