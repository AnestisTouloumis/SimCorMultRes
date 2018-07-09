#' Simulating Correlated Ordinal Responses Conditional on a Marginal Cumulative
#' Link Model Specification
#'
#' Simulates correlated ordinal responses assuming a cumulative link model for
#' the marginal probabilities.
#'
#' The formulae are easier to read from either the Vignette or the Reference
#' Manual (both available
#' \href{https://CRAN.R-project.org/package=SimCorMultRes}{here}).
#'
#' The assumed marginal cumulative link model is \deqn{Pr(Y_{it}\le j
#' |x_{it})=F(\beta_{tj0} +\beta^{'}_{t} x_{it})} where \eqn{F} is the
#' cumulative distribution function determined by \code{link}. For subject
#' \eqn{i}, \eqn{Y_{it}} is the \eqn{t}-th ordinal response and \eqn{x_{it}} is
#' the associated covariates vector. Finally, \eqn{\beta_{tj0}} is the
#' \eqn{j}-th category-specific intercept at the \eqn{t}-th measurement
#' occasion and \eqn{\beta_{tj}} is the \eqn{j}-th category-specific regression
#' parameter vector at the \eqn{t}-th measurement occasion.
#'
#' The ordinal response \eqn{Y_{it}} is obtained by extending the approach of
#' \cite{McCullagh (1980)} as suggested in \cite{Touloumis (2016)}.
#'
#' When \eqn{\beta_{tj0}=\beta_{j0}} for all \eqn{t}, then \code{intercepts}
#' should be provided as a numerical vector. Otherwise, \code{intercepts} must
#' be a numerical matrix such that row \eqn{t} contains the category-specific
#' intercepts at the \eqn{t}-th measurement occasion.
#'
#' \code{betas} should be provided as a numeric vector only when
#' \eqn{\beta_{t}=\beta} for all \eqn{t}. Otherwise, \code{betas} must be
#' provided as a numeric matrix with \code{clsize} rows such that the
#' \eqn{t}-th row contains the value of \eqn{\beta_{t}}. In either case,
#' \code{betas} should reflect the order of the terms implied by
#' \code{xformula}.
#'
#' The appropriate use of \code{xformula} is \code{xformula = ~ covariates},
#' where \code{covariates} indicate the linear predictor as in other marginal
#' regression models.
#'
#' The optional argument \code{xdata} should be provided in ``long'' format.
#'
#' The NORTA method is the default option for simulating the latent random
#' vectors denoted by \eqn{e^{O1}_{it}} in \cite{Touloumis (2016)}. To import
#' simulated values for the latent random vectors without utilizing the NORTA
#' method, the user can employ the \code{rlatent} argument. In this case,
#' element (\eqn{i,t}) of \code{rlatent} represents the realization of
#' \eqn{e^{O1}_{it}}.
#'
#' @param clsize integer indicating the common cluster size.
#' @param intercepts numerical vector or matrix containing the intercepts of
#' the marginal cumulative link model.
#' @param betas numerical vector or matrix containing the value of the marginal
#' regression parameter vector associated with the covariates (i.e., excluding
#' \code{intercepts}).
#' @param xformula formula expression as in other marginal regression models
#' but without including a response variable.
#' @param xdata optional data frame containing the variables provided in
#' \code{xformula}.
#' @param link character string indicating the link function in the marginal
#' cumulative link model. Options include \code{"probit"}, \code{"logit"},
#' \code{"cloglog"} or \code{"cauchit"}. Required when \code{rlatent = NULL}.
#' @param cor.matrix matrix indicating the correlation matrix of the
#' multivariate normal distribution when the NORTA method is employed
#' (\code{rlatent = NULL}).
#' @param rlatent matrix with \code{clsize} columns containing realizations of
#' the latent random vectors when the NORTA method is not preferred. See
#' details for more info.
#' @return Returns a list that has components: \item{Ysim}{the simulated
#' ordinal responses. Element (\eqn{i},\eqn{t}) represents the realization of
#' \eqn{Y_{it}}.} \item{simdata}{a data frame that includes the simulated
#' response variables (y), the covariates specified by \code{xformula},
#' subjects' identities (id) and the corresponding measurement occasions
#' (time).} \item{rlatent}{the latent random variables denoted by
#' \eqn{e^{O1}_{it}} in \cite{Touloumis (2016)}.}
#' @author Anestis Touloumis
#' @seealso \code{\link{rmult.bcl}} for simulating correlated nominal
#' responses, \code{\link{rmult.crm}} for simulating correlated ordinal
#' responses and \code{\link{rbin}} for simulating correlated binary responses.
#' @references Cario, M. C. and Nelson, B. L. (1997) \emph{Modeling and
#' generating random vectors with arbitrary marginal distributions and
#' correlation matrix}. Technical Report, Department of Industrial Engineering
#' and Management Sciences, Northwestern University, Evanston, Illinois.
#'
#' Li, S. T. and Hammond, J. L. (1975) Generation of pseudorandom numbers with
#' specified univariate distributions and correlation coefficients. \emph{IEEE
#' Transactions on Systems, Man and Cybernetics} \bold{5}, 557--561.
#'
#' McCullagh, P. (1980) Regression models for ordinal data. \emph{Journal of
#' the Royal Statistical Society B} \bold{42}, 109--142.
#'
#' Touloumis, A. (2016) Simulating Correlated Binary and Multinomial Responses
#' under Marginal Model Specification: The SimCorMultRes Package. \emph{The R
#' Journal} \bold{8}, 79--91.
#'
#' Touloumis, A., Agresti, A. and Kateri, M. (2013) GEE for multinomial
#' responses using a local odds ratios parameterization. \emph{Biometrics}
#' \bold{69}, 633--640.
#' @examples
#' ## See Example 3.2 in the Vignette.
#' set.seed(12345)
#' N <- 500
#' clsize <- 4
#' intercepts <- c(-1.5, -0.5, 0.5, 1.5)
#' betas <- matrix(c(1, 2, 3, 4), 4, 1)
#' x <- rep(rnorm(N), each = clsize)
#' cor.matrix <- toeplitz(c(1, 0.85, 0.5, 0.15))
#' CorOrdRes <- rmult.clm(clsize = clsize, intercepts = intercepts, betas = betas,
#'     xformula = ~x, cor.matrix = cor.matrix, link = "probit")
#' head(CorOrdRes$simdata, n = 8)
#'
#' ## Same sampling scheme except that the parameter vector is now
#' ## time-stationary.
#' set.seed(12345)
#' x <- rep(rnorm(N), each = clsize)
#' CorOrdRes <- rmult.clm(clsize = clsize, betas = 1, xformula = ~x, cor.matrix = toeplitz(c(1,
#'     0.85, 0.5, 0.15)), intercepts = c(-1.5, -0.5, 0.5, 1.5), link = "probit")
#' ## Fit a GEE model (Touloumis et al., 2013) to estimate the regression
#' ## coefficients.
#' library(multgee)
#' fitmod <- ordLORgee(y ~ x, id = id, repeated = time, link = "probit", data = CorOrdRes$simdata)
#' coef(fitmod)
#'
#' @export
rmult.clm <- function(clsize = clsize, intercepts = intercepts, betas = betas,
    xformula = formula(xdata), xdata = parent.frame(), link = "logit",
    cor.matrix = cor.matrix, rlatent = NULL) {
    if (all.equal(clsize, as.integer(clsize)) != TRUE | clsize < 2)
        stop("'clsize' must be a positive integer greater than or equal to two")
    if (
      !(is.vector(intercepts) & !is.list(intercepts)) & !is.matrix(intercepts))
        stop("'intercepts' must be a vector or a matrix")
    if (!is.numeric(intercepts))
        stop("'intercepts' must be numeric")
    if (is.vector(intercepts) & !is.list(intercepts)) {
        if (length(intercepts) == 1)
            stop("'intercepts' must have at least 2 elements")
        if (any(diff(intercepts) <= 0))
            stop("'intercepts' must be increasing")
        ncategories <- length(intercepts) + 1
        intercepts <- matrix(intercepts, clsize, ncategories - 1, TRUE)
        intercepts <- cbind(-Inf, intercepts, Inf)
    } else {
        ncategories <- ncol(intercepts) + 1
        intercepts <- cbind(-Inf, intercepts, Inf)
        for (i in 1:clsize) {
            if (any(diff(intercepts[i, ]) <= 0))
                stop("'intercepts' must be increasing at each row")
        }
    }
    if (!(is.vector(betas) & !is.list(betas)) & !is.matrix(betas))
        stop("'betas' must be a vector or a matrix")
    if (!is.numeric(betas))
        stop("'betas' must be numeric")
    if (is.vector(betas) & !is.list(betas)) {
        betas <- rep(betas, clsize)
    } else {
        if (nrow(betas) != clsize)
            stop("The number of rows in 'betas' should be equal to 'clsize'")
        betas <- c(t(betas))
    }
    lpformula <- as.formula(xformula)
    if (length(paste0(attr(terms(lpformula), "variables"))) == 1)
        stop("No covariates were found in 'formula' ")
    if (attr(terms(lpformula), "intercept") == 0) {
      lpformula <- update(lpformula, ~. + 1)
    }
    Xmat <- model.matrix(lpformula, data = xdata)
    Xmat <- matrix(Xmat[, -1], ncol = ncol(Xmat) - 1)
    if (length(betas) != (clsize) * ncol(Xmat))
        stop("The length of 'betas' does not match with the number of covariates")
    lin.pred <- matrix(betas, nrow = nrow(Xmat), ncol = ncol(Xmat),
                       byrow = TRUE) * Xmat
    lin.pred <- matrix(rowSums(lin.pred), ncol = clsize, byrow = TRUE)
    lin.pred <- as.matrix(lin.pred)
    R <- nrow(lin.pred)
    if (is.null(rlatent)) {
        if (length(link) != 1)
            stop("The length of 'link' must be one")
        links <- c("probit", "logit", "cloglog", "cauchit")
        if (!is.element(link, links))
            stop("'link' must be 'probit','logit','cloglog' and/or 'cauchit'")
        distr <- switch(link, probit = "qnorm", logit = "qlogis",
                        cloglog = "qgumbel", cauchit = "qcauchy")
        if (!is.numeric(cor.matrix))
            stop("'cor.matrix' must be numeric")
        if (!is.matrix(cor.matrix))
            stop("'cor.matrix' must be matrix")
        if (ncol(cor.matrix) != clsize | nrow(cor.matrix) != clsize)
            stop("'cor.matrix' must be a ", clsize, "x", clsize, " matrix")
        if (!isSymmetric(cor.matrix))
            stop("'cor.matrix' must be a symmetric matrix")
        if (any(diag(cor.matrix) != 1))
            stop("the diagonal elements of 'cor.matrix' must be one")
        if (any(cor.matrix > 1) | any(cor.matrix < -1))
            stop("all the elements of 'cor.matrix' must be on [-1,1]")
        if (any(
          eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <= 0))
            stop("'cor.matrix' must be positive definite")
        if (length(distr) == 1)
            err <- rnorta(R, cor.matrix, rep(distr, clsize))
        if (distr == "qgumbel")
            err <- -err
    } else {
        if (!is.matrix(rlatent))
            stop("'rlatent' must be a matrix")
        if (!is.numeric(rlatent))
            stop("'rlatent' must be numeric")
        if (ncol(rlatent) != clsize | nrow(rlatent) != R)
            stop("'rlatent' must be a ", R, "x", clsize, " matrix")
        cor.matrix <- NULL
        err <- rlatent
    }
    U <- -lin.pred + err
    Ysim <- matrix(0, R, clsize)
    for (i in 1:clsize) Ysim[, i] <-
      cut(U[, i], intercepts[i, ], labels = FALSE)
    id <- rep(1:R, each = clsize)
    time <- rep(1:clsize, R)
    y <- c(t(Ysim))
    rownames(Ysim) <- rownames(err) <- paste("i", 1:R, sep = "=")
    colnames(Ysim) <- colnames(err) <- paste("t", 1:clsize, sep = "=")
    simdata <- data.frame(y, model.frame(formula = lpformula, data = xdata),
        id, time)
    list(Ysim = Ysim, simdata = simdata, rlatent = err)
}
