#' Simulating Correlated Ordinal Responses Conditional on a Marginal
#' Continuation-Ratio Model Specification
#'
#' Simulates correlated ordinal responses assuming a continuation-ratio model
#' for the marginal probabilities.
#'
#' The formulae are easier to read from either the Vignette or the Reference
#' Manual (both available
#' \href{https://CRAN.R-project.org/package=SimCorMultRes}{here}).
#'
#' The assumed marginal continuation-ratio model is \deqn{Pr(Y_{it}=j |Y_{it}
#' \ge j,x_{it})=F(\beta_{tj0} +\beta^{'}_{t} x_{it})} where \eqn{F} is the
#' cumulative distribution function determined by \code{link}. For subject
#' \eqn{i}, \eqn{Y_{it}} is the \eqn{t}-th multinomial response and
#' \eqn{x_{it}} is the associated covariates vector. Finally, \eqn{\beta_{tj0}}
#' is the \eqn{j}-th category-specific intercept at the \eqn{t}-th measurement
#' occasion and \eqn{\beta_{tj}} is the \eqn{j}-th category-specific regression
#' parameter vector at the \eqn{t}-th measurement occasion.
#'
#' The ordinal response \eqn{Y_{it}} is determined by extending the latent
#' variable threshold approach of \cite{Tutz (1991)} as suggested in
#' \cite{Touloumis (2016)}.
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
#' vectors denoted by \eqn{e^{O2}_{itj}} in \cite{Touloumis (2016)}. In this
#' case, the algorithm forces \code{cor.matrix} to respect the local
#' independence assumption. To import simulated values for the latent random
#' vectors without utilizing the NORTA method, the user can employ the
#' \code{rlatent} argument. In this case, row \eqn{i} corresponds to subject
#' \eqn{i} and columns
#' \eqn{(t-1)*\code{ncategories}+1,...,t*\code{ncategories}} should contain the
#' realization of \eqn{e^{O2}_{it1},...,e^{O2}_{itJ}}, respectively, for
#' \eqn{t=1,\ldots,\code{clsize}}.
#'
#' @param clsize integer indicating the common cluster size.
#' @param intercepts numerical vector or matrix containing the intercepts of
#' the marginal continuation-ratio model.
#' @param betas numerical vector or matrix containing the value of the marginal
#' regression parameter vector associated with the covariates (i.e., excluding
#' \code{intercepts}).
#' @param xformula formula expression as in other marginal regression models
#' but without including a response variable.
#' @param xdata optional data frame containing the variables provided in
#' \code{xformula}.
#' @param link character string indicating the link function of the marginal
#' continuation-ratio model. Options include \code{"probit"}, \code{"logit"},
#' \code{"cloglog"} or \code{"cauchit"}. Required when \code{rlatent = NULL}.
#' @param cor.matrix matrix indicating the correlation matrix of the
#' multivariate normal distribution when the NORTA method is employed
#' (\code{rlatent = NULL}).
#' @param rlatent matrix with \code{clsize} rows and \code{ncategories} columns
#' containing realizations of the latent random vectors when the NORTA method
#' is not employed. See details for more info.
#' @return Returns a list that has components: \item{Ysim}{the simulated
#' ordinal responses. Element (\eqn{i},\eqn{t}) represents the realization of
#' \eqn{Y_{it}}.} \item{simdata}{a data frame that includes the simulated
#' response variables (y), the covariates specified by \code{xformula},
#' subjects' identities (id) and the corresponding measurement occasions
#' (time).} \item{rlatent}{the latent random variables denoted by
#' \eqn{e^{O2}_{it}} in \cite{Touloumis (2016)}.}
#' @author Anestis Touloumis
#' @seealso \code{\link{rmult.bcl}} for simulating correlated nominal
#' responses, \code{\link{rmult.clm}} for simulating correlated ordinal
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
#' Touloumis, A. (2016) Simulating Correlated Binary and Multinomial Responses
#' under Marginal Model Specification: The SimCorMultRes Package. \emph{The R
#' Journal (forthcoming)}.
#'
#' Tutz, G. (1991) Sequential models in categorical regression.
#' \emph{Computational Statistics & Data Analysis} \bold{11}, 275--295.
#' @examples
#' ## See Example 3.3 in the Vignette.
#' set.seed(1)
#' N <- 500
#' clsize <- 4
#' intercepts <- c(-1.5, -0.5, 0.5, 1.5)
#' betas <- 1
#' x <- rnorm(N * clsize)
#' ncategories <- 5
#' cor.matrix <- diag(1, (ncategories - 1) * clsize) + kronecker(toeplitz(c(0,
#'     rep(0.24, ncategories - 2))), matrix(1, clsize, clsize))
#' CorOrdRes <- rmult.crm(clsize = clsize, intercepts = intercepts, betas = betas,
#'     xformula = ~x, cor.matrix = cor.matrix, link = "probit")
#' head(CorOrdRes$Ysim)
#'
#'@export
rmult.crm <- function(clsize = clsize, intercepts = intercepts, betas = betas,
    xformula = formula(xdata), xdata = parent.frame(), link = "logit",
    cor.matrix = cor.matrix, rlatent = NULL) {
    if (all.equal(clsize, as.integer(clsize)) != TRUE | clsize < 2)
        stop("'clsize' must be a positive integer greater than or equal to two")
    if (!(is.atomic(betas) || is.list(betas)) & !is.matrix(betas))
        stop("'betas' must be a vector or a matrix")
    if (!is.numeric(betas))
        stop("'betas' must be numeric")
    if ((is.atomic(betas) || is.list(betas))) {
        betas <- rep(betas, clsize)
    } else {
        if (nrow(betas) != clsize)
            stop("The number of rows in 'betas' should be equal to 'clsize'")
        betas <- c(t(betas))
    }
    lpformula <- as.formula(xformula)
    if (length(paste0(attr(terms(lpformula), "variables"))) == 1)
        stop("No covariates were found in 'formula' ")
    Xmat <- model.matrix(lpformula, data = xdata)
    if (attr(terms(lpformula), "intercept") == 1) {
        Xmat <- matrix(Xmat[, -1], ncol = ncol(Xmat) - 1)
    } else {
        att <- attr(Xmat, "assign")
        factor.columns <- unique(att[duplicated(att)])
        if (length(factor.columns) >= 1) {
            Xmat <- Xmat[, -match(factor.columns[1], att)]
            Xmat <- data.matrix(Xmat)
        }
    }
    if (length(betas) != (clsize) * ncol(Xmat))
        stop("The length of 'betas' does not match with the number of covariates")
    lin.pred <- matrix(betas, nrow = nrow(Xmat), ncol = ncol(Xmat),
                       byrow = TRUE) * Xmat
    lin.pred <- matrix(rowSums(lin.pred), ncol = clsize, byrow = TRUE)
    lin.pred <- as.matrix(lin.pred)
    R <- nrow(lin.pred)
    if (
      !(is.vector(intercepts) & !is.list(intercepts)) & !is.matrix(intercepts))
        stop("'intercepts' must be a vector or a matrix")
    if (!is.numeric(intercepts))
        stop("'intercepts' must be numeric")
    if (is.vector(intercepts) & !is.list(intercepts)) {
        if (length(intercepts) == 1)
            stop("'intercepts' must have at least 2 elements")
        ncategories <- length(intercepts) + 1
        intercepts <- matrix(intercepts, R, clsize * (ncategories - 1),
            byrow = TRUE)
    } else {
        ncategories <- ncol(intercepts) + 1
        intercepts <- matrix(t(intercepts), R, clsize * (ncategories -
            1), byrow = TRUE)
    }
    lin.pred.extended <- t(apply(lin.pred, 1,
                                 function(x) rep(x, each = ncategories - 1)))
    if (is.null(rlatent)) {
        if (length(link) != 1)
            stop("The length of 'link' must be one")
        links <- c("probit", "logit", "cloglog", "cauchit")
        if (!is.element(link, links))
            stop("'link' must be 'probit','logit','cloglog' and/or 'cauchit'")
        distr <- switch(link, probit = "qnorm", logit = "qlogis",
                        cloglog = "qgumbel",
            cauchit = "qcauchy")
        if (!is.numeric(cor.matrix))
            stop("'cor.matrix' must be numeric")
        if (!is.matrix(cor.matrix))
            stop("'cor.matrix' must be matrix")
        cor.matrix <- as.matrix(cor.matrix)
        dimcor <- clsize * (ncategories - 1)
        if (ncol(cor.matrix) != dimcor)
            stop("'cor.matrix' must be a ", dimcor, "x", dimcor, " matrix")
        if (!isSymmetric(cor.matrix))
            stop("'cor.matrix' must be symmetric")
        for (i in 1:clsize) {
            diag.index <- 1:(ncategories - 1) + (i - 1) * (ncategories -
                1)
            cor.matrix[diag.index, diag.index] <- diag(1, ncategories -
                1)
        }
        if (any(cor.matrix > 1) | any(cor.matrix < -1))
            stop("all the elements of 'cor.matrix' must be on [-1,1]")
        if (any(
          eigen(cor.matrix, symmetric = TRUE, only.values = TRUE)$values <= 0))
            stop("'cor.matrix' must be positive definite")
        err <- rnorta(R, cor.matrix, rep(distr, ncol(cor.matrix)))
        index <- distr == "qgumbel"
        if (distr == "qgumbel")
            err <- -err
    } else {
        if (!is.matrix(rlatent))
            stop("'rlatent' must be a matrix")
        if (!is.numeric(rlatent))
            stop("'rlatent' must be numeric")
        if (ncol(rlatent) != ncol(lin.pred.extended) | nrow(rlatent) !=
            R)
            stop("'rlatent' must be a ", R, "x", ncol(lin.pred.extended),
                " matrix")
        cor.matrix <- NULL
        err <- rlatent
    }
    U <- err - lin.pred.extended
    Ysim <- matrix(as.numeric(t(U <= intercepts)), R * clsize, ncategories -
        1, TRUE)
    for (i in 1:(ncategories - 1)) Ysim[, i] <- ifelse(Ysim[, i] == 1,
        i, ncategories)
    Ysim <- apply(Ysim, 1, min)
    Ysim <- matrix(Ysim, R, clsize, byrow = TRUE)
    id <- rep(1:R, each = clsize)
    time <- rep(1:clsize, R)
    y <- c(t(Ysim))
    rownames(Ysim) <- rownames(err) <- paste("i", 1:R, sep = "=")
    colnames(Ysim) <- paste("t", 1:clsize, sep = "=")
    colnames(err) <- paste("t=", rep(1:clsize, each = ncategories - 1),
        " & j=", rep(1:(ncategories - 1), clsize), sep = "")
    simdata <- data.frame(y, model.frame(formula = lpformula, data = xdata),
        id, time)
    list(Ysim = Ysim, simdata = simdata, rlatent = err)
}
