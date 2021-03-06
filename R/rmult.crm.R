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
#' continuation-ratio model. Options include \code{'probit'}, \code{'logit'},
#' \code{'cloglog'} or \code{'cauchit'}. Required when \code{rlatent = NULL}.
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
#' responses, \code{\link{rmult.clm}} and \code{\link{rmult.acl}} for simulating
#' correlated ordinal responses and \code{\link{rbin}} for simulating
#' correlated binary responses.
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
#' sample_size <- 500
#' cluster_size <- 4
#' beta_intercepts <- c(-1.5, -0.5, 0.5, 1.5)
#' beta_coefficients <- 1
#' x <- rnorm(sample_size * cluster_size)
#' categories_no <- 5
#' identity_matrix <- diag(1, (categories_no - 1) * cluster_size)
#' equicorrelation_matrix <- toeplitz(c(0, rep(0.24, categories_no - 2)))
#' ones_matrix <- matrix(1, cluster_size, cluster_size)
#' latent_correlation_matrix <- identity_matrix +
#'   kronecker(equicorrelation_matrix, ones_matrix)
#' simulated_ordinal_dataset <- rmult.crm(clsize = cluster_size,
#'   intercepts = beta_intercepts, betas = beta_coefficients, xformula = ~x,
#'   cor.matrix = latent_correlation_matrix, link = "probit")
#' head(simulated_ordinal_dataset$Ysim)
#' @export
rmult.crm <- function(clsize = clsize, intercepts = intercepts, betas = betas, # nolint
                      xformula = formula(xdata), xdata = parent.frame(), link = "logit", # nolintr
                      cor.matrix = cor.matrix, rlatent = NULL) { # nolint
  check_cluster_size(clsize)
  beta_coefficients <- check_betas(betas, clsize)
  linear_predictor_formula <- check_xformula(xformula)
  if (!is.environment(xdata)) xdata <- data.frame(na.omit(xdata))
  linear_predictor <- create_linear_predictor(
    beta_coefficients, clsize, linear_predictor_formula, xdata, "rmult.clm"
  )
  sample_size <- nrow(linear_predictor)
  beta_intercepts <- check_intercepts(
    intercepts, clsize, "rmult.crm", sample_size
  )
  categories_no <- ncol(beta_intercepts) / clsize + 1
  linear_predictor_extended <- t(apply(linear_predictor, 1, function(x)
    rep(x, each = categories_no - 1)))
  simulated_latent_responses <- create_rlatent(
    rlatent, sample_size, link, clsize, cor.matrix, "rmult.crm", categories_no
  )
  simulated_ordinal_responses <- apply_threshold(
    linear_predictor_extended, simulated_latent_responses, clsize, "rmult.crm",
    beta_intercepts, categories_no
  )
  create_output(
    simulated_ordinal_responses, sample_size, clsize,
    simulated_latent_responses, linear_predictor_formula, xdata, "rmult.crm",
    categories_no
  )
}
