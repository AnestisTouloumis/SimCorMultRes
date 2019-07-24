#' Simulating Correlated Ordinal Responses Conditional on a Marginal
#' Adjacent-Category Logit Model Specification
#'
#' Simulates correlated ordinal responses assuming an adjacent-category logit
#' model for the marginal probabilities.
#'
#' The formulae are easier to read from either the Vignette or the Reference
#' Manual (both available
#' \href{https://CRAN.R-project.org/package=SimCorMultRes}{here}).
#'
#' The assumed marginal adjacent-category logit model is \deqn{log
#' \frac{Pr(Y_{it}=j |x_{it})}{Pr(Y_{it}=j+1 |x_{it})}=\beta_{tj0}
#' + \beta^{'}_{t} x_{it}}
#' For subject \eqn{i}, \eqn{Y_{it}} is the \eqn{t}-th ordinal response
#' and \eqn{x_{it}} is the associated covariates vector. Also \eqn{\beta_{tj0}}
#' is the \eqn{j}-th category-specific intercept at the \eqn{t}-th measurement
#' occasion and \eqn{\beta_{t}} is the regression
#' parameter vector at the \eqn{t}-th measurement occasion.
#'
#' The ordinal response \eqn{Y_{it}} is obtained by utilizing the threshold
#' approach described in the Vignette. This approach is based on the connection
#' between baseline-category and adjacent-category logit models.
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
#' vectors denoted by \eqn{e^{O3}_{itj}} in the Vignette. To import
#' simulated values for the latent random vectors without utilizing the NORTA
#' method, the user can employ the \code{rlatent} argument. In this case,
#' row \eqn{i} corresponds to subject \eqn{i} and columns
#' \eqn{(t-1)*\code{ncategories}+1,...,t*\code{ncategories}} should contain the
#' realization of \eqn{e^{O3}_{it1},...,e^{O3}_{itJ}}, respectively, for
#' \eqn{t=1,\ldots,\code{clsize}}.
#'
#'
#' @param clsize integer indicating the common cluster size.
#' @param intercepts numerical vector or matrix containing the intercepts of
#' the marginal adjacent-category logit model.
#' @param betas numerical vector or matrix containing the value of the marginal
#' regression parameter vector.
#' @param xformula formula expression as in other marginal regression models
#' but without including a response variable.
#' @param xdata optional data frame containing the variables provided in
#' \code{xformula}.
#' @param cor.matrix matrix indicating the correlation matrix of the
#' multivariate normal distribution when the NORTA method is employed
#' (\code{rlatent = NULL}).
#' @param rlatent matrix with \code{(clsize * ncategories)} columns containing
#' realizations of the latent random vectors when the NORTA method is not
#' preferred. See details for more info.
#' @return Returns a list that has components: \item{Ysim}{the simulated
#' nominal responses. Element (\eqn{i},\eqn{t}) represents the realization of
#' \eqn{Y_{it}}.} \item{simdata}{a data frame that includes the simulated
#' response variables (y), the covariates specified by \code{xformula},
#' subjects' identities (id) and the corresponding measurement occasions
#' (time).} \item{rlatent}{the latent random variables denoted by
#' \eqn{e^{O3}_{itj}} in the Vignette.}
#' @author Anestis Touloumis
#' @seealso \code{\link{rbin}} for simulating correlated binary responses,
#' \code{\link{rmult.clm}} and \code{\link{rmult.crm}} for simulating
#' correlated ordinal responses, and \code{\link{rmult.bcl}} for simulating
#' nominal responses.
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
#'
#' Touloumis, A., Agresti, A. and Kateri, M. (2013) GEE for multinomial
#' responses using a local odds ratios parameterization. \emph{Biometrics}
#' \bold{69}, 633--640.
#'
#' @examples
#' ## See Example 3.4 in the Vignette.
#' beta_intercepts <- c(3, 1, 2)
#' beta_coefficients <- c(1, 1)
#' sample_size <- 500
#' cluster_size <- 3
#' set.seed(321)
#' x1 <- rep(rnorm(sample_size), each = cluster_size)
#' x2 <- rnorm(sample_size * cluster_size)
#' xdata <- data.frame(x1, x2)
#' identity_matrix <- diag(4)
#' equicorrelation_matrix <- toeplitz(c(1, rep(0.95, cluster_size - 1)))
#' latent_correlation_matrix <- kronecker(
#'   equicorrelation_matrix,
#'   identity_matrix
#' )
#' simulated_ordinal_dataset <- rmult.acl(
#'   clsize = cluster_size,
#'   intercepts = beta_intercepts, betas = beta_coefficients,
#'   xformula = ~ x1 + x2, xdata = xdata,
#'   cor.matrix = latent_correlation_matrix
#' )
#' suppressPackageStartupMessages(library("multgee"))
#' ordinal_gee_model <- ordLORgee(y ~ x1 + x2,
#'   data = simulated_ordinal_dataset$simdata, id = id, repeated = time,
#'   LORstr = "time.exch", link = "acl"
#' )
#' round(coef(ordinal_gee_model), 2)
#' @export
rmult.acl <- function(clsize = clsize, intercepts = intercepts, betas = betas, # nolint
                      xformula = formula(xdata), xdata = parent.frame(), cor.matrix = cor.matrix, # nolint
                      rlatent = NULL) {
  check_cluster_size(clsize)
  beta_intercepts <- check_intercepts(intercepts, clsize, "rmult.acl")
  categories_no <- ncol(beta_intercepts) + 1
  beta_coefficients <- check_betas(betas, clsize)
  betas_bcl <- create_betas_acl2bcl(
    beta_intercepts, categories_no, beta_coefficients
  )
  rmult.bcl(
    clsize, categories_no, betas_bcl, xformula, xdata, cor.matrix, rlatent
  )
}
