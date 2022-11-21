#' Bivariate NORTA Generated Correlation
#'
#' Simulated dataset to understand
#'
#' @format ## `simulation`
#' A data frame with 100 rows and 4 columns:
#' \describe{
#'   \item{rho}{numeric indicating the value of the correlation parameter.}
#'   \item{normal}{numeric indicating the simulated average of the correlation parameter with
#'   normal margins.}
#'   \item{logistic}{numeric indicating the simulated average of the correlation parameter with
#'   logistic margins.}
#'   \item{gumbel}{numeric indicating the simulated average of the correlation parameter with
#'   gumbel margins.}
#' }
#' @examples
#' simulation |>
#'   plot(rho - normal ~ rho, data = _, type = "l", col = "blue",
#'        ylim = c(0, 0.016),
#'        ylab = "Difference between true and simulated correlation values",
#'        xlab = "Correlation parameter")
#' simulation |>
#'  points(rho - logistic ~ rho, data = _, type = "l", col = "red")
#' simulation |>
#'   points(rho - gumbel ~ rho, data = _, type = "l", col = "grey")
#' legend("topright", legend = c("Normal", "Logistic", "Gumbel"),
#'        col = c("blue", "red", "grey"), pch = "l" )
"simulation"
