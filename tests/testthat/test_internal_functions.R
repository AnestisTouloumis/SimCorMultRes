check_cluster_size

test_that("checking cluster size", {
  expect_silent(check_cluster_size(2))
  expect_error(check_cluster_size(1),
               "'clsize' must be a positive integer greater than or equal to two") # nolintr
  expect_error(check_cluster_size(2.5),
               "'clsize' must be a positive integer greater than or equal to two") # nolintr
})

test_that("checking categories", {
  expect_silent(check_ncategories(3))
  expect_error(check_ncategories(2),
               "'ncategories' must be numeric greater or equal to three")
  expect_error(check_ncategories(5.5),
               "'ncategories' must be a positive integer")
})


test_that("creating marginal distributions", {
  expect_equal(create_distribution("probit"), "qnorm")
  expect_equal(create_distribution("logit"), "qlogis")
  expect_equal(create_distribution("cloglog"), "qgumbel")
  expect_equal(create_distribution("cauchit"), "qcauchy")
})


test_that("checking correlation matrix", {
  cluster_size <- 5
  categories_no <- 3
  correlation_matrix_1 <- matrix(1, cluster_size, cluster_size)
  correlation_matrix_2 <- matrix(1,  cluster_size * categories_no,
                                 cluster_size * categories_no)
  correlation_matrix_3 <- matrix(1,  cluster_size * (categories_no - 1),
                                 cluster_size * (categories_no - 1))
  expect_error(check_correlation_matrix(correlation_matrix_1, cluster_size,
                                        rfctn = "rbin"),
               "'cor_matrix' must be positive definite")
  expect_error(check_correlation_matrix(correlation_matrix_1, cluster_size,
                                        rfctn = "rmult.clm"),
               "'cor_matrix' must be positive definite")
  expect_error(check_correlation_matrix(correlation_matrix_2, cluster_size,
                                        rfctn = "rmult.bcl", categories_no),
               "'cor_matrix' must be positive definite")
  expect_error(check_correlation_matrix(correlation_matrix_1, cluster_size,
                                        rfctn = "rmult.bcl", categories_no),
               "'cor_matrix' must be a 15x15 matrix")
  expect_error(check_correlation_matrix(correlation_matrix_3, cluster_size,
                                        rfctn = "rmult.crm", categories_no),
               "'cor_matrix' must be positive definite")
  expect_error(check_correlation_matrix(correlation_matrix_2, cluster_size,
                                        rfctn = "rmult.crm", categories_no),
               "'cor_matrix' must be a 10x10 matrix")
  correlation_matrix_4 <- diag(1, cluster_size)
  correlation_matrix_5 <- diag(1,  cluster_size * categories_no)
  correlation_matrix_6 <- diag(1,  cluster_size * (categories_no - 1))
  expect_equal(check_correlation_matrix(correlation_matrix_4, cluster_size,
                                        rfctn = "rbin"),
               correlation_matrix_4)
  expect_equal(check_correlation_matrix(correlation_matrix_4, cluster_size,
                                        rfctn = "rmult.clm"),
               correlation_matrix_4)
  expect_equal(check_correlation_matrix(correlation_matrix_5, cluster_size,
                                        rfctn = "rmult.bcl", categories_no),
               correlation_matrix_5)
  expect_equal(check_correlation_matrix(correlation_matrix_6, cluster_size,
                                        rfctn = "rmult.crm", categories_no),
               correlation_matrix_6)
})

test_that("checking covariates formula", {
  xformula <- ~ -1
  expect_error(check_xformula(xformula),
               "No covariates were found in 'formula'")
  xformula <-  ~ x1 + x2
  expect_equal(check_xformula(xformula), as.formula(~ x1 + x2))
  xformula <- y ~ x1 + x2
  expect_equal(check_xformula(xformula), as.formula(y ~ x1 + x2))
  xformula <- y ~ x1 + x2 - 1
  expect_equal(check_xformula(xformula), as.formula(y ~ x1 + x2))
})

test_that("checking intercepts", {
  cluster_size <- 5
  rfctn <- "rbin"
  intercepts <- c(-Inf, 3, Inf)
  expect_error(check_intercepts(intercepts, cluster_size, rfctn),
               "'intercepts' must not be Inf or -Inf")
  intercepts <- c(2, 3)
  expect_error(check_intercepts(intercepts, cluster_size, rfctn),
               "'intercepts' must have either one or 5 elements")
  intercepts <- 2
  ans <- rep(intercepts, cluster_size)
  ans <- cbind(-Inf, ans, Inf)
  colnames(ans) <- c("","intercepts","")
  expect_equal(check_intercepts(intercepts, cluster_size, rfctn), ans)
  intercepts <- 1:cluster_size
  ans <- cbind(-Inf, intercepts, Inf)
  colnames(ans) <- c("","intercepts","")
  expect_equal(check_intercepts(intercepts, cluster_size, rfctn), ans)
  rfctn <- "rmult.clm"
  intercepts <- c(-Inf, 3, Inf)
  expect_error(check_intercepts(intercepts, cluster_size, rfctn),
               "'intercepts' must not be Inf or -Inf")
  intercepts <- c(2, 3)
  ans <- matrix(intercepts, nrow = cluster_size, ncol = 2, byrow = TRUE)
  ans <- cbind(-Inf, ans, Inf)
  colnames(ans) <- NULL
  expect_equal(check_intercepts(intercepts, cluster_size, rfctn), ans)
})
