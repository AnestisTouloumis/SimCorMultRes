test_that("checking cluster size", {
  expect_silent(check_cluster_size(2))
  expect_error(
    check_cluster_size(1),
    "'clsize' must be a positive integer greater than or equal to two"
  )
  expect_error(
    check_cluster_size(2.5),
    "'clsize' must be a positive integer greater than or equal to two"
  )
})


test_that("checking categories", {
  expect_silent(check_ncategories(3))
  expect_error(
    check_ncategories(2),
    "'ncategories' must be numeric greater or equal to three"
  )
  expect_error(
    check_ncategories(5.5),
    "'ncategories' must be a positive integer"
  )
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
  correlation_matrix_2 <- matrix(
    1, cluster_size * categories_no, cluster_size * categories_no
  )
  correlation_matrix_3 <- matrix(
    1, cluster_size * (categories_no - 1),
    cluster_size * (categories_no - 1)
  )
  correlation_matrix_4 <- diag(1, cluster_size)
  correlation_matrix_5 <- diag(1, cluster_size * categories_no)
  correlation_matrix_6 <- diag(1, cluster_size * (categories_no - 1))
  expect_equal(
    check_correlation_matrix(
      correlation_matrix_4, cluster_size, "rbin"
    ),
    correlation_matrix_4
  )
  expect_equal(
    check_correlation_matrix(
      correlation_matrix_4, cluster_size, "rmult.clm"
    ),
    correlation_matrix_4
  )
  expect_equal(
    check_correlation_matrix(
      correlation_matrix_5, cluster_size, "rmult.bcl", categories_no
    ),
    correlation_matrix_5
  )
  expect_equal(
    check_correlation_matrix(
      correlation_matrix_6, cluster_size, "rmult.crm", categories_no
    ),
    correlation_matrix_6
  )
  expect_error(
    check_correlation_matrix(
      correlation_matrix_1, cluster_size, "rbin"
    ),
    "'cor_matrix' must be positive definite"
  )
  expect_error(
    check_correlation_matrix(
      correlation_matrix_1, cluster_size, "rmult.clm"
    ),
    "'cor_matrix' must be positive definite"
  )
  expect_error(
    check_correlation_matrix(
      correlation_matrix_2, cluster_size, "rmult.bcl", categories_no
    ),
    "'cor_matrix' must be positive definite"
  )
  expect_error(
    check_correlation_matrix(
      correlation_matrix_1, cluster_size, "rmult.bcl", categories_no
    ),
    "'cor_matrix' must be a 15x15 matrix"
  )
  expect_error(
    check_correlation_matrix(
      correlation_matrix_3, cluster_size, "rmult.crm", categories_no
    ),
    "'cor_matrix' must be positive definite"
  )
  expect_error(
    check_correlation_matrix(
      correlation_matrix_2, cluster_size, "rmult.crm", categories_no
    ),
    "'cor_matrix' must be a 10x10 matrix"
  )
})


test_that("checking covariates formula", {
  xformula <- ~ x1 + x2
  expect_equal(check_xformula(xformula), as.formula(~ x1 + x2))
  xformula <- y ~ x1 + x2
  expect_equal(check_xformula(xformula), as.formula(y ~ x1 + x2))
  xformula <- y ~ x1 + x2 - 1
  expect_equal(check_xformula(xformula), as.formula(y ~ x1 + x2))
  xformula <- ~ -1
  expect_error(
    check_xformula(xformula),
    "No covariates were found in 'formula'"
  )
  xformula <- y ~
  expect_error(
    check_xformula(xformula),
    "No covariates were found in 'formula'"
  )
})


test_that("checking intercepts", {
  cluster_size <- 4
  rfctn <- "rbin"
  intercepts <- 2
  ans <- rep(intercepts, cluster_size)
  ans <- cbind(-Inf, ans, Inf)
  colnames(ans) <- c("", "intercepts", "")
  expect_equal(check_intercepts(intercepts, cluster_size, rfctn), ans)
  intercepts <- 1:cluster_size
  ans <- cbind(-Inf, intercepts, Inf)
  colnames(ans) <- c("", "intercepts", "")
  expect_equal(check_intercepts(intercepts, cluster_size, rfctn), ans)
  intercepts <- c(-Inf, 3, Inf)
  expect_error(
    check_intercepts(intercepts, cluster_size, rfctn),
    "'intercepts' must not be Inf or -Inf"
  )
  intercepts <- c(2, 3)
  expect_error(
    check_intercepts(intercepts, cluster_size, rfctn),
    "'intercepts' must have either one or 4 elements"
  )
  rfctn <- "rmult.clm"
  intercepts <- c(2, 3)
  ans <- matrix(intercepts, nrow = cluster_size, ncol = 2, byrow = TRUE)
  ans <- cbind(-Inf, ans, Inf)
  colnames(ans) <- NULL
  expect_equal(check_intercepts(intercepts, cluster_size, rfctn), ans)
  intercepts <- matrix(c(2, 3), nrow = cluster_size, ncol = 2, byrow = TRUE)
  ans <- cbind(-Inf, intercepts, Inf)
  colnames(ans) <- NULL
  expect_equal(check_intercepts(intercepts, cluster_size, rfctn), ans)
  intercepts <- matrix(c(3, 2, 1, 2, 3, 4, 5, 6), cluster_size, 2, TRUE)
  intercepts <- c(-Inf, 3, Inf)
  expect_error(
    check_intercepts(intercepts, cluster_size, rfctn),
    "'intercepts' must not be Inf or -Inf"
  )
  intercepts <- matrix(c(3, 2, 1, 2, 3, 4, 5, 6), cluster_size, 2, TRUE)
  expect_error(
    check_intercepts(intercepts, cluster_size, rfctn),
    "'intercepts' must be increasing at each row"
  )
  intercepts <- c(4, 2)
  expect_error(
    check_intercepts(intercepts, cluster_size, rfctn),
    "'intercepts' must be increasing"
  )
  intercepts <- matrix(1:8, nrow = 2, ncol = cluster_size, byrow = TRUE)
  expect_error(check_intercepts(intercepts, cluster_size, rfctn))
  rfctn <- "rmult.crm"
  intercepts <- c(-1.5, -0.5, 0.5, 1.5)
  categories_no <- 5
  ans <- matrix(
    rep(intercepts, categories_no), categories_no,
    (categories_no - 1) * cluster_size, TRUE
  )
  colnames(ans) <- NULL
  expect_equal(
    check_intercepts(intercepts, cluster_size, rfctn, categories_no),
    ans
  )
  intercepts <- matrix(
    runif(
      cluster_size * (categories_no - 1), 0, 1
    ),
    cluster_size, categories_no - 1
  )
  ans <- matrix(
    c(t(intercepts)), categories_no,
    (categories_no - 1) * cluster_size, TRUE
  )
  colnames(ans) <- NULL
  expect_equal(
    check_intercepts(intercepts, cluster_size, rfctn, categories_no), ans
  )
})


test_that("checking betas", {
  cluster_size <- 2
  betas <- c(1, 2)
  expect_equal(check_betas(betas, cluster_size), rep(betas, cluster_size))
  betas <- matrix(1:4, cluster_size, 2, TRUE)
  expect_equal(check_betas(betas, cluster_size), c(1, 2, 3, 4))
  betas <- matrix(1:6, 3, cluster_size, TRUE)
  expect_error(
    check_betas(betas, cluster_size),
    "The number of rows in 'betas' should be equal to 'clsize'"
  )
})
