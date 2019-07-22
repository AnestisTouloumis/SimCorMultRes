test_that("rnorta", {
  set.seed(1)
  sample_size <- 100
  latent_correlation_matrix <- toeplitz(c(1, rep(0.8, 2)))
  common_marginals <- rep("qlogis", 3)
  sim_logistic <- rnorta(R = sample_size,
                         cor.matrix = latent_correlation_matrix,
                         distr = common_marginals)
  set.seed(1)
  sim_normal <- rsmvnorm(R = sample_size,
                         cor.matrix = latent_correlation_matrix)
  raw_code <- qlogis(pnorm(sim_normal))
  expect_equal(sim_logistic, raw_code)
})


test_that("rsmvnorm", {
  set.seed(1)
  sample_size <- 100 # nolint
  correlation_matrix <- toeplitz(c(1, 0.4)) # nolint
  sim_bivariate_normal <- rsmvnorm(R = sample_size,
                                 cor.matrix = correlation_matrix)
  set.seed(1)
  p <- ncol(correlation_matrix)
  raw_code <- matrix(rnorm(sample_size * p), sample_size, p) %*%
    chol(correlation_matrix)
  expect_equal(sim_bivariate_normal, raw_code)
})


test_that("rnorta sample size", {
  R <- 0
  expect_error(if (all.equal(R, as.integer(R)) != TRUE | R < 1)
    stop("'R' must be a positive integer"))
  R <- 3.4
  expect_error(
    if (all.equal(R, as.integer(R)) != TRUE | R < 1)
      stop("'R' must be a positive integer"))
  R <- -3
  expect_error(
    if (all.equal(R, as.integer(R)) != TRUE | R < 1)
      stop("'R' must be a positive integer"))
          })
