test_that("norta", {
  set.seed(1)
  R <- 100 # nolint
  latentcorrelation <- toeplitz(c(1, rep(0.8, 2)))
  commonmarginals <- rep("qlogis", 3)
  simlogistic <- rnorta(R = R, cor.matrix = latentcorrelation,
                        distr = commonmarginals)
  set.seed(1)
  simnormal <- rsmvnorm(R = R, cor.matrix = latentcorrelation)
  raw_code <- qlogis(pnorm(simnormal))
  expect_equal(simlogistic, raw_code)
})


test_that("rsmvnorm", {
  set.seed(1)
  R <- 100 # nolint
  cor.matrix <- toeplitz(c(1, 0.4)) # nolint
  simbivariatenormal <- rsmvnorm(R = R, cor.matrix = cor.matrix)
  set.seed(1)
  p <- ncol(cor.matrix)
  raw_code <- matrix(rnorm(R * p), R, p) %*% chol(cor.matrix)
  expect_equal(simbivariatenormal, raw_code)
})
