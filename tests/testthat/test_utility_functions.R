test_that("norta", {
  set.seed(1)
  sample_size <- 100
  latentcorrelation <- toeplitz(c(1, rep(0.8, 2)))
  commonmarginals <- rep("qlogis", 3)
  simlogistic <- rnorta(R = sample_size, cor.matrix = latentcorrelation,
                        distr = commonmarginals)
  set.seed(1)
  simnormal <- rsmvnorm(R = sample_size, cor.matrix = latentcorrelation)
  raw_code <- qlogis(pnorm(simnormal))
  expect_equal(simlogistic, raw_code)
})


test_that("rsmvnorm", {
  set.seed(1)
  sample_size <- 100 # nolint
  cor_matrix <- toeplitz(c(1, 0.4)) # nolint
  simbivariatenormal <- rsmvnorm(R = sample_size, cor.matrix = cor_matrix)
  set.seed(1)
  p <- ncol(cor_matrix)
  raw_code <- matrix(rnorm(sample_size * p), sample_size, p) %*%
    chol(cor_matrix)
  expect_equal(simbivariatenormal, raw_code)
})
