test_that("norta",{
  set.seed(1)
  R <- 100
  LatentCorrelation <- toeplitz(c(1, rep(0.8, 2)))
  CommonMarginals <- rep('qlogis', 3)
  SimLogistic <- rnorta(R = R, cor.matrix = LatentCorrelation, distr = CommonMarginals)
  set.seed(1)
  SimNormal <- rsmvnorm(R = R, cor.matrix = LatentCorrelation)
  raw_code <- qlogis(pnorm(SimNormal))
  expect_equal(SimLogistic, raw_code)
})


test_that("rsmvnorm",{
  set.seed(1)
  R <- 100
  cor.matrix <- toeplitz(c(1, 0.4))
  SimBivariateNormal <- rsmvnorm(R = R, cor.matrix = cor.matrix)
  set.seed(1)
  p <- ncol(cor.matrix)
  raw_code <- matrix(rnorm(R * p), R, p) %*% chol(cor.matrix)
  expect_equal(SimBivariateNormal, raw_code)
})
