test_that("rbin constant intercepts",{
  set.seed(1)
  N <- 10
  clsize <- 4
  intercepts <- 0
  betas <- 0.2
  cor.matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
  x <- rep(rnorm(N), each = clsize)
  CorBinRes <- rbin(clsize = clsize, intercepts = intercepts, betas = betas,
                    xformula = ~x, cor.matrix = cor.matrix, link = 'probit')
  Ysim <- as.numeric(c(t(CorBinRes$rlatent)) <= intercepts + betas * x)
  expect_equal(c(t(CorBinRes$Ysim)), Ysim)
})

test_that("rbin varying intercepts",{
  set.seed(1)
  N <- 10
  clsize <- 4
  intercepts <- c(0, 0.1, 0.2, 0)
  betas <- 0.2
  cor.matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
  x <- rep(rnorm(N), each = clsize)
  CorBinRes <- rbin(clsize = clsize, intercepts = intercepts, betas = betas,
                    xformula = ~x, cor.matrix = cor.matrix, link = 'probit')
  intercepts <- rep(intercepts, N)
  Ysim <- as.numeric(c(t(CorBinRes$rlatent)) <= intercepts + betas * x)
  expect_equal(c(t(CorBinRes$Ysim)), Ysim)
})


