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


test_that("rmult.bcl constant betas",{
  betas <- c(1, 3, 2, 1.25, 3.25, 1.75, 0.75, 2.75, 2.25, 0, 0, 0)
  N <- 10
  ncategories <- 4
  clsize <- 3
  set.seed(1)
  x1 <- rep(rnorm(N), each = clsize)
  x2 <- rnorm(N * clsize)
  xdata <- data.frame(x1, x2)
  cor.matrix <- kronecker(toeplitz(c(1, rep(0.95, clsize - 1))), diag(ncategories))
  CorNomRes <- rmult.bcl(clsize = clsize, ncategories = ncategories, betas = betas,
                         xformula = ~x1 + x2, xdata = xdata, cor.matrix = cor.matrix)
  Xmat <- model.matrix( ~x1 + x2, data = xdata)
  Xmat <- apply(Xmat, 2, function(x) rep(x, each = ncategories))
  lin_pred <- matrix(betas, nrow = nrow(Xmat), ncol = ncol(Xmat),
                     byrow = TRUE) * Xmat
  lin_pred <- rowSums(lin_pred) + c(t(CorNomRes$rlatent))
  lin_pred <- matrix(lin_pred, N * clsize, ncategories, TRUE)
  Ysim <- apply(lin_pred, 1, which.max)
  expect_equal(c(t(CorNomRes$Ysim)), Ysim)
})


test_that("rmult.bcl varying betas",{
  N <- 10
  ncategories <- 4
  clsize <- 3
  betas <- matrix(c(1, 3, 2, 1.25, 3.25, 1.75, 0.75, 2.75, 2.25, 0, 0, 0,
                    2, 1, 2, 1.15, 3.75, 1.25, 0.45, 2.65, 2.85, 0, 0, 0,
                    1, 2, 2, 1.75, 3.15, 1.35, 0.55, 2.75, 2.95, 0, 0, 0),
                  nrow = clsize, byrow=TRUE)
  set.seed(1)
  x1 <- rep(rnorm(N), each = clsize)
  x2 <- rnorm(N * clsize)
  xdata <- data.frame(x1, x2)
  cor.matrix <- kronecker(toeplitz(c(1, rep(0.95, clsize - 1))), diag(ncategories))
  CorNomRes <- rmult.bcl(clsize = clsize, ncategories = ncategories, betas = betas,
                         xformula = ~x1 + x2, xdata = xdata, cor.matrix = cor.matrix)
  Xmat <- model.matrix( ~x1 + x2, data = xdata)
  Xmat <- apply(Xmat, 2, function(x) rep(x, each = ncategories))
  lin_pred <- matrix(c(t(betas)), nrow = nrow(Xmat), ncol = ncol(Xmat),
                     byrow = TRUE) * Xmat
  lin_pred <- rowSums(lin_pred) + c(t(CorNomRes$rlatent))
  lin_pred <- matrix(lin_pred, N * clsize, ncategories, TRUE)
  Ysim <- apply(lin_pred, 1, which.max)
  expect_equal(c(t(CorNomRes$Ysim)), Ysim)
})
