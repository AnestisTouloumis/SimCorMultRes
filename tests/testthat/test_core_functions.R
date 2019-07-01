test_that("rbin constant intercepts", {
  set.seed(1)
  N <- 10
  clsize <- 4
  intercepts <- 0
  betas <- 0.2
  cor.matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
  x <- rep(rnorm(N), each = clsize)
  corbinres <- rbin(clsize = clsize, intercepts = intercepts, betas = betas,
                    xformula = ~x, cor.matrix = cor.matrix, link = "probit")
  Ysim <- as.numeric(c(t(corbinres$rlatent)) <= intercepts + betas * x)
  expect_equal(c(t(corbinres$Ysim)), Ysim)
})

test_that("rbin varying intercepts", {
  set.seed(1)
  N <- 10
  clsize <- 4
  intercepts <- c(0, 0.1, 0.2, 0)
  betas <- 0.2
  cor.matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
  x <- rep(rnorm(N), each = clsize)
  corbinres <- rbin(clsize = clsize, intercepts = intercepts, betas = betas,
                    xformula = ~x, cor.matrix = cor.matrix, link = "probit")
  intercepts <- rep(intercepts, N)
  Ysim <- as.numeric(c(t(corbinres$rlatent)) <= intercepts + betas * x)
  expect_equal(c(t(corbinres$Ysim)), Ysim)
})


test_that("rmult.bcl constant betas", {
  betas <- c(1, 3, 2, 1.25, 3.25, 1.75, 0.75, 2.75, 2.25, 0, 0, 0)
  N <- 10
  ncategories <- 4
  clsize <- 3
  set.seed(1)
  x1 <- rep(rnorm(N), each = clsize)
  x2 <- rnorm(N * clsize)
  xdata <- data.frame(x1, x2)
  cor.matrix <- kronecker(toeplitz(c(1, rep(0.95, clsize - 1))),
                          diag(ncategories))
  cornomres <- rmult.bcl(clsize = clsize, ncategories = ncategories,
                         betas = betas, xformula = ~x1 + x2, xdata = xdata,
                         cor.matrix = cor.matrix)
  Xmat <- model.matrix( ~x1 + x2, data = xdata)
  Xmat <- apply(Xmat, 2, function(x) rep(x, each = ncategories))
  lin_pred <- matrix(betas, nrow = nrow(Xmat), ncol = ncol(Xmat),
                     byrow = TRUE) * Xmat
  lin_pred <- rowSums(lin_pred) + c(t(cornomres$rlatent))
  lin_pred <- matrix(lin_pred, N * clsize, ncategories, TRUE)
  Ysim <- apply(lin_pred, 1, which.max)
  expect_equal(c(t(cornomres$Ysim)), Ysim)
})


test_that("rmult.bcl varying betas", {
  N <- 10
  ncategories <- 4
  clsize <- 3
  betas <- matrix(c(1, 3, 2, 1.25, 3.25, 1.75, 0.75, 2.75, 2.25, 0, 0, 0,
                    2, 1, 2, 1.15, 3.75, 1.25, 0.45, 2.65, 2.85, 0, 0, 0,
                    1, 2, 2, 1.75, 3.15, 1.35, 0.55, 2.75, 2.95, 0, 0, 0),
                  nrow = clsize, byrow = TRUE)
  set.seed(1)
  x1 <- rep(rnorm(N), each = clsize)
  x2 <- rnorm(N * clsize)
  xdata <- data.frame(x1, x2)
  cor.matrix <- kronecker(toeplitz(c(1, rep(0.95, clsize - 1))),
                          diag(ncategories))
  cornomres <- rmult.bcl(clsize = clsize, ncategories = ncategories,
                         betas = betas, xformula = ~x1 + x2, xdata = xdata,
                         cor.matrix = cor.matrix)
  Xmat <- model.matrix( ~x1 + x2, data = xdata)
  Xmat <- apply(Xmat, 2, function(x) rep(x, each = ncategories))
  lin_pred <- matrix(c(t(betas)), nrow = nrow(Xmat), ncol = ncol(Xmat),
                     byrow = TRUE) * Xmat
  lin_pred <- rowSums(lin_pred) + c(t(cornomres$rlatent))
  lin_pred <- matrix(lin_pred, N * clsize, ncategories, TRUE)
  Ysim <- apply(lin_pred, 1, which.max)
  expect_equal(c(t(cornomres$Ysim)), Ysim)
})



test_that("rmult.clm varying betas", {
  set.seed(1)
  N <- 10
  clsize <- 4
  intercepts <- c(-1.5, -0.5, 0.5, 1.5)
  betas <- matrix(c(1, 2, 3, 4), 4, 1)
  x <- rep(rnorm(N), each = clsize)
  cor.matrix <- toeplitz(c(1, 0.85, 0.5, 0.15))
  corordres <- rmult.clm(clsize = clsize, intercepts = intercepts,
                         betas = betas, xformula = ~x, cor.matrix = cor.matrix,
                         link = "probit")
  U <-  c(t(corordres$rlatent)) - c(betas) * x
  intercepts <- c(-Inf, intercepts, Inf)
  Ysim <- cut(U, intercepts, labels = FALSE)
  expect_equal(c(t(corordres$Ysim)), Ysim)
})



test_that("rmult.clm constant betas", {
  set.seed(1)
  N <- 10
  clsize <- 4
  intercepts <- c(-1.5, -0.5, 0.5, 1.5)
  betas <- 1
  x <- rep(rnorm(N), each = clsize)
  cor.matrix <- toeplitz(c(1, 0.85, 0.5, 0.15))
  corordres <- rmult.clm(clsize = clsize, intercepts = intercepts,
                         betas = betas, xformula = ~x, cor.matrix = cor.matrix,
                         link = "probit")
  U <-  c(t(corordres$rlatent)) - c(betas) * x
  intercepts <- c(-Inf, intercepts, Inf)
  Ysim <- cut(U, intercepts, labels = FALSE)
  expect_equal(c(t(corordres$Ysim)), Ysim)
})


test_that("rmult.acl constant betas", {
  intercepts <- c(3, 1, 2)
  betas <- c(1, 1)
  N <- 10
  clsize <- 3
  set.seed(321)
  x1 <- rep(rnorm(N), each = clsize)
  x2 <- rnorm(N * clsize)
  xdata <- data.frame(x1, x2)
  set.seed(1)
  cor.matrix <- kronecker(toeplitz(c(1, rep(0.95, clsize - 1))), diag(4))
  corordres <- rmult.acl(clsize = clsize, intercepts = intercepts,
                         betas = betas, xformula = ~ x1 + x2, xdata = xdata,
                         cor.matrix = cor.matrix)
  intercepts <- rev(cumsum(rev(c(intercepts, 0))))
  betas_bcl <- c(intercepts[1], 3, 3, intercepts[2], 2, 2,
                 intercepts[3], 1, 1, intercepts[4], 0, 0)
  set.seed(1)
  cornomres <- rmult.bcl(clsize = clsize, ncategories = 4, betas = betas_bcl,
                         xformula = ~x1 + x2, xdata = xdata,
                         cor.matrix = cor.matrix)
  expect_equal(c(t(corordres$Ysim)), c(t(cornomres$Ysim)))
})


test_that("rmult.crm constant betas", {
  N <- 10
  clsize <- 4
  intercepts <- c(-1.5, -0.5, 0.5, 1.5)
  betas <- 1
  x <- rnorm(N * clsize)
  ncategories <- 5
  cor.matrix <- diag(1, (ncategories - 1) * clsize) +
    kronecker(toeplitz(c(0, rep(0.24, ncategories - 2))),
              matrix(1, clsize, clsize))
  corordres <- rmult.crm(clsize = clsize, intercepts = intercepts,
                         betas = betas, xformula = ~x, cor.matrix = cor.matrix,
                         link = "probit")
  U <- c(t(corordres$rlatent)) - rep(x, each = ncategories - 1)
  Ysim <- matrix(as.numeric(t(U <= intercepts)), N * clsize, ncategories -
                   1, TRUE)
  for (i in 1:(ncategories - 1)) Ysim[, i] <- ifelse(Ysim[, i] ==
                                                       1, i, ncategories)
  Ysim <- apply(Ysim, 1, min)
  expect_equal(c(t(corordres$Ysim)), Ysim)
})
