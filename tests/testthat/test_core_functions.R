test_that("rbin constant beta_intercepts", {
  set.seed(1)
  sample_size <- 10
  cluster_size <- 4
  beta_intercepts <- 0
  beta_coefficients <- 0.2
  latent_correlation_matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
  x <- rep(rnorm(sample_size), each = cluster_size)
  simulated_binary_responses <-
    rbin(clsize = cluster_size, intercepts = beta_intercepts,
         betas = beta_coefficients, xformula = ~ x,
         cor.matrix = latent_correlation_matrix, link = "probit")
  y_sim <-
    as.numeric(c(t(simulated_binary_responses$rlatent)) <=
                 beta_intercepts + beta_coefficients * x)
  expect_equal(c(t(simulated_binary_responses$Ysim)), y_sim)
})

test_that("rbin varying beta_intercepts", {
  set.seed(1)
  sample_size <- 10
  cluster_size <- 4
  beta_intercepts <- c(0, 0.1, 0.2, 0)
  beta_coefficients <- 0.2
  latent_correlation_matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
  x <- rep(rnorm(sample_size), each = cluster_size)
  simulated_binary_responses <-
    rbin(clsize = cluster_size, intercepts = beta_intercepts,
         betas = beta_coefficients, xformula = ~ x,
         cor.matrix = latent_correlation_matrix, link = "probit")
  beta_intercepts <- rep(beta_intercepts, sample_size)
  y_sim <-
    as.numeric(c(t(simulated_binary_responses$rlatent)) <=
                 beta_intercepts + beta_coefficients * x)
  expect_equal(c(t(simulated_binary_responses$Ysim)), y_sim)
})


test_that("rmult.bcl constant beta_coefficients", {
  beta_coefficients <- c(1, 3, 2, 1.25, 3.25, 1.75, 0.75, 2.75, 2.25, 0, 0, 0)
  sample_size <- 10
  categories_no <- 4
  cluster_size <- 3
  set.seed(1)
  x1 <- rep(rnorm(sample_size), each = cluster_size)
  x2 <- rnorm(sample_size * cluster_size)
  xdata <- data.frame(x1, x2)
  latent_correlation_matrix <- kronecker(
    toeplitz(c(1, rep(0.95, cluster_size - 1))), diag(categories_no)
    )
  simulated_nominal_responses <-
    rmult.bcl(clsize = cluster_size, ncategories = categories_no,
              betas = beta_coefficients, xformula = ~ x1 + x2, xdata = xdata,
              cor.matrix = latent_correlation_matrix)
  xmat <- model.matrix(~x1 + x2, data = xdata)
  xmat <- apply(xmat, 2, function(x) rep(x, each = categories_no))
  lin_pred <- matrix(beta_coefficients, nrow = nrow(xmat), ncol = ncol(xmat),
                     byrow = TRUE) * xmat
  lin_pred <- rowSums(lin_pred) + c(t(simulated_nominal_responses$rlatent))
  lin_pred <- matrix(lin_pred, sample_size * cluster_size, categories_no, TRUE)
  y_sim <- apply(lin_pred, 1, which.max)
  expect_equal(c(t(simulated_nominal_responses$Ysim)), y_sim)
})


test_that("rmult.bcl varying beta_coefficients", {
  sample_size <- 10
  categories_no <- 4
  cluster_size <- 3
  beta_coefficients <- matrix(
    c(1, 3, 2, 1.25, 3.25, 1.75, 0.75, 2.75, 2.25, 0, 0, 0,
      2, 1, 2, 1.15, 3.75, 1.25, 0.45, 2.65, 2.85, 0, 0, 0,
      1, 2, 2, 1.75, 3.15, 1.35, 0.55, 2.75, 2.95, 0, 0, 0),
    nrow = cluster_size, byrow = TRUE)
  set.seed(1)
  x1 <- rep(rnorm(sample_size), each = cluster_size)
  x2 <- rnorm(sample_size * cluster_size)
  xdata <- data.frame(x1, x2)
  latent_correlation_matrix <- kronecker(
    toeplitz(c(1, rep(0.95, cluster_size - 1))), diag(categories_no))
  simulated_nominal_responses <- rmult.bcl(clsize = cluster_size,
                                           ncategories = categories_no,
                                           betas = beta_coefficients,
                                           xformula = ~ x1 + x2, xdata = xdata,
                                           cor.matrix =
                                             latent_correlation_matrix)
  xmat <- model.matrix(~x1 + x2, data = xdata)
  xmat <- apply(xmat, 2, function(x) rep(x, each = categories_no))
  lin_pred <- matrix(c(t(beta_coefficients)), nrow = nrow(xmat),
                     ncol = ncol(xmat), byrow = TRUE) * xmat
  lin_pred <- rowSums(lin_pred) + c(t(simulated_nominal_responses$rlatent))
  lin_pred <- matrix(lin_pred, sample_size * cluster_size, categories_no, TRUE)
  y_sim <- apply(lin_pred, 1, which.max)
  expect_equal(c(t(simulated_nominal_responses$Ysim)), y_sim)
})



test_that("rmult.clm varying beta_coefficients", {
  set.seed(1)
  sample_size <- 10
  cluster_size <- 4
  beta_intercepts <- c(-1.5, -0.5, 0.5, 1.5)
  beta_coefficients <- matrix(c(1, 2, 3, 4), 4, 1)
  x <- rep(rnorm(sample_size), each = cluster_size)
  latent_correlation_matrix <- toeplitz(c(1, 0.85, 0.5, 0.15))
  simulated_ordinal_responses <-
    rmult.clm(clsize = cluster_size, intercepts = beta_intercepts,
              betas = beta_coefficients, xformula = ~ x,
              cor.matrix = latent_correlation_matrix, link = "probit")
  u_sim <-  c(t(simulated_ordinal_responses$rlatent)) - c(beta_coefficients) * x
  beta_intercepts <- c(-Inf, beta_intercepts, Inf)
  y_sim <- cut(u_sim, beta_intercepts, labels = FALSE)
  expect_equal(c(t(simulated_ordinal_responses$Ysim)), y_sim)
})



test_that("rmult.clm constant beta_coefficients", {
  set.seed(1)
  sample_size <- 10
  cluster_size <- 4
  beta_intercepts <- c(-1.5, -0.5, 0.5, 1.5)
  beta_coefficients <- 1
  x <- rep(rnorm(sample_size), each = cluster_size)
  latent_correlation_matrix <- toeplitz(c(1, 0.85, 0.5, 0.15))
  simulated_ordinal_responses <-
    rmult.clm(clsize = cluster_size, intercepts = beta_intercepts,
              betas = beta_coefficients, xformula = ~ x,
              cor.matrix = latent_correlation_matrix, link = "probit")
  u_sim <-  c(t(simulated_ordinal_responses$rlatent)) - c(beta_coefficients) * x
  beta_intercepts <- c(-Inf, beta_intercepts, Inf)
  y_sim <- cut(u_sim, beta_intercepts, labels = FALSE)
  expect_equal(c(t(simulated_ordinal_responses$Ysim)), y_sim)
})


test_that("rmult.acl constant beta_coefficients", {
  beta_intercepts <- c(3, 1, 2)
  beta_coefficients <- c(1, 1)
  sample_size <- 10
  cluster_size <- 3
  set.seed(321)
  x1 <- rep(rnorm(sample_size), each = cluster_size)
  x2 <- rnorm(sample_size * cluster_size)
  xdata <- data.frame(x1, x2)
  set.seed(1)
  latent_correlation_matrix <-
    kronecker(toeplitz(c(1, rep(0.95, cluster_size - 1))), diag(4))
  simulated_ordinal_responses <-
    rmult.acl(clsize = cluster_size, intercepts = beta_intercepts,
              betas = beta_coefficients, xformula = ~ x1 + x2, xdata = xdata,
              cor.matrix = latent_correlation_matrix)
  beta_intercepts <- rev(cumsum(rev(c(beta_intercepts, 0))))
  beta_coefficients_bcl <-
    c(beta_intercepts[1], 3, 3,
      beta_intercepts[2], 2, 2,
      beta_intercepts[3], 1, 1,
      beta_intercepts[4], 0, 0)
  set.seed(1)
  simulated_nominal_responses <-
    rmult.bcl(clsize = cluster_size, ncategories = 4,
              betas = beta_coefficients_bcl, xformula = ~ x1 + x2,
              xdata = xdata, cor.matrix = latent_correlation_matrix)
  expect_equal(c(t(simulated_ordinal_responses$Ysim)),
               c(t(simulated_nominal_responses$Ysim)))
})


test_that("rmult.crm constant beta_coefficients", {
  sample_size <- 10
  cluster_size <- 4
  beta_intercepts <- c(-1.5, -0.5, 0.5, 1.5)
  beta_coefficients <- 1
  x <- rnorm(sample_size * cluster_size)
  categories_no <- 5
  latent_correlation_matrix <-
    diag(1, (categories_no - 1) * cluster_size) +
    kronecker(toeplitz(c(0, rep(0.24, categories_no - 2))),
              matrix(1, cluster_size, cluster_size))
  simulated_ordinal_responses <-
    rmult.crm(clsize = cluster_size, intercepts = beta_intercepts,
              betas = beta_coefficients, xformula = ~ x,
              cor.matrix = latent_correlation_matrix, link = "probit")
  u_sim <- c(t(simulated_ordinal_responses$rlatent)) -
    rep(x, each = categories_no - 1)
  y_sim <-
    matrix(as.numeric(t(u_sim <= beta_intercepts)), sample_size * cluster_size,
           categories_no - 1, TRUE)
  for (i in 1:(categories_no - 1))
    y_sim[, i] <- ifelse(y_sim[, i] == 1, i, categories_no)
  y_sim <- apply(y_sim, 1, min)
  expect_equal(c(t(simulated_ordinal_responses$Ysim)), y_sim)
})
