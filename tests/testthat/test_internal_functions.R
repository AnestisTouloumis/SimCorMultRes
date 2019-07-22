check_cluster_size

test_that("checking cluster size", {
  expect_silent(check_cluster_size(2))
  expect_error(check_cluster_size(1))
  expect_error(check_cluster_size(2.5))
})

test_that("checking categories", {
  expect_silent(check_ncategories(3))
  expect_error(check_ncategories(2))
  expect_error(check_ncategories(2.5))
})


test_that("creating marginal distributions", {
  expect_equal(create_distribution("probit"), "qnorm")
  expect_equal(create_distribution("logit"), "qlogis" )
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
                                        rfctn = "rbin"))
  expect_error(check_correlation_matrix(correlation_matrix_1, cluster_size,
                                        rfctn = "rmult.clm"))
  expect_error(check_correlation_matrix(correlation_matrix_2, cluster_size,
                                        rfctn = "rmult.bcl", categories_no))
  expect_error(check_correlation_matrix(correlation_matrix_3, cluster_size,
                                        rfctn = "rmult.crm", categories_no))
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
