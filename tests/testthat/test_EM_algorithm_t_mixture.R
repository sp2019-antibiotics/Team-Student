context("Input and output of EM.algorithm.t.mixture()")
library("t.mix")

test_that("Inputs are correctly checked", {
  expect_error(EM.algorithm.t.mixture())
  expect_error(EM.algorithm.t.mixture(NULL))
  expect_error(EM.algorithm.t.mixture("a"))
  expect_error(EM.algorithm.t.mixture(matrix()))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8), NULL))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8), "a"))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8), matrix()))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),c(1,1)))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8), "a"))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12,  c(1,1,1,1)))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12,  matrix(c(1,1,1,1,1,1,1),7), NULL))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12,  matrix(c(1,1,1,1,1,1,1),7), "a"))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12,  matrix(c(1,1,1,1,1,1,1),7), 0.8,NULL))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12,  matrix(c(1,1,1,1,1,1,1),7), 0.8, "a"))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12,  matrix(c(1,1,1,1,1,1,1),7), 0.8, 0.6, NULL))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12,  matrix(c(1,1,1,1,1,1,1),7), 0.8, 0.6, 1000, "a"))
  expect_error(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12,  matrix(c(1,1,1,1,1,1,1),7), 0.8, 0.6, 1000, "Nelder-Mead", NULL))
})

test_that("Output has correct format", {
  expect_true(dim(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12, matrix(c(1,1,1,1,1,1,1),7), 0.8, 0.6, 1000, "Nelder-Mead", 0.001)$parameters)[2]==4)
  expect_true(dim(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12, matrix(c(1,1,1,1,1,1,1),7), 0.8, 0.6, 1000, "Nelder-Mead", 0.001)$parameters)[1]==dim(matrix(c(1,1,1,1,1,1,1),7))[1])
  expect_true(is.numeric(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12, matrix(c(1,1,1,1,1,1,1),7), 0.8, 0.6, 1000, "Nelder-Mead", 0.001)$parameters))
  expect_true(is.numeric(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12, matrix(c(1,1,1,1,1,1,1),7), 0.8, 0.6, 1000, "Nelder-Mead", 0.001)$log_likelihood))
  expect_true(is.numeric(EM.algorithm.t.mixture(c(122,0,44,999,0,123,8),6:12, matrix(c(1,1,1,1,1,1,1),7), 0.8, 0.6, 1000, "Nelder-Mead", 0.001)$iterations))
})


