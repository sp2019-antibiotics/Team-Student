context("Input and output of weighted.quantile()")
library("t.mix")

test_that("Inputs are correctly checked", {
  expect_error(weighted.quantile())
  expect_error(weighted.quantile(NULL))
  expect_error(weighted.quantile("a"))
  expect_error(weighted.quantile(matrix()))
  expect_error(weighted.quantile(6:12, NULL))
  expect_error(weighted.quantile(6:12, "a"))
  expect_error(weighted.quantile(6:12, c(-1,-1)))
  expect_error(weighted.quantile(6:12, c(122,0,44,999,0,123,8), NULL))
  expect_error(weighted.quantile(6:12, c(122,0,44,999,0,123,8), "a"))
  expect_error(weighted.quantile(6:12, c(122,0,44,999,0,123,8), 2))
  expect_error(weighted.quantile(6:12, c(122,0,44,999,0,123,8), c(0.5,2)))
})

test_that("Output has correct format", {
  expect_true(is.numeric(weighted.quantile(6:12, c(122,0,44,999,0,123,8))))
})

