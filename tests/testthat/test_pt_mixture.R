context("Input and output of pt.mixture()")
library("t.mix")

test_that("Inputs are correctly checked", {
  expect_error(pt.mixture())
  expect_error(pt.mixture(NULL))
  expect_error(pt.mixture("a"))
  expect_error(pt.mixture(matrix()))
  expect_error(pt.mixture(c(1,2,3)))
  expect_error(pt.mixture(c(1,2,3), NULL))
  expect_error(pt.mixture(c(1,2,3), "a"))
  expect_error(pt.mixture(c(1,2,3), c(-1,0,0)))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), NULL))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), "a"))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(1)))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), NULL))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), "a"))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), c(-1,-1,-1)))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), c(1)))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), c(2,2,2), NULL))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), c(2,2,2), "a"))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), c(2,2,2), c(1)))
  expect_error(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), c(2,2,2), c(-1,-1,-1)))
})

test_that("Output has correct format", {
  expect_true(length(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), c(2,2,2), c(10,10,10)))==length(c(1,2,3)))
  expect_true(length(pt.mixture(1,1,4,3,10))==length(1))
  expect_true(is.numeric(pt.mixture(c(1,2,3), c(0.4,0.1,0.5), c(2,3,4), c(2,2,2), c(10,10,10))))
})


