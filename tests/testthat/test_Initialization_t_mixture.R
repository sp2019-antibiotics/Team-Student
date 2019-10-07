context("Input and output of Initialization.t.mixture()")
library("t.mix")

test_that("Inputs are correctly checked", {
  expect_error(Initialization.t.mixture())
  expect_error(Initialization.t.mixture(NULL))
  expect_error(Initialization.t.mixture("a"))
  expect_error(Initialization.t.mixture(matrix()))
  expect_error(Initialization.t.mixture(matrix(-1)))
  expect_error(Initialization.t.mixture(c(-1,0)))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8),NULL))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8),"a"))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8),c(-2,-2)))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8),2))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8), 6:12, NULL))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8), 6:12, "a"))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8), 6:12, -4))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8), 6:12, 1, "a"))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8), 6:12, 1, NULL))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8), 6:12, 1, -1))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8),6:12, 1, 2, NULL))
  expect_error(Initialization.t.mixture(c(122,0,44,999,0,123,8), 6:12,1, 2, "a"))

})

test_that("Output has correct format", {
  expect_true(is.numeric(Initialization.t.mixture(c(122,4,44,999,6,123,8), 6:12,1, 2, F)))
  expect_true(dim(Initialization.t.mixture(c(122,4,44,999,6,123,8), 6:12,3, 2, F))[2]==3)
  expect_true(dim(Initialization.t.mixture(c(122,4,44,999,6,123,8), 6:12,3, 2, F))[1]==length(c(122,4,44,999,6,123,8)))
  expect_true(dim(Initialization.t.mixture(c(122,0,44,999,6,123,0), 6:12,3, 2, F))[1]==length(c(122,44,999,6,123)))
})



