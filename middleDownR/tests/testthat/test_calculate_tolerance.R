library(testthat)
library(middleDownR)

test_that("calculate_tolerance computes ppm", {
  expect_equal(calculate_tolerance(1000, 10), 0.01)
})
