library(testthat)
library(middleDownR)

test_that("find_matches_with_intensity identifies matches", {
  ions_df <- data.frame(`b_ions_charge 1` = rep(100, 7))
  scan_data <- data.frame(mz = c(100, 101), intensity = c(50, 40))
  matches <- find_matches_with_intensity(ions_df, scan_data, ppm_tolerance = 20, sequence = "ACDEFGHI")
  expect_equal(length(matches), 7)
})
