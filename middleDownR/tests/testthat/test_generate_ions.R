library(testthat)
library(middleDownR)

test_that("generate_ions_for_sequence returns expected elements", {
  mods <- data.frame(position = integer(), variable = double())
  ions <- generate_ions_for_sequence("ACDE", mods)
  expect_true("b_ions_charge 1" %in% names(ions))
  expect_equal(length(ions[["b_ions_charge 1"]]), nchar("ACDE") - 1)
})
