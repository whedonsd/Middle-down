library(testthat)
library(middleDownR)

test_that("load_config reads YAML", {
  tmp <- tempfile(fileext = ".yml")
  writeLines("max_mods: 3", tmp)
  cfg <- load_config(tmp)
  expect_equal(cfg$max_mods, 3)
})
