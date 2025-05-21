context("assemble_ions_df")

# Minimal version of assemble_ions_df taken from the notebook
assemble_ions_df <- function(df, target_scan_number, target_PID) {
  abc_types <- c("a", "b", "c", "a_nl_O", "b_nl_O", "c_nl_O",
                 "a_nl_N", "b_nl_N", "c_nl_N", "a_nl_Kme3", "b_nl_Kme3",
                 "c_nl_Kme3", "a_nl_pST", "b_nl_pST", "c_nl_pST",
                 "a_nl_Ksu", "b_nl_Ksu", "c_nl_Ksu")
  yz_types <- c("y", "z", "y_nl_O", "z_nl_O", "y_nl_N", "z_nl_N",
                "y_nl_Kme3", "z_nl_Kme3", "y_nl_pST", "z_nl_pST",
                "y_nl_Ksu", "z_nl_Ksu")

  a_b_c_ions_df <- df[df$scan_number == as.numeric(target_scan_number) &
                        df$PID == as.numeric(target_PID) &
                        df$ion_type %in% abc_types, ]

  y_z_ions_df <- df[df$scan_number == as.numeric(target_scan_number) &
                      df$PID == as.numeric(target_PID) &
                      df$ion_type %in% yz_types, ]

  y_z_ions_df$ion_position <- nchar(y_z_ions_df$peptide) - y_z_ions_df$ion_position + 1

  list(a_b_c_ions_df = a_b_c_ions_df, y_z_ions_df = y_z_ions_df)
}


test_that("assemble_ions_df splits ions and converts positions", {
  df <- data.frame(
    scan_number = rep(1, 5),
    PID = rep(1, 5),
    ion_type = c("a", "b", "c", "y", "z"),
    ion_position = 1:5,
    peptide = rep("PEPTIDE", 5),
    stringsAsFactors = FALSE
  )

  res <- assemble_ions_df(df, 1, 1)
  abc <- res$a_b_c_ions_df
  yz <- res$y_z_ions_df

  expect_equal(sort(abc$ion_type), c("a", "b", "c"))
  expect_equal(abc$ion_position, 1:3)

  expect_equal(yz$ion_type, c("y", "z"))
  expect_equal(yz$ion_position, c(4, 3))
})
