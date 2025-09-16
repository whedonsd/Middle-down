#' Assemble ion data frames for a given scan/PID combination
#'
#' @param df Data frame containing ion observations.
#' @param target_scan_number Scan number of interest.
#' @param target_PID Peptide identifier of interest.
#'
#' @return List containing `a_b_c_ions_df` and `y_z_ions_df` data frames.
#' @export
assemble_ions_df <- function(df, target_scan_number, target_PID) {
  scan_num <- as.numeric(target_scan_number)
  pid_num <- as.numeric(target_PID)

  abc_types <- c(
    "a", "b", "c",
    "a_nl_O", "b_nl_O", "c_nl_O",
    "a_nl_N", "b_nl_N", "c_nl_N",
    "a_nl_Kme3", "b_nl_Kme3", "c_nl_Kme3",
    "a_nl_pST", "b_nl_pST", "c_nl_pST",
    "a_nl_Ksu", "b_nl_Ksu", "c_nl_Ksu"
  )

  yz_types <- c(
    "y", "z",
    "y_nl_O", "z_nl_O",
    "y_nl_N", "z_nl_N",
    "y_nl_Kme3", "z_nl_Kme3",
    "y_nl_pST", "z_nl_pST",
    "y_nl_Ksu", "z_nl_Ksu"
  )

  abc_filter <- df$scan_number == scan_num & df$PID == pid_num & df$ion_type %in% abc_types
  yz_filter <- df$scan_number == scan_num & df$PID == pid_num & df$ion_type %in% yz_types

  a_b_c_ions_df <- df[abc_filter, , drop = FALSE]
  y_z_ions_df <- df[yz_filter, , drop = FALSE]

  if (nrow(y_z_ions_df) > 0) {
    y_z_ions_df$ion_position <- nchar(y_z_ions_df$peptide) - y_z_ions_df$ion_position + 1
  }

  list(a_b_c_ions_df = a_b_c_ions_df, y_z_ions_df = y_z_ions_df)
}

#' Check for disambiguating ions before a modification site
#'
#' @param nearest_before Position of the nearest candidate before the modification.
#' @param mod_position Position of the modification under evaluation.
#' @param a_b_c_ions_df Data frame of a/b/c ion observations.
#' @param y_z_ions_df Data frame of y/z ion observations.
#'
#' @return Logical indicating whether a disambiguating ion is present.
#' @export
check_unambiguity_before <- function(nearest_before, mod_position, a_b_c_ions_df, y_z_ions_df) {
  if (is.na(nearest_before)) {
    return(TRUE)
  }

  abc_unambiguous <- any(a_b_c_ions_df$ion_position >= nearest_before &
                           a_b_c_ions_df$ion_position < mod_position, na.rm = TRUE)
  yz_unambiguous <- any(y_z_ions_df$ion_position > nearest_before &
                          y_z_ions_df$ion_position <= mod_position, na.rm = TRUE)

  abc_unambiguous || yz_unambiguous
}

#' Check for disambiguating ions after a modification site
#'
#' @inheritParams check_unambiguity_before
#'
#' @return Logical indicating whether a disambiguating ion is present.
#' @export
check_unambiguity_after <- function(nearest_after, mod_position, a_b_c_ions_df, y_z_ions_df) {
  if (is.na(nearest_after)) {
    return(TRUE)
  }

  abc_unambiguous <- any(a_b_c_ions_df$ion_position < nearest_after &
                           a_b_c_ions_df$ion_position >= mod_position, na.rm = TRUE)
  yz_unambiguous <- any(y_z_ions_df$ion_position <= nearest_after &
                          y_z_ions_df$ion_position > mod_position, na.rm = TRUE)

  abc_unambiguous || yz_unambiguous
}

#' Apply ambiguity checks for all modification columns
#'
#' @param df Data frame containing modification metadata and ion matches.
#' @param max_mods Maximum number of modification columns to evaluate.
#'
#' @return Data frame with ambiguity assessment columns populated (1 for TRUE, 0 for FALSE).
#' @export
apply_ambiguity_checks <- function(df, max_mods) {
  if (nrow(df) == 0) {
    return(df)
  }

  for (i in seq_len(nrow(df))) {
    ions <- assemble_ions_df(df, df$scan_number[i], df$PID[i])

    for (j in seq_len(max_mods)) {
      mod_pos_col <- paste0("mod_pos_", j)
      nearest_before_col <- paste0("nearest_before_mod_", j)
      nearest_after_col <- paste0("nearest_after_mod_", j)

      mod_position <- df[[mod_pos_col]][i]
      if (is.na(mod_position)) {
        next
      }

      before_value <- if (check_unambiguity_before(
        df[[nearest_before_col]][i],
        mod_position,
        ions$a_b_c_ions_df,
        ions$y_z_ions_df
      )) 1 else 0

      after_value <- if (check_unambiguity_after(
        df[[nearest_after_col]][i],
        mod_position,
        ions$a_b_c_ions_df,
        ions$y_z_ions_df
      )) 1 else 0

      df[[paste0("mod_pos_", j, "_abc_before_unambiguous")]][i] <- before_value
      df[[paste0("mod_pos_", j, "_yz_before_unambiguous")]][i] <- before_value
      df[[paste0("mod_pos_", j, "_abc_after_unambiguous")]][i] <- after_value
      df[[paste0("mod_pos_", j, "_yz_after_unambiguous")]][i] <- after_value
    }
  }

  df
}
