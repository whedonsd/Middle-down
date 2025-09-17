#' Modification column helpers
#'
#' Tools for expanding stored modification summaries into position/mass columns used by later localisation steps.
#' @keywords internal
NULL

#' Populate modification columns from summary strings
#'
#' @param df Data frame containing a `modifications` column.
#' @param max_mods Maximum number of modification columns to fill.
#'
#' @return Data frame with modification position and mass columns populated.
#' @export
fill_modifications_direct <- function(df, max_mods) {
  for (row in seq_len(nrow(df))) {
    modification_str <- df$modifications[row]
    if (is.na(modification_str)) {
      next
    }

    mod_parts <- strsplit(modification_str, ";", fixed = TRUE)[[1]]
    mods <- lapply(mod_parts, function(x) strsplit(x, ":", fixed = TRUE)[[1]])

    for (i in seq_along(mods)) {
      if (i > max_mods) {
        break
      }
      mod_pair <- mods[[i]]
      df[row, paste0("mod_pos_", i)] <- as.numeric(mod_pair[1])
      df[row, paste0("mod_mass_", i)] <- as.numeric(mod_pair[2])
    }
  }
  df
}

#' Add a column to a data frame if it is absent
#'
#' @param df Data frame to modify.
#' @param column_name Column name to add when missing.
#'
#' @return Data frame guaranteed to contain the requested column.
#' @export
add_column_if_not_exists <- function(df, column_name) {
  if (!column_name %in% names(df)) {
    df[[column_name]] <- NA
  }
  df
}
