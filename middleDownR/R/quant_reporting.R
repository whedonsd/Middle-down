#' Quantitation and reporting utilities
#'
#' Functions for translating numeric results into reportable labels and validating TMT reporter completeness.
#' @keywords internal
NULL

#' Map numeric modification/site values to descriptive labels
#'
#' @param input_value Numeric value to look up.
#' @param ref_df Reference data frame with columns `num_value` and `text_string`.
#'
#' @return Character label corresponding to the numeric value, or an empty string when unmatched.
#' @export
check_against_reference <- function(input_value, ref_df) {
  if (is.na(input_value)) {
    return("")
  }

  matches <- ref_df[["text_string"]][ref_df$num_value == input_value]
  if (length(matches) == 0) {
    ""
  } else {
    as.character(matches[1])
  }
}

#' Evaluate whether a row contains a minimum number of TMT channels
#'
#' @param data Data frame (typically a single-row tibble) containing intensity columns.
#' @param col_range_start Name of the first TMT intensity column.
#' @param col_range_end Name of the last TMT intensity column.
#' @param TMT_count Minimum number of channels that must contain signal.
#'
#' @return Logical value indicating whether the criterion is satisfied.
#' @export
TMT_check <- function(data, col_range_start, col_range_end, TMT_count) {
  columns <- names(data)
  start_idx <- match(col_range_start, columns)
  end_idx <- match(col_range_end, columns)

  if (is.na(start_idx) || is.na(end_idx)) {
    stop("Specified columns were not found in the data frame.")
  }

  if (start_idx > end_idx) {
    indices <- end_idx:start_idx
  } else {
    indices <- start_idx:end_idx
  }

  subset_df <- data[, indices, drop = FALSE]
  columns_check_result <- vapply(subset_df, function(x) any(x > 0, na.rm = TRUE), logical(1))
  sum(columns_check_result) >= TMT_count
}
