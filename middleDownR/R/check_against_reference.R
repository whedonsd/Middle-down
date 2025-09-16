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
