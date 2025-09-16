#' Retrieve the first available ratio value based on rank order
#'
#' @param rank_order Numeric rank for the peptide-spectrum match.
#' @param average_12binary_ratio Numeric ratio comparing ranks 1 and 2.
#' @param average_13binary_ratio Numeric ratio comparing ranks 1 and 3.
#' @param average_ratio Fallback average ratio across ions.
#'
#' @return Numeric value representing the selected ratio.
#' @export
retrieve_value <- function(rank_order, average_12binary_ratio, average_13binary_ratio, average_ratio) {
  candidates <- switch(
    as.character(rank_order),
    `1` = c(average_13binary_ratio, average_12binary_ratio, average_ratio, 1),
    `2` = c(average_12binary_ratio, average_ratio, 1),
    `3` = c(average_13binary_ratio, average_ratio, 1),
    numeric(0)
  )

  for (value in candidates) {
    if (!is.null(value) && !is.na(value) && !is.nan(value)) {
      return(value)
    }
  }

  NA_real_
}
