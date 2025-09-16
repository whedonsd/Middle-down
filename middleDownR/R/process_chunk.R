#' Convert a list chunk of ion matches into a data frame
#'
#' @param chunk Nested list of matches for a subset of scans.
#'
#' @return Data frame containing ion match details for the chunk.
#' @export
process_chunk <- function(chunk) {
  matches_list_df_chunk <- data.frame(
    scan_number = integer(),
    PID = integer(),
    ion_type = character(),
    ion_position = integer(),
    ion_charge = integer(),
    mz = numeric(),
    intensity = numeric(),
    peptide = character(),
    modifications = character(),
    stringsAsFactors = FALSE
  )

  rows_to_add <- vector("list", length = 0L)

  for (scan_number in names(chunk)) {
    for (PID in names(chunk[[scan_number]])) {
      hit_data <- chunk[[scan_number]][[PID]]

      peptide_sequence <- hit_data[["peptide_sequence"]]
      modifications <- hit_data[["modifications"]]
      ion_keys <- setdiff(names(hit_data), c("peptide_sequence", "modifications"))

      for (ion_key in ion_keys) {
        ion_data <- hit_data[[ion_key]]

        ion_type <- sub("^(.*?)_ions.*$", "\\1", ion_key)
        parts <- strsplit(ion_key, "_pos_")[[1]]
        ion_position <- as.integer(strsplit(parts[2], "_")[[1]][1])

        charge_pattern <- "_charge[ .]?([0-9]+)_"
        charge_matches <- gregexpr(charge_pattern, ion_key)[[1]]
        if (length(charge_matches) > 0 && charge_matches[1] != -1) {
          matches <- regmatches(ion_key, gregexpr(charge_pattern, ion_key))
          ion_charge <- as.integer(sub(charge_pattern, "\\1", matches[[1]]))
        } else {
          ion_charge <- NA_integer_
        }

        rows_to_add[[length(rows_to_add) + 1L]] <- data.frame(
          scan_number = as.integer(scan_number),
          PID = as.integer(PID),
          ion_type = ion_type,
          ion_position = ion_position,
          ion_charge = ion_charge,
          mz = ion_data$mz,
          intensity = ion_data$intensity,
          peptide = peptide_sequence,
          modifications = modifications,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows_to_add) == 0) {
    return(matches_list_df_chunk)
  }

  do.call(rbind, rows_to_add)
}
