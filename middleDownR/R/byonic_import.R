#' Byonic result ingestion helpers
#'
#' Parsers that reshape Byonic modification summaries and organise spectrum-level lists prior to spectrum matching.
#' @keywords internal
NULL

#' Parse modification positions and masses from a Byonic string
#'
#' @param mod_str Character string describing modifications in "position/mass" format.
#'
#' @return Tibble with `position` and `variable` columns.
#' @export
parse_position_mass <- function(mod_str) {
  mod_list <- stringr::str_split(mod_str, ";\\s*")[[1]]
  mods_df <- purrr::map_dfr(mod_list, ~ {
    if (.x != "") {
      position <- as.integer(stringr::str_extract(.x, "\\d+"))
      variable <- as.numeric(stringr::str_extract(.x, "(?<=/)[\\d\\.]+"))
      tibble::tibble(position = position, variable = variable)
    } else {
      tibble::tibble(position = NA_integer_, variable = NA_real_)
    }
  })

  if (all(c(36, 37) %in% mods_df$position)) {
    sum_mass <- sum(mods_df$variable[mods_df$position %in% c(36, 37)], na.rm = TRUE)
    mods_df <- mods_df[!mods_df$position %in% c(37), ]

    if (any(mods_df$position == 36)) {
      mods_df <- dplyr::mutate(mods_df, variable = dplyr::if_else(position == 36, sum_mass, variable))
      mods_df <- dplyr::mutate(mods_df, position = dplyr::if_else(position == 36, 36L, position))
    } else {
      mods_df <- dplyr::bind_rows(mods_df, tibble::tibble(position = 36L, variable = sum_mass))
    }
  }

  mods_df
}

#' Re-index Byonic hits by scan and peptide identifier
#'
#' @param input_list List of Byonic identifications containing `scan_number` and `PID` entries.
#'
#' @return Nested list organised as `scan_number` -> `PID` -> hit data.
#' @export
process_Byonic_list <- function(input_list) {
  empty_list <- list()

  for (i in seq_along(input_list)) {
    scan_number <- as.character(input_list[[i]]$scan_number)
    PID <- as.character(input_list[[i]]$PID)
    scan_data <- input_list[[i]]
    empty_list[[as.character(scan_number)]][[PID]] <- scan_data
  }

  empty_list
}
