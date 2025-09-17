#' Ion intensity aggregation helpers
#'
#' Shared utilities for merging hit tables, extracting ion intensities, and choosing representative quantitative values.
#' @keywords internal
NULL

#' Merge hit data frames while preserving character columns
#'
#' @param appended_df Data frame to append.
#' @param main_df Existing data frame to extend.
#'
#' @return Combined data frame with consistent column classes.
#' @export
conditional_merge_hit_df <- function(appended_df, main_df) {
  if (nrow(appended_df) == 0) {
    return(main_df)
  }

  character_cols <- c("a_ions_shared", "b_ions_shared", "y_ions_shared", "c_ions_shared", "z_ions_shared")

  for (col in intersect(character_cols, names(appended_df))) {
    appended_df[[col]] <- as.character(appended_df[[col]])
  }

  for (col in intersect(character_cols, names(main_df))) {
    main_df[[col]] <- as.character(main_df[[col]])
  }

  combined <- rbind(main_df, appended_df)
  rownames(combined) <- NULL
  combined
}

#' Extract intensities for a set of ions within a grouped sublist
#'
#' @param sublist Data frame containing ion observations for a single scan/PID.
#' @param ions_list Character vector of ion identifiers of interest.
#'
#' @return Numeric vector of matched intensities.
#' @export
extract_intensities_from_sublist <- function(sublist, ions_list) {
  if (is.null(sublist) || nrow(sublist) == 0) {
    return(numeric(0))
  }

  matches <- sublist[sublist$ion_type_position_charge %in% ions_list, , drop = FALSE]
  if (nrow(matches) == 0) {
    return(numeric(0))
  }

  matches$intensity
}

#' Extract intensities for a scan/PID across grouped ion lists
#'
#' @param scan_number Scan number to locate.
#' @param PID Peptide identifier to locate.
#' @param ions_list Character vector of ion identifiers.
#' @param grouped_list List of data frames grouped by scan and PID.
#'
#' @return Numeric vector of matched intensities.
#' @export
extract_intensities <- function(scan_number, PID, ions_list, grouped_list) {
  if (length(grouped_list) == 0) {
    return(numeric(0))
  }

  match_index <- which(vapply(
    grouped_list,
    function(x) {
      any(x$scan_number == scan_number & x$PID == PID)
    },
    logical(1)
  ))

  if (length(match_index) == 0) {
    return(numeric(0))
  }

  extract_intensities_from_sublist(grouped_list[[match_index[1]]], ions_list)
}

#' Extract ion identities from a grouped sublist
#'
#' @inheritParams extract_intensities_from_sublist
#'
#' @return Character vector of ion identifiers present in the sublist.
#' @export
extract_ion_identities_from_sublist <- function(sublist, ions_list) {
  if (is.null(sublist) || nrow(sublist) == 0) {
    return(character(0))
  }

  matches <- sublist[sublist$ion_type_position_charge %in% ions_list, , drop = FALSE]
  if (nrow(matches) == 0) {
    return(character(0))
  }

  matches$ion_type_position_charge
}

#' Extract ordered ion identifiers for a scan/PID across grouped ion lists
#'
#' @inheritParams extract_intensities
#'
#' @return Character vector of matched ion identifiers.
#' @export
extract_ordered_ion_lists <- function(scan_number, PID, ions_list, grouped_list) {
  if (length(grouped_list) == 0) {
    return(character(0))
  }

  match_index <- which(vapply(
    grouped_list,
    function(x) {
      any(x$scan_number == scan_number & x$PID == PID)
    },
    logical(1)
  ))

  if (length(match_index) == 0) {
    return(character(0))
  }

  extract_ion_identities_from_sublist(grouped_list[[match_index[1]]], ions_list)
}

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
