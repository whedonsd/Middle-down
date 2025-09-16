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
