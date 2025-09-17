#' Spectra processing helpers
#'
#' Utilities for turning mzML scans into analysis-ready data structures and for retrieving targeted signal intensities.
#' @keywords internal
NULL

#' Process scans from an mzML file
#'
#' @param Scan_list Numeric or character vector of scan numbers to extract.
#' @param mz An `mzR` object representing an mzML file.
#'
#' @return List of data frames containing m/z and intensity columns for each scan.
#' @export
process_scans_from_mzMLgz <- function(Scan_list, mz) {
  mz_list <- list()

  for (i in seq_along(Scan_list)) {
    scan_number <- as.integer(Scan_list[[i]])
    scan_data <- as.data.frame(mzR::peaks(mz, scan = scan_number))
    mz_list[[as.character(scan_number)]] <- scan_data
  }

  mz_list
}

#' Locate intensities for target m/z values within a scan
#'
#' @param scan_number Scan number to inspect.
#' @param mz mzR object containing spectra.
#' @param target_mz Numeric vector of target m/z values.
#' @param ppm_tolerance Mass tolerance in parts-per-million.
#'
#' @return Numeric vector of intensities corresponding to the closest matches.
#' @export
find_intensities_for_scan <- function(scan_number, mz, target_mz, ppm_tolerance) {
  scan_number_int <- as.integer(scan_number)
  scan_data <- mzR::peaks(mz, scan = scan_number_int)

  vapply(target_mz, function(target) {
    tolerance <- calculate_tolerance(target, ppm_tolerance)
    within_tolerance <- which(scan_data[, 1] >= target - tolerance & scan_data[, 1] <= target + tolerance)

    if (length(within_tolerance) == 0) {
      return(NA_real_)
    }

    distances <- abs(scan_data[within_tolerance, 1] - target)
    closest_index <- within_tolerance[which.min(distances)]
    scan_data[closest_index, 2]
  }, numeric(1))
}
