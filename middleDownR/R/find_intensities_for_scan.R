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
