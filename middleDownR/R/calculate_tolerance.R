#' Calculate mass tolerance
#'
#' @param mz Numeric mass-to-charge value.
#' @param ppm Parts-per-million tolerance.
#'
#' @return Numeric tolerance value.
#' @examples
#' calculate_tolerance(1000, 10)
#' @export
calculate_tolerance <- function(mz, ppm) {
  as.numeric(mz) * ppm / 1e6
}
