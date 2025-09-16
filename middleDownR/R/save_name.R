#' Construct an output filename from stored parameters
#'
#' @param stored_parameters Named list containing relevant parameters.
#'
#' @return Character string representing the formatted file name.
#' @export
save_name <- function(stored_parameters) {
  paste0(
    "Byonic_",
    "max_mods_", stored_parameters$max_mods, "_",
    "PEP.2D_", stored_parameters$PEP, "_",
    "min_delta_mod_", stored_parameters$min_delta_mod, "_",
    "ppm_tolerance_", stored_parameters$ppm_tolerance, "_",
    stored_parameters$aa_tag, "-tag"
  )
}
