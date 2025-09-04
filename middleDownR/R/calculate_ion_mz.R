#' Calculate m/z for ions
#'
#' @param peptide_mass Peptide mass.
#' @param charge Charge state.
#'
#' @return Calculated m/z value.
#' @keywords internal
calculate_ion_mz <- function(peptide_mass, charge) {
  (peptide_mass + atomic_masses["h"] * charge) / charge
}
