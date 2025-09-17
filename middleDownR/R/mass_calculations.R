#' Mass calculation utilities
#'
#' Constants and helper functions for computing peptide, fragment, and tolerance values required during ion generation and matching.
#' @keywords internal
NULL

#' Monoisotopic masses of amino acids
#'
#' These values are used for fragment mass calculations.
#' @keywords internal
amino_acid_masses <- c(
  A = 71.03711, R = 156.10111, N = 114.04293, D = 115.02694,
  C = 103.00919, E = 129.04259, Q = 128.05858, G = 57.02146,
  H = 137.05891, I = 113.08406, L = 113.08406, K = 128.09496,
  M = 131.04049, F = 147.06841, P = 97.05276, S = 87.03203,
  T = 101.04768, W = 186.07931, Y = 163.06333, V = 99.06841
)

#' Atomic masses for elements (monoisotopic)
#' @keywords internal
atomic_masses <- c(
  C = 12.0, c = 13.0033548378, H = 1.007825, h = 1.007276,
  N = 14.003074, n = 15.0001088982, O = 15.994915,
  P = 30.973762, S = 31.972071
)

#' Adjustments for ion types
#' @keywords internal
ion_type_adjustments <- list(
  a = -(atomic_masses["C"] + atomic_masses["O"] + atomic_masses["H"]),
  b = 0,
  y = atomic_masses["H"] * (2) + atomic_masses["O"],
  c = atomic_masses["H"] + (atomic_masses["H"] * 2 + atomic_masses["N"]),
  z = atomic_masses["O"] - atomic_masses["N"],
  Ion_AmmoniaLossMass = atomic_masses["H"] * 3 + atomic_masses["N"],
  Ion_WaterLossMass = atomic_masses["H"] * 2 + atomic_masses["O"]
)

#' Calculate fragment mass from a peptide sequence
#'
#' @param sequence Amino acid sequence.
#' @param modifications Optional data frame with columns `position` and `variable`.
#'
#' @return Fragment mass.
#' @keywords internal
#' @export
calculate_fragment_mass_from_sequence <- function(sequence, modifications = NULL) {
  mass <- sum(sapply(strsplit(sequence, "")[[1]], function(aa) {
    if (aa %in% names(amino_acid_masses)) {
      amino_acid_masses[[aa]]
    } else {
      stop(paste("Undefined amino acid or modification:", aa))
    }
  }))

  if (!is.null(modifications) && nrow(modifications) > 0) {
    for (i in seq_len(nrow(modifications))) {
      mod <- modifications[i, ]
      if (mod$position <= nchar(sequence)) {
        mass <- mass + as.numeric(mod$variable)
        if (abs(mod$variable - 42.046950) < .Machine$double.eps^0.5) {
          mass <- mass - atomic_masses["h"]
        }
      }
    }
  }
  mass
}

#' Calculate m/z for ions
#'
#' @param peptide_mass Peptide mass.
#' @param charge Charge state.
#'
#' @return Calculated m/z value.
#' @keywords internal
#' @export
calculate_ion_mz <- function(peptide_mass, charge) {
  (peptide_mass + atomic_masses["h"] * charge) / charge
}

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
