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
