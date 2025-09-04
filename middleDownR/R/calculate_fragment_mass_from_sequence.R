#' Calculate fragment mass from a peptide sequence
#'
#' @param sequence Amino acid sequence.
#' @param modifications Optional data frame with columns `position` and `variable`.
#'
#' @return Fragment mass.
#' @keywords internal
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
