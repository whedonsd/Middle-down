#' Nearest-site annotation helpers
#'
#' Routines that compute proximal candidate residues for each modification and record the rule set used.
#' @keywords internal
NULL

#' Locate nearest candidate residues for modifications at K/R sites
#'
#' @param peptide Character vector of peptide sequences.
#' @param mod_position Numeric vector of modification positions.
#'
#' @return Data frame with nearest residues before and after the modification.
#' @export
find_nearest_KR_df <- function(peptide, mod_position) {
  n <- length(peptide)
  nearest_before <- rep(NA_real_, n)
  nearest_after <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    aa_sequence <- strsplit(peptide[i], "")[[1]]
    kr_positions <- which(aa_sequence %in% c("K", "R"))
    current_position <- mod_position[min(i, length(mod_position))]
    prior <- kr_positions[kr_positions < current_position]
    after <- kr_positions[kr_positions > current_position]
    if (length(prior) > 0) {
      nearest_before[i] <- max(prior)
    }
    if (length(after) > 0) {
      nearest_after[i] <- min(after)
    }
  }

  data.frame(nearest_before = nearest_before, nearest_after = nearest_after)
}

#' Locate nearest candidate residues for modifications at K sites
#'
#' @inheritParams find_nearest_KR_df
#'
#' @return Data frame with nearest K residues before and after the modification.
#' @export
find_nearest_K_df <- function(peptide, mod_position) {
  n <- length(peptide)
  nearest_before <- rep(NA_real_, n)
  nearest_after <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    aa_sequence <- strsplit(peptide[i], "")[[1]]
    k_positions <- which(aa_sequence %in% c("K"))
    current_position <- mod_position[min(i, length(mod_position))]
    prior <- k_positions[k_positions < current_position]
    after <- k_positions[k_positions > current_position]
    if (length(prior) > 0) {
      nearest_before[i] <- max(prior)
    }
    if (length(after) > 0) {
      nearest_after[i] <- min(after)
    }
  }

  data.frame(nearest_before = nearest_before, nearest_after = nearest_after)
}

#' Locate nearest candidate residues for modifications at S/T sites
#'
#' @inheritParams find_nearest_KR_df
#'
#' @return Data frame with nearest S/T residues before and after the modification.
#' @export
find_nearest_ST_df <- function(peptide, mod_position) {
  n <- length(peptide)
  nearest_before <- rep(NA_real_, n)
  nearest_after <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    aa_sequence <- strsplit(peptide[i], "")[[1]]
    st_positions <- which(aa_sequence %in% c("S", "T"))
    current_position <- mod_position[min(i, length(mod_position))]
    prior <- st_positions[st_positions < current_position]
    after <- st_positions[st_positions > current_position]
    if (length(prior) > 0) {
      nearest_before[i] <- max(prior)
    }
    if (length(after) > 0) {
      nearest_after[i] <- min(after)
    }
  }

  data.frame(nearest_before = nearest_before, nearest_after = nearest_after)
}

#' Apply nearest-residue lookups across all modifications in a data frame
#'
#' @param df Data frame containing peptide sequences and modification columns.
#' @param max_mods Maximum number of modification columns to process.
#'
#' @return Data frame with nearest residue columns populated.
#' @export
apply_nearest_function_directly <- function(df, max_mods) {
  if (nrow(df) == 0) {
    return(df)
  }

  for (i in seq_len(nrow(df))) {
    peptide <- df$peptide[i]
    mod_positions <- vapply(seq_len(max_mods), function(n) df[[paste0("mod_pos_", n)]][i], numeric(1))
    mod_masses <- vapply(seq_len(max_mods), function(n) df[[paste0("mod_mass_", n)]][i], numeric(1))

    for (j in seq_len(max_mods)) {
      mod_position <- mod_positions[j]
      mod_mass <- mod_masses[j]

      if (is.na(mod_position) || is.na(mod_mass)) {
        next
      }

      if (mod_mass %in% c(14.015650, 28.031300)) {
        result <- find_nearest_KR_df(peptide, mod_position)
        nearest_fn <- "KR"
      } else if (mod_mass %in% c(79.966331)) {
        result <- find_nearest_ST_df(peptide, mod_position)
        nearest_fn <- "ST"
      } else if (mod_mass %in% c(42.010565, 42.046950, 56.026215, 72.02113, 86.036779, 100.01604)) {
        result <- find_nearest_K_df(peptide, mod_position)
        nearest_fn <- "K"
      } else {
        result <- data.frame(nearest_before = NA_real_, nearest_after = NA_real_)
        nearest_fn <- "NA"
      }

      df[[paste0("nearest_before_mod_", j)]][i] <- result$nearest_before[1]
      df[[paste0("nearest_after_mod_", j)]][i] <- result$nearest_after[1]
      df[[paste0("nearest_fn_mod_", j)]][i] <- nearest_fn
    }
  }

  df
}
