#' Ion generation and matching pipeline
#'
#' Functions that create theoretical ion series, locate corresponding experimental peaks, and reshape match results.
#' @keywords internal
NULL

#' Generate ion masses for a peptide sequence
#'
#' @param sequence Amino acid sequence.
#' @param modifications_df Data frame of modifications with columns `position` and `variable`.
#'
#' @return List of ion m/z values by ion type and charge state.
#' @export
generate_ions_for_sequence <- function(sequence, modifications_df) {
  ions_list <- list()
  sequence_length <- nchar(sequence)

  for (charge in 1:5) {
    a_ions <- numeric(sequence_length - 1)
    a_nl_O_ions <- numeric(sequence_length - 1)
    a_nl_N_ions <- numeric(sequence_length - 1)
    a_nl_Kme3_ions <- numeric(sequence_length - 1)
    a_nl_pST_ions <- numeric(sequence_length - 1)
    a_nl_Ksu_ions <- numeric(sequence_length - 1)
    b_ions <- numeric(sequence_length - 1)
    b_nl_O_ions <- numeric(sequence_length - 1)
    b_nl_N_ions <- numeric(sequence_length - 1)
    b_nl_Kme3_ions <- numeric(sequence_length - 1)
    b_nl_pST_ions <- numeric(sequence_length - 1)
    b_nl_Ksu_ions <- numeric(sequence_length - 1)
    y_ions <- numeric(sequence_length - 1)
    y_nl_O_ions <- numeric(sequence_length - 1)
    y_nl_N_ions <- numeric(sequence_length - 1)
    y_nl_Kme3_ions <- numeric(sequence_length - 1)
    y_nl_pST_ions <- numeric(sequence_length - 1)
    y_nl_Ksu_ions <- numeric(sequence_length - 1)
    c_ions <- numeric(sequence_length - 1)
    c_nl_O_ions <- numeric(sequence_length - 1)
    c_nl_N_ions <- numeric(sequence_length - 1)
    c_nl_Kme3_ions <- numeric(sequence_length - 1)
    c_nl_pST_ions <- numeric(sequence_length - 1)
    c_nl_Ksu_ions <- numeric(sequence_length - 1)
    z_ions <- numeric(sequence_length - 1)
    z_nl_O_ions <- numeric(sequence_length - 1)
    z_nl_N_ions <- numeric(sequence_length - 1)
    z_nl_Kme3_ions <- numeric(sequence_length - 1)
    z_nl_pST_ions <- numeric(sequence_length - 1)
    z_nl_Ksu_ions <- numeric(sequence_length - 1)

    for (i in 1:(sequence_length - 1)) {
      a_seq <- substr(sequence, 1, i)
      b_seq <- substr(sequence, 1, i)
      y_seq <- substr(sequence, sequence_length - i + 1, sequence_length)
      c_seq <- substr(sequence, 1, i)
      z_seq <- substr(sequence, sequence_length - i + 1, sequence_length)

      a_mods <- modifications_df[modifications_df$position <= i, ]
      b_mods <- modifications_df[modifications_df$position <= i, ]
      c_mods <- modifications_df[modifications_df$position <= i, ]

      y_mods <- modifications_df
      if (nrow(y_mods) > 0) {
        y_mods$position <- sequence_length - y_mods$position + 1
        y_mods <- y_mods[y_mods$position <= i, ]
      }
      z_mods <- modifications_df
      if (nrow(z_mods) > 0) {
        z_mods$position <- sequence_length - z_mods$position + 1
        z_mods <- z_mods[z_mods$position <= i, ]
      }

      if (any(abs(a_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        specific_Kme3_loss <- 59.073499
        a_mass_specific_mod_loss <- a_mass - specific_Kme3_loss
        a_nl_Kme3_ions[i] <- calculate_ion_mz(a_mass_specific_mod_loss, charge)
      }
      if (any(abs(b_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        specific_Kme3_loss <- 59.073499
        b_mass_specific_mod_loss <- b_mass - specific_Kme3_loss
        b_nl_Kme3_ions[i] <- calculate_ion_mz(b_mass_specific_mod_loss, charge)
      }
      if (any(abs(y_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        specific_Kme3_loss <- 59.073499
        y_mass_specific_mod_loss <- y_mass - specific_Kme3_loss
        y_nl_Kme3_ions[i] <- calculate_ion_mz(y_mass_specific_mod_loss, charge)
      }
      if (any(abs(c_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        specific_Kme3_loss <- 59.073499
        c_mass_specific_mod_loss <- c_mass - specific_Kme3_loss
        c_nl_Kme3_ions[i] <- calculate_ion_mz(c_mass_specific_mod_loss, charge)
      }
      if (any(abs(z_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        specific_Kme3_loss <- 59.073499
        z_mass_specific_mod_loss <- z_mass - specific_Kme3_loss
        z_nl_Kme3_ions[i] <- calculate_ion_mz(z_mass_specific_mod_loss, charge)
      }
      if (any(abs(a_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        specific_pST_loss <- 97.976896
        a_mass_specific_mod_loss <- a_mass - specific_pST_loss
        a_nl_pST_ions[i] <- calculate_ion_mz(a_mass_specific_mod_loss, charge)
      }
      if (any(abs(b_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        specific_pST_loss <- 97.976896
        b_mass_specific_mod_loss <- b_mass - specific_pST_loss
        b_nl_pST_ions[i] <- calculate_ion_mz(b_mass_specific_mod_loss, charge)
      }
      if (any(abs(y_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        specific_pST_loss <- 97.976896
        y_mass_specific_mod_loss <- y_mass - specific_pST_loss
        y_nl_pST_ions[i] <- calculate_ion_mz(y_mass_specific_mod_loss, charge)
      }
      if (any(abs(c_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        specific_pST_loss <- 97.976896
        c_mass_specific_mod_loss <- c_mass - specific_pST_loss
        c_nl_pST_ions[i] <- calculate_ion_mz(c_mass_specific_mod_loss, charge)
      }
      if (any(abs(z_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        specific_pST_loss <- 97.976896
        z_mass_specific_mod_loss <- z_mass - specific_pST_loss
        z_nl_pST_ions[i] <- calculate_ion_mz(z_mass_specific_mod_loss, charge)
      }
      if (any(abs(a_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        specific_Ksu_loss <- 101.02387
        a_mass_specific_mod_loss <- a_mass - specific_Ksu_loss
        a_nl_Ksu_ions[i] <- calculate_ion_mz(a_mass_specific_mod_loss, charge)
      }
      if (any(abs(b_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        specific_Ksu_loss <- 101.02387
        b_mass_specific_mod_loss <- b_mass - specific_Ksu_loss
        b_nl_Ksu_ions[i] <- calculate_ion_mz(b_mass_specific_mod_loss, charge)
      }
      if (any(abs(y_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        specific_Ksu_loss <- 101.02387
        y_mass_specific_mod_loss <- y_mass - specific_Ksu_loss
        y_nl_Ksu_ions[i] <- calculate_ion_mz(y_mass_specific_mod_loss, charge)
      }
      if (any(abs(c_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        specific_Ksu_loss <- 101.02387
        c_mass_specific_mod_loss <- c_mass - specific_Ksu_loss
        c_nl_Ksu_ions[i] <- calculate_ion_mz(c_mass_specific_mod_loss, charge)
      }
      if (any(abs(z_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        specific_Ksu_loss <- 101.02387
        z_mass_specific_mod_loss <- z_mass - specific_Ksu_loss
        z_nl_Ksu_ions[i] <- calculate_ion_mz(z_mass_specific_mod_loss, charge)
      }

      a_mass <- calculate_fragment_mass_from_sequence(a_seq, a_mods) + ion_type_adjustments$a
      b_mass <- calculate_fragment_mass_from_sequence(b_seq, b_mods) + ion_type_adjustments$b
      y_mass <- calculate_fragment_mass_from_sequence(y_seq, y_mods) + ion_type_adjustments$y
      c_mass <- calculate_fragment_mass_from_sequence(c_seq, c_mods) + ion_type_adjustments$c
      z_mass <- calculate_fragment_mass_from_sequence(z_seq, z_mods) + ion_type_adjustments$z

      if (grepl("[KRQN]", a_seq)) {
        a_mass_N_loss <- a_mass - ion_type_adjustments$Ion_AmmoniaLossMass
        a_nl_N_ions[i] <- calculate_ion_mz(a_mass_N_loss, charge)
      }
      if (grepl("[STED]", a_seq)) {
        a_mass_O_loss <- a_mass - ion_type_adjustments$Ion_WaterLossMass
        a_nl_O_ions[i] <- calculate_ion_mz(a_mass_O_loss, charge)
      }
      if (grepl("[KRQN]", b_seq)) {
        b_mass_N_loss <- b_mass - ion_type_adjustments$Ion_AmmoniaLossMass
        b_nl_N_ions[i] <- calculate_ion_mz(b_mass_N_loss, charge)
      }
      if (grepl("[STED]", b_seq)) {
        b_mass_O_loss <- b_mass - ion_type_adjustments$Ion_WaterLossMass
        b_nl_O_ions[i] <- calculate_ion_mz(b_mass_O_loss, charge)
      }
      if (grepl("[KRQN]", y_seq)) {
        y_mass_N_loss <- y_mass - ion_type_adjustments$Ion_AmmoniaLossMass
        y_nl_N_ions[i] <- calculate_ion_mz(y_mass_N_loss, charge)
      }
      if (grepl("[STED]", y_seq)) {
        y_mass_O_loss <- y_mass - ion_type_adjustments$Ion_WaterLossMass
        y_nl_O_ions[i] <- calculate_ion_mz(y_mass_O_loss, charge)
      }
      if (grepl("[KRQN]", c_seq)) {
        c_mass_N_loss <- c_mass - ion_type_adjustments$Ion_AmmoniaLossMass
        c_nl_N_ions[i] <- calculate_ion_mz(c_mass_N_loss, charge)
      }
      if (grepl("[STED]", c_seq)) {
        c_mass_O_loss <- c_mass - ion_type_adjustments$Ion_WaterLossMass
        c_nl_O_ions[i] <- calculate_ion_mz(c_mass_O_loss, charge)
      }
      if (grepl("[KRQN]", z_seq)) {
        z_mass_N_loss <- z_mass - ion_type_adjustments$Ion_AmmoniaLossMass
        z_nl_N_ions[i] <- calculate_ion_mz(z_mass_N_loss, charge)
      }
      if (grepl("[STED]", z_seq)) {
        z_mass_O_loss <- z_mass - ion_type_adjustments$Ion_WaterLossMass
        z_nl_O_ions[i] <- calculate_ion_mz(z_mass_O_loss, charge)
      }

      a_ions[i] <- calculate_ion_mz(a_mass, charge)
      b_ions[i] <- calculate_ion_mz(b_mass, charge)
      y_ions[i] <- calculate_ion_mz(y_mass, charge)
      c_ions[i] <- calculate_ion_mz(c_mass, charge)
      z_ions[i] <- calculate_ion_mz(z_mass, charge)
    }

    ions_list[[paste("a_ions_charge", charge)]] <- a_ions
    ions_list[[paste("a_nl_O_ions_charge", charge)]] <- a_nl_O_ions
    ions_list[[paste("a_nl_N_ions_charge", charge)]] <- a_nl_N_ions
    ions_list[[paste("a_nl_Kme3_ions_charge", charge)]] <- a_nl_Kme3_ions
    ions_list[[paste("a_nl_pST_ions_charge", charge)]] <- a_nl_pST_ions
    ions_list[[paste("a_nl_Ksu_ions_charge", charge)]] <- a_nl_Ksu_ions
    ions_list[[paste("b_ions_charge", charge)]] <- b_ions
    ions_list[[paste("b_nl_O_ions_charge", charge)]] <- b_nl_O_ions
    ions_list[[paste("b_nl_N_ions_charge", charge)]] <- b_nl_N_ions
    ions_list[[paste("b_nl_Kme3_ions_charge", charge)]] <- b_nl_Kme3_ions
    ions_list[[paste("b_nl_pST_ions_charge", charge)]] <- b_nl_pST_ions
    ions_list[[paste("b_nl_Ksu_ions_charge", charge)]] <- b_nl_Ksu_ions
    ions_list[[paste("y_ions_charge", charge)]] <- y_ions
    ions_list[[paste("y_nl_O_ions_charge", charge)]] <- y_nl_O_ions
    ions_list[[paste("y_nl_N_ions_charge", charge)]] <- y_nl_N_ions
    ions_list[[paste("y_nl_Kme3_ions_charge", charge)]] <- y_nl_Kme3_ions
    ions_list[[paste("y_nl_pST_ions_charge", charge)]] <- y_nl_pST_ions
    ions_list[[paste("y_nl_Ksu_ions_charge", charge)]] <- y_nl_Ksu_ions
    ions_list[[paste("c_ions_charge", charge)]] <- c_ions
    ions_list[[paste("c_nl_O_ions_charge", charge)]] <- c_nl_O_ions
    ions_list[[paste("c_nl_N_ions_charge", charge)]] <- c_nl_N_ions
    ions_list[[paste("c_nl_Kme3_ions_charge", charge)]] <- c_nl_Kme3_ions
    ions_list[[paste("c_nl_pST_ions_charge", charge)]] <- c_nl_pST_ions
    ions_list[[paste("c_nl_Ksu_ions_charge", charge)]] <- c_nl_Ksu_ions
    ions_list[[paste("z_ions_charge", charge)]] <- z_ions
    ions_list[[paste("z_nl_O_ions_charge", charge)]] <- z_nl_O_ions
    ions_list[[paste("z_nl_N_ions_charge", charge)]] <- z_nl_N_ions
    ions_list[[paste("z_nl_Kme3_ions_charge", charge)]] <- z_nl_Kme3_ions
    ions_list[[paste("z_nl_pST_ions_charge", charge)]] <- z_nl_pST_ions
    ions_list[[paste("z_nl_Ksu_ions_charge", charge)]] <- z_nl_Ksu_ions
  }

  ions_list
}

#' Find ion matches with intensity information
#'
#' @param ions_df Data frame of theoretical ion m/z values.
#' @param scan_data Data frame with columns `mz` and `intensity`.
#' @param ppm_tolerance Mass tolerance in ppm.
#' @param sequence Peptide sequence for positional filtering.
#'
#' @return List of matched ions with m/z and intensity.
#' @export
find_matches_with_intensity <- function(ions_df, scan_data, ppm_tolerance, sequence) {
  matches <- list()
  sequence_length <- nchar(sequence)
  max_intensity <- max(scan_data[, "intensity"], na.rm = TRUE)

  ion_columns <- names(ions_df)

  for (ion_column in ion_columns) {
    charge_state <- as.numeric(gsub(".*charge[ .](\\d+)$", "\\1", ion_column))
    ion_vector <- ions_df[[ion_column]]

    for (i in seq_along(ion_vector)) {
      ion_mass <- ion_vector[i]

      if ((charge_state == 1 && i > (sequence_length / 4)) ||
          (charge_state == 2 && i > (sequence_length / 3)) ||
          (charge_state == 3 && i > (sequence_length / 3)) ||
          (charge_state == 4 && i < 4) ||
          (charge_state == 5 && i < (sequence_length / 3)) ||
          (charge_state == 6 && i < (sequence_length / 2))) {
        next
      }

      ion_tolerance <- calculate_tolerance(ion_mass, ppm_tolerance)
      all_matching_indices <- which(abs(scan_data[, "mz"] - ion_mass) <= ion_tolerance)

      if (length(all_matching_indices) > 0) {
        distances <- abs(scan_data[all_matching_indices, "mz"] - ion_mass)
        closest_match_index <- all_matching_indices[which.min(distances)]

        expected_isotopic_mz <- scan_data[closest_match_index, "mz"] + (1 / charge_state)
        isotopic_tolerance <- calculate_tolerance(expected_isotopic_mz, ppm_tolerance)
        isotopic_peak_indices <- which(abs(scan_data[, "mz"] - expected_isotopic_mz) <= isotopic_tolerance)

        if (length(isotopic_peak_indices) > 0) {
          distances <- abs(scan_data[isotopic_peak_indices, "mz"] - expected_isotopic_mz)
          nearest_peak_index <- isotopic_peak_indices[which.min(distances)]

          if (scan_data[nearest_peak_index, "intensity"] <= scan_data[closest_match_index, "intensity"]) {
            match_key <- paste(ion_column,
                               "pos",
                               i,
                               "mz",
                               round(ion_mass, 4),
                               "charge",
                               charge_state,
                               sep = "_")
            matches[[match_key]] <- list(
              mz = scan_data[closest_match_index, "mz"],
              intensity = scan_data[closest_match_index, "intensity"]
            )
          }
        }
      }
    }
  }

  matches
}

#' Match theoretical ions against experimental scans
#'
#' @param list1_chunk Nested list of Byonic hits for a subset of scans.
#' @param list2 List of scan data frames keyed by scan number.
#'
#' @return Nested list of matches with metadata.
#' @export
list_match <- function(list1_chunk, list2) {
  matches_data <- list()

  for (scan_number in names(list1_chunk)) {
    scan_data_list <- list1_chunk[[scan_number]]

    for (PID in names(scan_data_list)) {
      hit_data <- scan_data_list[[PID]]
      sequence <- hit_data[["sequence"]]
      modifications_df <- as.data.frame(hit_data[["modifications_df"]])

      scan_data <- list2[[scan_number]]
      ions_df <- generate_ions_for_sequence(sequence, modifications_df)
      matches_df <- find_matches_with_intensity(ions_df, scan_data, ppm_tolerance, sequence)

      matches_df$peptide_sequence <- sequence
      modifications_str <- apply(modifications_df, 1, function(x) paste(x["position"], x["variable"], sep = ":"))
      matches_df$modifications <- paste(modifications_str, collapse = "; ")

      if (is.null(matches_data[[scan_number]])) {
        matches_data[[scan_number]] <- list()
      }
      matches_data[[scan_number]][[PID]] <- matches_df
    }
  }

  matches_data
}

#' Convert a list chunk of ion matches into a data frame
#'
#' @param chunk Nested list of matches for a subset of scans.
#'
#' @return Data frame containing ion match details for the chunk.
#' @export
process_chunk <- function(chunk) {
  matches_list_df_chunk <- data.frame(
    scan_number = integer(),
    PID = integer(),
    ion_type = character(),
    ion_position = integer(),
    ion_charge = integer(),
    mz = numeric(),
    intensity = numeric(),
    peptide = character(),
    modifications = character(),
    stringsAsFactors = FALSE
  )

  rows_to_add <- vector("list", length = 0L)

  for (scan_number in names(chunk)) {
    for (PID in names(chunk[[scan_number]])) {
      hit_data <- chunk[[scan_number]][[PID]]

      peptide_sequence <- hit_data[["peptide_sequence"]]
      modifications <- hit_data[["modifications"]]
      ion_keys <- setdiff(names(hit_data), c("peptide_sequence", "modifications"))

      for (ion_key in ion_keys) {
        ion_data <- hit_data[[ion_key]]

        ion_type <- sub("^(.*?)_ions.*$", "\\1", ion_key)
        parts <- strsplit(ion_key, "_pos_")[[1]]
        ion_position <- as.integer(strsplit(parts[2], "_")[[1]][1])

        charge_pattern <- "_charge[ .]?([0-9]+)_"
        charge_matches <- gregexpr(charge_pattern, ion_key)[[1]]
        if (length(charge_matches) > 0 && charge_matches[1] != -1) {
          matches <- regmatches(ion_key, gregexpr(charge_pattern, ion_key))
          ion_charge <- as.integer(sub(charge_pattern, "\\1", matches[[1]]))
        } else {
          ion_charge <- NA_integer_
        }

        rows_to_add[[length(rows_to_add) + 1L]] <- data.frame(
          scan_number = as.integer(scan_number),
          PID = as.integer(PID),
          ion_type = ion_type,
          ion_position = ion_position,
          ion_charge = ion_charge,
          mz = ion_data$mz,
          intensity = ion_data$intensity,
          peptide = peptide_sequence,
          modifications = modifications,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows_to_add) == 0) {
    return(matches_list_df_chunk)
  }

  do.call(rbind, rows_to_add)
}
