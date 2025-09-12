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
  
  # Identify ion type columns in the dataframe
  ion_columns <- names(ions_df)
  
  for (ion_column in ion_columns) {
    # Determine the charge state based on the column name
    charge_state <- as.numeric(gsub(".*charge[ .](\\d+)$", "\\1", ion_column))
    ion_vector <- ions_df[[ion_column]]
    
    for (i in seq_along(ion_vector)) {
      ion_mass <- ion_vector[i]

      # Early skip conditions based on charge state and position in sequence
      if ((charge_state == 1 && i > (sequence_length / 4)) ||
          (charge_state == 2 && i > (sequence_length / 3)) ||
          (charge_state == 3 && i > (sequence_length / 3)) ||
          (charge_state == 4 && i < 4) ||
          (charge_state == 5 && i < (sequence_length / 3)) ||
          (charge_state == 6 && i < (sequence_length / 2))) {
        next
      }
      
      ion_tolerance <- calculate_tolerance(ion_mass, ppm_tolerance)
      
      # Find all potential matching indices
      all_matching_indices <- which(abs(scan_data[, "mz"] - ion_mass) <= ion_tolerance)

      if (length(all_matching_indices) > 0) {
        # Calculate distances for all matches
        distances <- abs(scan_data[all_matching_indices, "mz"] - ion_mass)
        # Select the index with the smallest distance
        closest_match_index <- all_matching_indices[which.min(distances)]

        # Calculate expected isotopic m/z and its tolerance
        expected_isotopic_mz <- scan_data[closest_match_index, "mz"] + (1 / charge_state)
        isotopic_tolerance <- calculate_tolerance(expected_isotopic_mz, ppm_tolerance)
        isotopic_peak_indices <- which(abs(scan_data[, "mz"] - expected_isotopic_mz) <= isotopic_tolerance)

        if (length(isotopic_peak_indices) > 0) {
          distances <- abs(scan_data[isotopic_peak_indices, "mz"] - expected_isotopic_mz)
          nearest_peak_index <- isotopic_peak_indices[which.min(distances)]

          if (scan_data[nearest_peak_index, "intensity"] <= scan_data[closest_match_index, "intensity"]) {
            # Construct a key for this match, including the column name and ion position
            match_key <- paste(ion_column,
                               "pos",
                               i,
                               "mz",
                               round(ion_mass, 4),
                               "charge",
                               charge_state,
                               sep = "_")
          # Store match with its intensity
            matches[[match_key]] <- list(mz = scan_data[closest_match_index, "mz"],
                                         intensity = scan_data[closest_match_index, "intensity"]
                                         )
          }
        }
      }
    }
  }

  return(matches)
}
