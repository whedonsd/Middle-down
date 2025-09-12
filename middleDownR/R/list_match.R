# Function to find ions in the raw data that match the PSM called by Byonic
list_match <- function(list1_chunk, list2) {
  matches_data <- list()  # Initialize outside the loops
  
  # Iterate through the chunk of list1 using names for direct access
  for (scan_number in names(list1_chunk)) {
    # Direct access to each scan_number's data in list1
    scan_data_list <- list1_chunk[[scan_number]]
    
    # Iterate through hits within each scan_number
    for (PID in names(scan_data_list)) {
      hit_data <- scan_data_list[[PID]]
      sequence <- hit_data[["sequence"]]
      modifications_df <- as.data.frame(hit_data[["modifications_df"]])  
  
      # Retrieve m/z and intensity data from list2 using scan_number
      scan_data <- list2[[scan_number]]
      
      # Generate ions for the sequence with its modifications
      ions_df <- generate_ions_for_sequence(sequence, modifications_df)
      
      # Find matching ions in the mass spectrometry data
      matches_df <- find_matches_with_intensity(ions_df, scan_data, ppm_tolerance, sequence)
    
      # Enrich matches_df with additional data
      matches_df$peptide_sequence <- sequence
      modifications_str <- apply(modifications_df, 1, function(x) paste(x['position'], x['variable'], sep=":"))
      matches_df$modifications <- paste(modifications_str, collapse="; ")
      
      # Store the match data indexed by scan_number and PID
      if (!exists(scan_number, matches_data)) {
        matches_data[[scan_number]] <- list()
      }
      matches_data[[scan_number]][[PID]] <- matches_df
    }
  }
  
  return(matches_data)
}
