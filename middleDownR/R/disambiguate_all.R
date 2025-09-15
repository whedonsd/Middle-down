# Function to iterate over all scan_number and PID combinations and apply check_columns
# Applying disambiguate to all PTMs for a given peptide lets us establish whether the assignment has enough ions to be unambiguous
disambiguate_all <- function(data, col_range_start, col_range_end) {
  results <- data %>%
    group_by(scan_number, PID) %>%
    summarise(all_true = disambiguate(cur_data(), col_range_start, col_range_end), .groups = 'drop')
  
  return(results)
}
