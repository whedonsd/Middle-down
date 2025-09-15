# Simplified check_columns function
# Columns contain ions that fall between two sites that could have the same PTM (e.g. K and R for methyl, K and K for acetyl)
# This function checks whether at least one ion falls between the putative modifications site and the preceding and following candidate sites (e.g. for K14ac check K9 and K18)
disambiguate <- function(data, col_range_start, col_range_end) {
  # Assuming data is already filtered for a specific scan_number and PID
  
  # Dynamically select the range of columns if col_range_start and col_range_end are column names
  columns_to_check <- select(data, dplyr::all_of(col_range_start):dplyr::all_of(col_range_end))
  
  # Apply the check across the specified columns
  columns_check_result <- purrr::map_lgl(columns_to_check, function(column) {
    # Check if the column contains at least one "1"
    has_one <- any(column == 1, na.rm = TRUE)
    # If the column has at least one "1", it satisfies the condition
    if (has_one) {
      return(TRUE)
    } else {
      # If there are no "1"s, check if there's at least one NA (which also satisfies the condition)
      return(all(is.na(column)))
    }
  })
  
  # Return TRUE if all columns meet the criteria
  return(all(columns_check_result))
}
