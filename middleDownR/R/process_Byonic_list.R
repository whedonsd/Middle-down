# Create list of mz intensity values
process_Byonic_list <- function(input_list) {
  # Initialize list to store values
  empty_list <- list()

  # Iterate through Scan_list using retrieved scan_number values
  for(i in seq_along(input_list)) {
    scan_number <- as.character(input_list[[i]]$scan_number)
    PID <- as.character(input_list[[i]]$PID)
    scan_data <- input_list[[i]]
    empty_list[[as.character(scan_number)]][[PID]] <- scan_data
  }
  
  # Explicitly return the list
  return(empty_list)
}
