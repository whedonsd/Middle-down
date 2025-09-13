# Function to process scans from thermo .RAW data converted to mzML.gz
# Create list of mz intensity values
process_scans_from_mzMLgz <- function(Scan_list, mz) {
  # Initialize list to store values
  mz_list <- list()

  # Iterate through Scan_list using retrieved scan_number values
  for(i in seq_along(Scan_list)) {
    scan_number <- as.integer(Scan_list[[i]])
    scan_data <- as.data.frame(peaks(mz, scan = scan_number))
    mz_list[[as.character(scan_number)]] <- scan_data
  }
  
  # Explicitly return the list
  return(mz_list)
}
