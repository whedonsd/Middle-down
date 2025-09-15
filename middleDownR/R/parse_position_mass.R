# Define a function to parse position and mass values and create a dataframe
parse_position_mass <- function(mod_str) {
  mod_list <- str_split(mod_str, ";\\s*")[[1]]
  mods_df <- map_df(mod_list, ~ {
    if (.x != "") {
      # Extract only the numeric part of the position
      position <- as.integer(str_extract(.x, "\\d+"))
      variable <- as.numeric(str_extract(.x, "(?<=/)[\\d\\.]+"))
      tibble(position = position, variable = variable)
    } else {
      tibble(position = NA_integer_, variable = NA_real_)
    }
  })
  
  # Check if both positions 36 and 37 exist in the dataframe
  if(all(c(36, 37) %in% mods_df$position)) {
    # Calculate the sum of the masses for positions 36 and 37
    sum_mass <- sum(mods_df$variable[mods_df$position %in% c(36, 37)], na.rm = TRUE)
    
    # Remove the original rows for position 37
    mods_df <- mods_df[!mods_df$position %in% c(37), ]
    
    # Update the mass at position 36 to the new sum
    if(any(mods_df$position == 36)) {
      mods_df <- mods_df %>%
        mutate(variable = if_else(position == 36, sum_mass, variable))
      # Change the position of this updated mass from 36 to 37
      mods_df <- mods_df %>%
        mutate(position = if_else(position == 36, 36, position))
    } else {
      # If position 36 doesn't exist, just add a new row for position 37
      mods_df <- bind_rows(mods_df, tibble(position = 36, variable = sum_mass))
    }
  }
  
  return(mods_df)
}
