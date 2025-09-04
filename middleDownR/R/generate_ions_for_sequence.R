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
  
  for (charge in 1:5) { # Revised to +5 charge state max based on Byonic ion depiction
    a_ions <- numeric(sequence_length - 1)
    a_nl_O_ions <- numeric(sequence_length - 1)
    a_nl_N_ions <- numeric(sequence_length - 1)
    a_nl_Kme3_ions <- numeric(sequence_length - 1) # New vector for Kme3 neutral loss
    a_nl_pST_ions <- numeric(sequence_length - 1) # New vector for pST neutral loss
    a_nl_Ksu_ions <- numeric(sequence_length - 1) # New vector for Ksu neutral loss
    b_ions <- numeric(sequence_length - 1)
    b_nl_O_ions <- numeric(sequence_length - 1)
    b_nl_N_ions <- numeric(sequence_length - 1)
    b_nl_Kme3_ions <- numeric(sequence_length - 1) # New vector for Kme3 neutral loss
    b_nl_pST_ions <- numeric(sequence_length - 1) # New vector for pST neutral loss
    b_nl_Ksu_ions <- numeric(sequence_length - 1) # New vector for Ksu neutral losS
    y_ions <- numeric(sequence_length - 1)
    y_nl_O_ions <- numeric(sequence_length - 1)
    y_nl_N_ions <- numeric(sequence_length - 1)
    y_nl_Kme3_ions <- numeric(sequence_length - 1) # New vector for Kme3 neutral loss
    y_nl_pST_ions <- numeric(sequence_length - 1) # New vector for pST neutral loss
    y_nl_Ksu_ions <- numeric(sequence_length - 1) # New vector for Ksu neutral losS
    c_ions <- numeric(sequence_length - 1)
    c_nl_O_ions <- numeric(sequence_length - 1)
    c_nl_N_ions <- numeric(sequence_length - 1)
    c_nl_Kme3_ions <- numeric(sequence_length - 1) # New vector for Kme3 neutral loss
    c_nl_pST_ions <- numeric(sequence_length - 1) # New vector for pST neutral loss
    c_nl_Ksu_ions <- numeric(sequence_length - 1) # New vector for Ksu neutral losS
    z_ions <- numeric(sequence_length - 1)
    z_nl_O_ions <- numeric(sequence_length - 1)
    z_nl_N_ions <- numeric(sequence_length - 1)
    z_nl_Kme3_ions <- numeric(sequence_length - 1) # New vector for Kme3 neutral loss
    z_nl_pST_ions <- numeric(sequence_length - 1) # New vector for pST neutral loss
    z_nl_Ksu_ions <- numeric(sequence_length - 1) # New vector for Ksu neutral losS
    
    for (i in 1:(sequence_length - 1)) {
      a_seq <- substr(sequence, 1, i)
      b_seq <- substr(sequence, 1, i)
      y_seq <- substr(sequence, sequence_length - i + 1, sequence_length)
      c_seq <- substr(sequence, 1, i)
      z_seq <- substr(sequence, sequence_length - i + 1, sequence_length)
      
      # For b ions, use positions directly as they are from N-terminus
      a_mods <- modifications_df[modifications_df$position <= i,]
      b_mods <- modifications_df[modifications_df$position <= i,]
      c_mods <- modifications_df[modifications_df$position <= i,]
      
      # For y ions, adjust positions relative to the fragment's C-terminus
      y_mods <- modifications_df
      if (nrow(y_mods) > 0) {
        y_mods$position <- sequence_length - y_mods$position + 1
        y_mods <- y_mods[y_mods$position <= i,]
      }
      z_mods <- modifications_df
      if (nrow(z_mods) > 0) {
        z_mods$position <- sequence_length - z_mods$position + 1
        z_mods <- z_mods[z_mods$position <= i,]
      }

      # Check for specific modification by mass and calculate neutral loss
      if (any(abs(a_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Kme3_loss <- 59.073499  # Define this value
        a_mass_specific_mod_loss <- a_mass - specific_Kme3_loss
        a_nl_Kme3_ions[i] <- calculate_ion_mz(a_mass_specific_mod_loss, charge)
      }
      if (any(abs(b_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Kme3_loss <- 59.073499  # Define this value
        b_mass_specific_mod_loss <- b_mass - specific_Kme3_loss
        b_nl_Kme3_ions[i] <- calculate_ion_mz(b_mass_specific_mod_loss, charge)
      }
      if (any(abs(y_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Kme3_loss <- 59.073499  # Define this value
        y_mass_specific_mod_loss <- y_mass - specific_Kme3_loss
        y_nl_Kme3_ions[i] <- calculate_ion_mz(y_mass_specific_mod_loss, charge)
      }
      if (any(abs(c_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Kme3_loss <- 59.073499  # Define this value
        c_mass_specific_mod_loss <- c_mass - specific_Kme3_loss
        c_nl_Kme3_ions[i] <- calculate_ion_mz(c_mass_specific_mod_loss, charge)
      }
      if (any(abs(z_mods$variable - 42.046950) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Kme3_loss <- 59.073499  # Define this value
        z_mass_specific_mod_loss <- z_mass - specific_Kme3_loss
        z_nl_Kme3_ions[i] <- calculate_ion_mz(z_mass_specific_mod_loss, charge)
      }
      if (any(abs(a_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_pST_loss <- 97.976896  # Define this value
        a_mass_specific_mod_loss <- a_mass - specific_pST_loss
        a_nl_pST_ions[i] <- calculate_ion_mz(a_mass_specific_mod_loss, charge)
      }
      if (any(abs(b_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_pST_loss <- 97.976896  # Define this value
        b_mass_specific_mod_loss <- b_mass - specific_pST_loss
        b_nl_pST_ions[i] <- calculate_ion_mz(b_mass_specific_mod_loss, charge)
      }
      if (any(abs(y_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_pST_loss <- 97.976896  # Define this value
        y_mass_specific_mod_loss <- y_mass - specific_pST_loss
        y_nl_pST_ions[i] <- calculate_ion_mz(y_mass_specific_mod_loss, charge)
      }
      if (any(abs(c_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_pST_loss <- 97.976896  # Define this value
        c_mass_specific_mod_loss <- c_mass - specific_pST_loss
        c_nl_pST_ions[i] <- calculate_ion_mz(c_mass_specific_mod_loss, charge)
      }
      if (any(abs(z_mods$variable - 79.966331) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_pST_loss <- 97.976896  # Define this value
        z_mass_specific_mod_loss <- z_mass - specific_pST_loss
        z_nl_pST_ions[i] <- calculate_ion_mz(z_mass_specific_mod_loss, charge)
      }
      if (any(abs(a_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Ksu_loss <- 101.02387  # Define this value
        a_mass_specific_mod_loss <- a_mass - specific_Ksu_loss
        a_nl_Ksu_ions[i] <- calculate_ion_mz(a_mass_specific_mod_loss, charge)
      }
      if (any(abs(b_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Ksu_loss <- 101.02387  # Define this value
        b_mass_specific_mod_loss <- b_mass - specific_Ksu_loss
        b_nl_Ksu_ions[i] <- calculate_ion_mz(b_mass_specific_mod_loss, charge)
      }
      if (any(abs(y_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Ksu_loss <- 101.02387  # Define this value
        y_mass_specific_mod_loss <- y_mass - specific_Ksu_loss
        y_nl_Ksu_ions[i] <- calculate_ion_mz(y_mass_specific_mod_loss, charge)
      }
      if (any(abs(c_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Ksu_loss <- 101.02387  # Define this value
        c_mass_specific_mod_loss <- c_mass - specific_Ksu_loss
        c_nl_Ksu_ions[i] <- calculate_ion_mz(c_mass_specific_mod_loss, charge)
      }
      if (any(abs(z_mods$variable - 100.01604) < .Machine$double.eps^0.5)) {
        # Assuming a specific neutral loss for this modification
        specific_Ksu_loss <- 101.02387  # Define this value
        z_mass_specific_mod_loss <- z_mass - specific_Ksu_loss
        z_nl_Ksu_ions[i] <- calculate_ion_mz(z_mass_specific_mod_loss, charge)
      }
      
      # Calculate the mass of the b and y ion sequences
      a_mass <- calculate_fragment_mass_from_sequence(a_seq, a_mods) + ion_type_adjustments$a
      b_mass <- calculate_fragment_mass_from_sequence(b_seq, b_mods) + ion_type_adjustments$b
      y_mass <- calculate_fragment_mass_from_sequence(y_seq, y_mods) + ion_type_adjustments$y
      c_mass <- calculate_fragment_mass_from_sequence(c_seq, c_mods) + ion_type_adjustments$c
      z_mass <- calculate_fragment_mass_from_sequence(z_seq, z_mods) + ion_type_adjustments$z
      
      # Check for potential neutral losses and calculate accordingly
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
      # Check for potential neutral losses and calculate accordingly
      if (grepl("[KRQN]", y_seq)) {
        y_mass_N_loss <- y_mass - ion_type_adjustments$Ion_AmmoniaLossMass
        y_nl_N_ions[i] <- calculate_ion_mz(y_mass_N_loss, charge)
      }
      if (grepl("[STED]", y_seq)) {
        y_mass_O_loss <- y_mass - ion_type_adjustments$Ion_WaterLossMass
        y_nl_O_ions[i] <- calculate_ion_mz(y_mass_O_loss, charge)
      }
      # Check for potential neutral losses and calculate accordingly
      if (grepl("[KRQN]", c_seq)) {
        c_mass_N_loss <- c_mass - ion_type_adjustments$Ion_AmmoniaLossMass
        c_nl_N_ions[i] <- calculate_ion_mz(c_mass_N_loss, charge)
      }
      if (grepl("[STED]", c_seq)) {
        c_mass_O_loss <- c_mass - ion_type_adjustments$Ion_WaterLossMass
        c_nl_O_ions[i] <- calculate_ion_mz(c_mass_O_loss, charge)
      }
      # Check for potential neutral losses and calculate accordingly
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
  
  return(ions_list)
