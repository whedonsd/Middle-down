#' Merge hit data frames while preserving character columns
#'
#' @param appended_df Data frame to append.
#' @param main_df Existing data frame to extend.
#'
#' @return Combined data frame with consistent column classes.
#' @export
conditional_merge_hit_df <- function(appended_df, main_df) {
  if (nrow(appended_df) == 0) {
    return(main_df)
  }

  character_cols <- c("a_ions_shared", "b_ions_shared", "y_ions_shared", "c_ions_shared", "z_ions_shared")

  for (col in intersect(character_cols, names(appended_df))) {
    appended_df[[col]] <- as.character(appended_df[[col]])
  }

  for (col in intersect(character_cols, names(main_df))) {
    main_df[[col]] <- as.character(main_df[[col]])
  }

  combined <- rbind(main_df, appended_df)
  rownames(combined) <- NULL
  combined
}
