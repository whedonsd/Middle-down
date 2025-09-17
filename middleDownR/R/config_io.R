#' Configuration helpers
#'
#' Functions that read experiment configurations and derive standardised file-name suffixes used at the beginning of the workflow.
#' @keywords internal
NULL

#' Load configuration file
#'
#' Reads a YAML or JSON configuration file and returns a list of values.
#'
#' @param path Path to the configuration file.
#'
#' @return Named list of configuration values.
#' @examples
#' \dontrun{
#' cfg <- load_config("config.yml")
#' }
#' @export
load_config <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("yml", "yaml")) {
    yaml::read_yaml(path)
  } else if (ext == "json") {
    jsonlite::fromJSON(path)
  } else {
    stop("Unsupported configuration file format: ", ext)
  }
}

#' Construct an output filename from stored parameters
#'
#' @param stored_parameters Named list containing relevant parameters.
#'
#' @return Character string representing the formatted file name.
#' @export
save_name <- function(stored_parameters) {
  paste0(
    "Byonic_",
    "max_mods_", stored_parameters$max_mods, "_",
    "PEP.2D_", stored_parameters$PEP, "_",
    "min_delta_mod_", stored_parameters$min_delta_mod, "_",
    "ppm_tolerance_", stored_parameters$ppm_tolerance, "_",
    stored_parameters$aa_tag, "-tag"
  )
}
