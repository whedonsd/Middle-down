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
