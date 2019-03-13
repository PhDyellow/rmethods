#' Convert CSIRO CPR copepod CSV into R format
#'
#' Takes the CSIRO CPR copepod data in CSV form and converts it into
#' a form suitable for processing with Gradient Forest, then
#' exports to different formats. The species list is exported
#' to the same format.
#'
#' The magic numbers in this function are:
#' column indexes and names for non-Species columns
#'
#'
#' @param in_file Pathe to the original CSV file. Default "/vmshare/phd/data/CSIRO_cpr/cpr_copepods_matrix.csv"
#' @param date_format lubridate::parse_date_time() orders string for converting dates into standard R format. Defaults to "ymd HMS"
#'
#' @return A list containing:
#'   \item{data}{A site by species data_frame. Metadata columns "lat", "lon" and "date",
#'   and one column per species.}
#'   \item{sp_names}{Names of species columns in \code{data}}

convert_csiro_cpr_copepod <- function(in_file ="/vmshare/phd/data/CSIRO_cpr/cpr_copepods_matrix.csv",
                                      date_format = "ymd HMS"){
  raw_data <- tryCatch({
    data.table::fread(file= in_file, header = TRUE)
  }, warning = function(w) {
    stop("fread warning: ", w)
  })
  mod_data <- as.data.frame(raw_data)
  names(mod_data) <- make.names(names(raw_data))
  sp_names <- names(mod_data)[-(1:4)]
  names(mod_data)[1:4] <- c("row", "lat", "lon", "date")
  mod_data["row"] <- NULL
  mod_data$date <- tryCatch({
    lubridate::parse_date_time(mod_data$date, orders = date_format)
  }, warning = function(w) {
    stop("Date parse string does not match data (", date_format,  "). ", "Warn: ", w)
  })

  result <- list()
  result$data <- mod_data
  result$sp_names <- sp_names
  return(result)
}
