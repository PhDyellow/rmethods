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


#' Import BioORACLE layers from sdmpredictors library
#'
#' Given a list of environmental variables and a list of
#' statistical measures, returns all valid combinations
#' of environmental variable and statistical measure.
#'
#' Currently only sea surface measures are returned. ("ss" only, not "bdmean" "bdmax" or "bdmin")
#'
#' Available environmental variables are:
#' "depth", "temp", "nitrate", "silicate", "chlo", "iron",  "salinity", "curvel" (current velocity),
#' "icethick", "icecover" "dissox" (dissolved oxygen), "phosphate", "lightbot" (light at bottom),
#' "carbonphyto" (carbon phytoplankton), "pp" (primary production)
#'
#' Available statistical measures are:
#' "mean", "range", "min" (average annual min), "max", "ltmin" (long term minimum), "ltmax",
#'
#' @param env_vars character list of environmental variables. Default is
#' c("depth", "temp", "nitrate", "silicate", "chlo", "iron",  "salinity", "curvel")
#' @param env_modes character list of environmental statistics. Default is c("mean", "range")
#' @param data_dir location for sdm_predictors to save downloaded layers into. Defaults to "./data"
#'
#' @return raster brick of all selected environmental layers
import_biooracle_env <- function(env_vars = c("depth", "temp", "nitrate",
                                              "silicate", "chlo", "iron",
                                              "salinity", "curvel"),
                                 env_modes =c("mean", "range"),
                                 data_dir = "./data"){
  pairs <- data.table::as.data.table(data.table::merge(env_vars, env_modes, all=TRUE))
  bioOracle_names <- apply(pairs[pairs["x"] != "depth",], 1,  function(x){sprintf("BO2_%s%s_ss", x[1], x[2])})
  ## Add bathymetry separately
  if ("depth" %in% env_vars){
    bioOracle_names <- c(bioOracle_names, "MS_bathy_5m")
  }


  env_data <- sdmpredictors::load_layers(bioOracle_names, datadir = data_dir)
  return(env_data)
}
