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
#'
import_biooracle_env <- function(env_vars = c("depth", "temp", "nitrate",
                                              "silicate", "chlo", "iron",
                                              "salinity", "curvel"),
                                 env_modes =c("mean", "range"),
                                 data_dir = "./data"){
  if(class(env_vars) != "character" | class(env_modes) != "character"){
    stop("Parameters `env_vars` and `env_modes` must be character vectors.")
  }

  pairs <- merge(as.data.frame(env_vars, stringsAsFactors = FALSE), as.data.frame(env_modes,stringsAsFactors = FALSE), all=TRUE)
  bioOracle_names <- apply(pairs[pairs["env_vars"] != "depth",], 1,  function(x){sprintf("BO2_%s%s_ss", x[1], x[2])})
  ## Add bathymetry separately
  if ("depth" %in% env_vars){
    bioOracle_names <- c(bioOracle_names, "MS_bathy_5m")
  }


  env_data <- sdmpredictors::load_layers(bioOracle_names, datadir = data_dir)
  names(env_data) <- bioOracle_names
  return(env_data)
}


#' Import SeaAroundUs data
#'
#' Pull in SeaAroundUs data and clean up to site by species matrix
#'
#' For now, just an RDS file, function is not accessing SeaAroundUs servers or a local cache.
#'
#'
#'
#' @param data_dir Directory to cache SeaAroundUs data
#' @param sau_groups Character Vector of SeaAroundUs groups. Defaults to pelagic groups
#' @param extra_sau_sp Character Vector of additional SeaAroundUs species. Defaults to pelagic sharks.
#' @param years Integer range of years from 1950 to 2014. Defaults to 2010
import_seaaroundus_bio <- function(data_dir = "/vmshare/phd/data/SeaAroundUs/ausEEZcatches.rds",
                                   sau_groups = c("pelagiclg", "pelagicmd", "pelagicsm"),
                                   extra_sau_sp = c("Prionace glauca", "Carcharhinus longimanus",
                                                        "Sphyrna lewini", "Sphyrna zygaena", "Alopias vulpinus",
                                                        "Alopias superciliosus", "Lamna nasus",
                                                        "Emmelichthys nitidus nitidus"),
                                   years = 2010){

  if(strsplit(x =basename(data_dir), split = "\\.")[[1]][2] != "rds" | is.na(strsplit(x =basename(data_dir), split = "\\.")[[1]][2])){

    stop("import_seaaroundus_bio is still under development, and must have an RDS file.")
  }

  raw_data <- readRDS(data_dir)
  pelagic_sp_names <- unique(raw_data[raw_data$functional_group_name %in% sau_groups & grepl("^\\S+\\s{1}\\S+$", raw_data$taxon_scientific_name),"taxon_scientific_name"])
  pelagic_sp_total <- c(pelagic_sp_names, extra_sau_sp)

  site_sp <- data.table::dcast(raw_data[raw_data$taxon_scientific_name %in% pelagic_sp_total & raw_data$year %in% years, c("cell_id", "taxon_scientific_name", "catch_sum")], formula = cell_id ~ taxon_scientific_name, value.var = "catch_sum", fun.aggregate = mean)
  site_sp_lat <- cbind(sau_id_to_lat_lon(site_sp$cell_id), site_sp)
  site_sp_lat$cell_id <- NULL
  valid_sp_names <- make.names(names(site_sp_lat))
  names(site_sp_lat) <- valid_sp_names
  result <- list()
  result$data <- site_sp_lat
  result$sp_names <- valid_sp_names
  return(result)
}

#' helper function to convert SeaAroundUs IDs to lat and lon values
sau_id_to_lat_lon <- function(cellId, gridSize = 0.5, centred = TRUE){
  #need to reverse a cycle. every 720 steps is a latitude band
  lonPerLat <- 360/gridSize
  firstLat <- 90 - gridSize
  firstLon <- -180 + gridSize
  cellIdZeroed <- cellId - 1
  latSteps <- floor(cellIdZeroed / lonPerLat)
  lonSteps <- cellIdZeroed - latSteps*lonPerLat

  if (centred){
    data.frame(lon = firstLon + gridSize * lonSteps - (gridSize/2), lat = firstLat - gridSize * latSteps + (gridSize/2))
  } else {
    data.frame(lon = firstLon + gridSize * lonSteps, lat = firstLat - gridSize * latSteps)
  }
}
