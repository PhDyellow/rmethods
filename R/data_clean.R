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
#' @return A list containing: \describe{
#'   \item{\code{data}}{A site by species data_frame. Metadata columns "lat", "lon" and "date",
#'   and one column per species.}
#'   \item{\code{sp_names}}{Names of species columns in \code{data}}
#'  }

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
#' Only samples resolved down to the species level are included.
#'
#'
#' @param data_dir Directory to cache SeaAroundUs data
#' @param sau_groups Character Vector of SeaAroundUs groups. Defaults to pelagic groups
#' @param extra_sau_sp Character Vector of additional SeaAroundUs species. Defaults to pelagic sharks.
#' @param bouning_wkt Well Known Text string representation of a polygon or multipolygon. Calling st_as_text() on an sf
#' object will produce valid WKT.
#' @param chunk_size Integer indicating how many cells to fetch in a single call. Too small will slow down the download process, too large will cause server timeouts.
#' @param years Integer range of years from 1950 to 2014. Defaults to 2010
#'
import_seaaroundus_bio <- function(data_dir = "/vmshare/phd/data/SeaAroundUs/ausEEZcatches.rds",
                                   sau_groups = c("pelagiclg", "pelagicmd", "pelagicsm"),
                                   extra_sau_sp = c("Prionace glauca", "Carcharhinus longimanus",
                                                        "Sphyrna lewini", "Sphyrna zygaena", "Alopias vulpinus",
                                                        "Alopias superciliosus", "Lamna nasus",
                                                        "Emmelichthys nitidus nitidus"),
                                   bounding_wkt = NULL,
                                   chunk_size = 100,
                                   years = 2010){

  raw_data <- download_seaaroundus_bio(data_dir = data_dir, bounding_wkt = bounding_wkt,
                                       chunk_size = chunk_size, years = years)

  sp_group_names <- unique(raw_data[raw_data$functional_group_name %in% sau_groups & grepl("^\\S+\\s{1}\\S+$", raw_data$taxon_scientific_name),"taxon_scientific_name"])
  sp_total <- c(sp_group_names, extra_sau_sp)

  site_sp <- reshape2::dcast(
    raw_data[raw_data$taxon_scientific_name %in% sp_total & raw_data$year %in% years,
             c("cell_id", "taxon_scientific_name", "catch_sum")],
    formula = cell_id ~ taxon_scientific_name, value.var = "catch_sum",
    fun.aggregate = mean, fill = 0)
  site_sp_lat <- cbind(sau_id_to_lat_lon(site_sp$cell_id), site_sp)
  site_sp_lat$cell_id <- NULL
  valid_sp_names <- make.names(names(site_sp_lat))
  names(site_sp_lat) <- valid_sp_names
  result <- list()
  result$data <- site_sp_lat
  result$sp_names <- valid_sp_names
  return(result)
}

#' helper function to get SeaAroundUs cells
#'
#' either download or fetch from cache
#'
#'
#'
#' @param data_dir Directory to cache SeaAroundUs data
#' @param bouning_wkt Well Known Text string representation of a polygon or multipolygon. Calling st_as_text() on an sf
#' object will produce valid WKT.
#' @param chunk_size Integer indicating how many cells to fetch in a single call. Too small will slow down the download process, too large will cause server timeouts.
#' @param years Integer range of years from 1950 to 2014. Defaults to 2010
#'
#' @return SeaAroundUs formatted data.frame
#'
#' @importFrom foreach foreach %:% %dopar% %do%
download_seaaroundus_bio <- function(data_dir = "/vmshare/phd/data/SeaAroundUs/cache",
                                     bounding_wkt = NULL,
                                     chunk_size = 100,
                                     years = 2010){
  if(is.null(bounding_wkt)){
    stop("Please provide a WKT string as a boundary box")
  }

  if(dir.exists(data_dir) == FALSE){
    stop("SeaAroundUs cache dir does not exist.")

  }
  cells <- seaaroundus::getcells(shape = bounding_wkt, check_wkt = TRUE)

  cell_chunks <- split(cells, f = ceiling(seq_along(cells)/chunk_size))

  #fetch only first year, to identify and exclude land cells.

  sau_base <- foreach(year = years[1], .combine="rbind", .inorder = FALSE) %:%
    foreach(chunk=cell_chunks, .combine="rbind", .packages = c("seaaroundus"), .inorder = FALSE) %do% {
      ch_hash <- digest::digest(list(chunk, year))
      if(file.exists(sprintf("%s/chunk_%s.rds", data_dir, ch_hash)) == FALSE){
        message(sprintf("Downloading year %d, cells %s.", year, paste(as.character(chunk),collapse = ", ")))

        ##Don't crash the servers
        Sys.sleep(runif(1)/4)

        cell_data <- seaaroundus::getcelldata(year = year, cells = chunk)
        message("saving cells")
        saveRDS(cell_data, sprintf("%s/chunk_%s.rds", data_dir, ch_hash))
      } else {
        message(sprintf("Cache exists for Year %d, cells %s", year, paste(as.character(chunk),collapse = ", ")))

        cell_data <- readRDS(sprintf("%s/chunk_%s.rds", data_dir, ch_hash))
      }

      return(cell_data)

    }



  all_valid <- unique(sau_base$cell_id)

  cell_chunks <- split(all_valid, f = ceiling(seq_along(all_valid)/chunk_size))

  if (length(years) > 1){
    sau_all <- foreach(year = years[-1], .combine="rbind", .inorder = FALSE) %:%
      foreach(chunk=cell_chunks, .combine="rbind", .packages = c("seaaroundus"), .inorder = FALSE) %do% {
        ch_hash <- digest::digest(list(chunk, year))
        if(file.exists(sprintf("%s/chunk_%s.rds", data_dir, ch_hash)) == FALSE){
          message(sprintf("Downloading year %d, cells %s", year, paste(as.character(chunk),collapse = ", ")))

          ##Don't crash the servers
          Sys.sleep(runif(1)/4)

          cell_data <- seaaroundus::getcelldata(year = year, cells = chunk)
          message("saving cells")
          saveRDS(cell_data, sprintf("%s/chunk_%s.rds", data_dir, ch_hash))
        } else {
          message(sprintf("Cache exists for Year %d, cells %s", year, paste(as.character(chunk),collapse = ", ")))

          cell_data <- readRDS(sprintf("%s/chunk_%s.rds", data_dir, ch_hash))
        }

        return(cell_data)

      }

  } else {
    sau_all <- sau_base
  }

  return(sau_all)

}


#' helper function to convert SeaAroundUs IDs to lat and lon values
#' @export
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


#' Take log of selected env vars
#'
#' log_env_data safely applies a log transform to
#' specific columns in an environmental data.frame
#'
#' Generally, I want to transform everything except a few key cols
#' so log_env_data just does that, inverting and offsetting as needed
#'
#' @param dataset A data.frame to transform
#' @param exclude_cols A character vector of col names that will not be transformed
log_env_data <- function(dataset,
                         exclude_cols = NULL
){

  data_names <- names(dataset)

  log_names <- setdiff(data_names, exclude_cols)

  trans_data <- dataset

  for(i in log_names){
    if(all(dataset[[i]] < 0 )) {
      #invert then log
      trans_data[[i]] <- log(abs(dataset[[i]]))
    } else if (any(dataset[[i]] <= 0, na.rm = TRUE )){
      #shift above 0, by 10% of smallest value
      min_val <- min(dataset[[i]])
      smallest <- min(abs(dataset[[i]]))
      if(smallest == 0 ){
        smallest <- .Machine$double.eps
      } else {
        smallest <- smallest/10
      }
      trans_data[[i]] <- log(dataset[[i]] - min_val + smallest)
    } else {
      #Data is all positive
      trans_data[[i]] <- log(dataset[[i]])
    }

  }

  return(trans_data)

}

#' Detect outliers in environmental data.frame
#'
#' Returns a boolean vector of length \code{nrow(dataset)} showing which
#' observations to keep
#'
#' Outliers are specified as being more than \code{range}*\eqn{\sigma} from the mean
#' within a column. Any row with at least one outlier is excluded.
#'
#' @param dataset A data.frame of site by environmental obervations
#' @param range Keep points that lie within \code{range}*\eqn{\sigma} of the mean
#' @param exclude_cols Columns to ignore for identifying rows with outliers
outlier_rows_env <- function(dataset, range = 3, exclude_cols = NULL){

  if(class(dataset) != "data.frame"){
    stop("outlier_rows_env requires a data.frame")
  }

  data_names <- names(dataset)

  filter_names <- setdiff(data_names, exclude_cols)

  outlier_bool <- list()

  for(i in filter_names){
    outlier_bool[[i]] <- flag_upper_outliers(x = dataset[[i]], range = range) |
      flag_lower_outliers(x = dataset[[i]], range = range)
  }
  outlier_bool <- as.data.frame(outlier_bool)
  return(apply(outlier_bool, 1 ,any))

}


#' Remove outliers from a species data.frame
#'
#' Returns a boolean vector of length \code{nrow(dataset)} showing which
#' observations to keep
#'
#' Outliers are specified as being more than \code{range}*\eqn{\sigma} from the mean
#' within a column. Any row with at least one outlier is excluded.
#'
#' Species observations are zero-inflated, so zero's are ignored for calculating
#' the mean and standard deviation. Also, only very large counts above the mean
#' are excluded as outliers, small counts below the mean are always kept.
#'
#' @param dataset A data.frame of site by environmental obervations
#' @param range Keep points that lie within \code{range}*\eqn{\sigma} of the mean
#' @param exclude_cols Columns to ignore for identifying rows with outliers
#'
outlier_rows_sp <- function(dataset, range = 3, exclude_cols = NULL){

  if(class(dataset) != "data.frame"){
    stop("outlier_rows_sp requires a data.frame")
  }

  data_names <- names(dataset)

  filter_names <- setdiff(data_names, exclude_cols)

  outlier_bool <- list()

  for(i in filter_names){
    m <- mean(dataset[[i]][dataset[[i]] > 0])
    s <- sd(dataset[[i]][dataset[[i]] > 0])
    outlier_bool[[i]] <- flag_upper_outliers(x = dataset[[i]], range = range, m = m, s = s)
  }
  outlier_bool <- as.data.frame(outlier_bool)
  return(apply(outlier_bool, 1 ,any))

}



#' Flag upper outliers from a single vector
#'
#' helper function to flag upper outliers from a vector
flag_upper_outliers <- function(x, range, m = NULL, s = NULL){
  if(is.null(m)){
    m <- mean(x)
  }
  if(is.null(s)){
    s <- sd(x)
  }

  outliers <- x > m + range * s
  return(outliers)
}


#' Flag lower outliers from a single vector
#'
#' helper function to flag upper outliers from a vector
#'
flag_lower_outliers <- function(x, range, m = NULL, s = NULL){
  if(is.null(m)){
    m <- mean(x)
  }
  if(is.null(s)){
    s <- sd(x)
  }
  outliers <- x < m - range * s
  return(outliers)
}


#' Filter Rare species
#'
#' filters out species with fewer than n non-zero observations
#'
#'
#' @param dataset A data.frame to filter
#' @param n minimum number of non-zero observations to accept species
#' @param exclude_cols Columns to ignore during filtering, always return true
#'
#' @return boolean vector of columns with rare species
rare_sp_cols <- function(dataset, n = 5, exclude_cols = NULL){

  data_names <- names(dataset)

  filter_names <- setdiff(data_names, exclude_cols)

  always_keep <- data_names %in% exclude_cols

  rare_cols <- as.vector(apply(dataset > 0, 2, sum)) < n

  return(rare_cols & !always_keep)

}

#' Restrict biological data to a bounding polygon
#'
restrict_poly_sp <- function(){
  stop("not implemented yet")

  tempA <- st_as_sf(x = bioData[, .(lon, lat)],
                    coords = c("lon", "lat"),
                    crs = st_crs(boundingPoly))

  withinBoundingPoly <- as.vector(st_contains(boundingPoly,  tempA, sparse = FALSE))


  bioData <- bioData[withinBoundingPoly,]

}

#' Filter dataset for use in Gradient Forest
#'
#' Shapes and filters environmental and biological data
#' into a format suitalbe for passing into a gradient forest
#' fitting method
#'
#' If the env_data is in a data.frame format, then the coord_cols (or rows, if `coord_cols_env = c("")`)must match exactly between bio_data and env_data. Raster* objects assign bio_data rows to the nearest grid cell centre in env_data.
#'
#' @param env_data A Raster* object OR a data.frame containing environmental variables. Format should be site by environmental variables
#' @param bio_data A data.frame of site by species observations of biological data
#' @param coord_cols_bio A character vector of column names in bio_data that define longitude and latitude. When env_data is a raster, longitude must come first.
#' @param coord_cols_env A character vector of column names in env_data that define sites spatially. If env_data is a Raster* object, this parameter is ignored. If coord_cols_env = c(""), then the sites in env_data are assumed to line up with the sites in bio_data
#' @param env_range Number of standard deviations from the mean before a value is considered an outlier.
#' @param bio_range Number of standard deviations from the mean before a value is considered an outlier.
#' @param env_accept_outliers Environmental columns to accept all outliers
#' @param bio_accept_outliers Biological sample columns to accept all outliers
#' @param min_sp_samples Minimum number of non-zero observations needed for species to be accepted
#'
#' @return A list containing: \describe{
#'   \item{\code{data}}{A site by species and environment data.frame}
#'   \item{\code{sp_names}}{Names of species columns in \code{data}}
#'   \item{\code{env_names}}{Names of environmental columns in \code{data}}
#'   \item{\code{coord_names}}{Names of coordinate columns in \code{data}}
#'   }
clean_env_bio_gf <- function(env_data, bio_data,
                             coord_cols_bio, coord_cols_env = NULL,
                             env_range = 3, env_accept_outliers = NULL,
                             bio_range = 3, bio_accept_outliers = NULL,
                             min_sp_samples = 5){

  bio_cols <- names(bio_data)[!(names(bio_data) %in% coord_cols_bio)]
  env_cols <- names(env_data)[!(names(env_data) %in% coord_cols_env)]

  if(grepl(pattern = "Raster[LBS].+", x = class(env_data))){
    site_env <-  as.data.frame(raster::extract(env_data, bio_data[, coord_cols_bio]))
    site_sp_env <- cbind(bio_data, site_env)
  } else if(class(env_data) != "data.frame" ) {
    stop("clean_env_bio_gf: Only Raster* objects and data.frames are supported for env_data")
  } else if(is.null(coord_cols_env)){
    stop("clean_env_bio_gf: When env_data is a data.frame, coord_cols_env must be defined")
  } else if(all(c(coord_cols_env) == "")) {
    ##rows in env_data are assumed to line up with rows in bio_data
    site_sp_env <- cbind(bio_data, env_data)
  } else {
    ##Find sites in env_data that also appear in bio_data
    env_cols <- names(env_data)
    bio_cols <- names(bio_data)
    site_sp_env <- merge(x = env_data, y = bio_data,  by.x = coord_cols_env, by.y = coord_cols_bio)
  }


  site_sp_env_complete <- complete.cases(site_sp_env)

  ## env_cols already excludes the coord_cols_env
  env_outliers <- rmethods::outlier_rows_env(site_sp_env[site_sp_env_complete, env_cols], range = env_range, exclude_cols = env_accept_outliers)

  ## bio_cols already excludes the coord_cols_bio
  sp_outliers <- rmethods::outlier_rows_sp(site_sp_env[site_sp_env_complete, bio_cols][!env_outliers,], range = bio_range, exclude_cols = bio_accept_outliers)

  sp_rare <- rmethods::rare_sp_cols(site_sp_env[site_sp_env_complete, ]
                                    [!env_outliers, ][!sp_outliers,],
                                    n = min_sp_samples,
                                    exclude_cols =  c(env_cols, coord_cols_env, coord_cols_bio))

  sp_rare_cols <- names(site_sp_env[sp_rare])

  ## Logging can be done afterwards, by passing the merged frame into
  ## log_env_data and excluding the coord_cols and the sp_names cols
  ## site_env_logged <- log_env_data(site_env, exclude_cols = env_not_log)
  result <- list()
  result$data <- site_sp_env[site_sp_env_complete, !sp_rare][!env_outliers,][!sp_outliers, ]
  result$sp_names <- bio_cols[!(bio_cols %in% sp_names_rare)]
  result$env_names <- env_cols
  result$coord_names <- names(site_sp_env)[!names(site_sp_env) %in% c(result$sp_names, result$env_names, sp_rare_cols)]
  return(result)
}

