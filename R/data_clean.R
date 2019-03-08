#' Convert CSIRO CPR copepod CSV into R format
#'
#' Takes the CSIRO CPR copepod data in CSV form and converts it into
#' a form suitable for processing with Gradient Forest, then
#' exports to different formats. The species list is exported
#' to the same format.
#'
#' The magic numbers in this function are:
#' column indexes and names for non-Species columns
#' The date parsing string
#'
#'
#' @param inFile Pathe to the original CSV file. Default "/vmshare/phd/data/CSIRO_cpr/cpr_copepods_matrix.csv"
#' @param outFormat Type of output. Values are ("csv", "rds", "archivist"). Default "archivist"
#' @param outDir Directory of output file. Defaults to "results"
#' @param outFile Name of output file. Defaults to inFile basename
#' @param outFileSpecies Name of species list. Defaults to paste0(outFile, "_species")
#'
#' @return File name as char string for outFormats "csv" and "rds" and hash string for ("archivist")

convert_csiro_cpr_copepod <- function(inFile ="/vmshare/phd/data/CSIRO_cpr/cpr_copepods_matrix.csv",
                                      outFormat = "archivist", outDir = "results", outFile =NULL,
                                      outFileSpecies = NULL){

  if (!(outFormat[1] %in% c("csv", "rds", "archivist"))) {
    stop('Output format not recognised. Use one of ("csv", "rds", "archivist")')
  }

  if(is.null(outFile)){
    outFile <- unlist(strsplit(basename(inFile), split = '.', fixed = TRUE))[1]
  }

  if(is.null(outFileSpecies)){
    outFileSpecies <- paste0(outFile, "_species")
  }
  
  print(outFileSpecies)
  raw_data <- data.table::fread(file= inFile, header = TRUE)
  mod_data <- as.data.frame(raw_data)
  names(mod_data) <- make.names(names(raw_data))
  sp_names <- names(mod_data)[-(1:4)]
  names(mod_data)[1:4] <- c("row", "lat", "lon", "date")
  mod_data$date <- tryCatch({
    lubridate::parse_date_time(mod_data$date, orders = "ymd HMS")
  }, warning = function(w) {
    stop("Date parse string does not match data")
  }

  )


  result <- list()
  switch(outFormat[1],
        csv = {
           data.table::fwrite(x = mod_data, file = paste0(outDir, outFile, ".csv"))
           data.table::fwrite(x = sp_names, file = paste0(outDir, outFileSpecies, ".csv"))
           result$outFile <- paste0(outDir, outFile, ".csv")
           result$outFileSpecies <- paste0(outDir, outFileSpecies, ".csv")
  },

        rds = {
           saveRDS(object = mod_data, file = paste0(outDir, outFile, ".rds"))
           saveRDS(object = sp_names, file = paste0(outDir, outFileSpecies, ".rds"))
           result$outFile <- paste0(outDir, outFile, ".rds")
           result$outFileSpecies <- paste0(outDir, outFileSpecies, ".rds")
         },

        archivist = {
          result$outHash <- archivist::saveToLocalRepo(artifact = mod_data, repoDir = outDir, artifactName = outFile)
          result$outHashSpecies <- archivist::saveToLocalRepo(artifact = sp_names, repoDir = outDir, artifactName = outFileSpecies)
  },
         {
           stop('Unrecognised output format. Use one of ("csv", "rds", "archivist")')
         }
  )


}
