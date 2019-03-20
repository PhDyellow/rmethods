context("test-data_clean")

#libraries go here

test_data_dir <- "test_import"
dir.create(test_data_dir, showWarnings = FALSE)

test_that("convert_csiro_copepod() returns list of species names and data", {
  expect_named(convert_csiro_cpr_copepod(in_file ="/vmshare/phd/data/CSIRO_cpr/cpr_copepods_matrix.csv",
                                         date_format = "ymd HMS"),
               c("data", "sp_names")
               )
  expect_s3_class(convert_csiro_cpr_copepod(in_file ="/vmshare/phd/data/CSIRO_cpr/cpr_copepods_matrix.csv",
                                        date_format = "ymd HMS")[["data"]],
              "data.frame"
              )

  expect_type(convert_csiro_cpr_copepod(in_file ="/vmshare/phd/data/CSIRO_cpr/cpr_copepods_matrix.csv",
                                        date_format = "ymd HMS")[["sp_names"]],
              "character"
              )


})

test_that("import_seaaroundus_bio returns correctly", {
  expect_error(import_seaaroundus_bio(data_dir = test_data_dir,
                                       sau_groups = c("pelagiclg", "pelagicmd", "pelagicsm"),
                                       extra_sau_sp = c("Prionace glauca", "Carcharhinus longimanus",
                                                        "Sphyrna lewini", "Sphyrna zygaena", "Alopias vulpinus",
                                                        "Alopias superciliosus", "Lamna nasus",
                                                        "Emmelichthys nitidus nitidus"),
                                       years = 2010),
    "import_seaaroundus_bio is still under development, and must have an RDS file")
  expect_type(import_seaaroundus_bio(data_dir = "/vmshare/phd/data/SeaAroundUs/ausEEZcatches.rds",
                                     sau_groups = c("pelagiclg", "pelagicmd", "pelagicsm"),
                                     extra_sau_sp = c("Prionace glauca", "Carcharhinus longimanus",
                                                      "Sphyrna lewini", "Sphyrna zygaena", "Alopias vulpinus",
                                                      "Alopias superciliosus", "Lamna nasus",
                                                      "Emmelichthys nitidus nitidus"),
                                     years = 2010), "list")
  expect_s3_class(import_seaaroundus_bio(data_dir = "/vmshare/phd/data/SeaAroundUs/ausEEZcatches.rds",
                                     sau_groups = c("pelagiclg", "pelagicmd", "pelagicsm"),
                                     extra_sau_sp = c("Prionace glauca", "Carcharhinus longimanus",
                                                      "Sphyrna lewini", "Sphyrna zygaena", "Alopias vulpinus",
                                                      "Alopias superciliosus", "Lamna nasus",
                                                      "Emmelichthys nitidus nitidus"),
                                     years = 2010)[["data"]], "data.frame")
})
