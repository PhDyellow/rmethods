context("test-data_clean")

#libraries go here


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
