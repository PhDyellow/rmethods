context("test-data_clean")

#libraries go here
library(archivist)


test_that("convert_csiro_cpr_copepod() properly rejects bad outFormats", {
  expect_error(convert_csiro_cpr_copepod(outFormat = "bad"), "Output format not recognised")
})


local_test_dir <- "local_test"
local_repo <- "local_test_repo"
not_a_dir <- "not_a_dir"
dir.create(local_test_dir)
dir.create(local_repo)
archivist::createLocalRepo(local_repo)

test_that("convert_csiro_cpr_copepod() handles archivist repo folders properly", {
  expect_error(convert_csiro_cpr_copepod(outDir=not_a_dir), "There is no such repository a")
  expect_error(convert_csiro_cpr_copepod(outDir=local_test_dir), "is not a proper repository. There is neither backpack.db nor gallery")
  expect_equal(convert_csiro_cpr_copepod(outDir=local_repo), "e71963322b41ed50f50ec45978dbbeac")
})
unlink(local_test_dir, recursive=TRUE)
unlink(local_repo, recursive=TRUE)
