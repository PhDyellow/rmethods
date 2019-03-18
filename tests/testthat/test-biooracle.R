context("test-biooracle")


dir.create("data", showWarnings = FALSE)

test_that("import_biooracle_env() returns appropriate raster bricks",{
  expect_error(import_biooracle_env(env_vars = as.data.frame( c("temp", "salinity")), env_modes = c("mean", "range"), data_dir = "data"),
               "Parameters `env_vars` and `env_modes` must be character vectors.")
  expect_error(import_biooracle_env(env_vars =  c("temp", "salinity"), env_modes = as.data.frame(c("mean", "range")), data_dir = "data"),
               "Parameters `env_vars` and `env_modes` must be character vectors.")
  expect_s4_class(import_biooracle_env(env_vars = c("temp", "salinity"), env_modes = c("mean", "range"), data_dir = "data"),
                  "RasterStack")
  expect_named(import_biooracle_env(env_vars = c("temp", "salinity"), env_modes = c("mean", "range"), data_dir = "data"),
               c("BO2_tempmean_ss", "BO2_salinitymean_ss", "BO2_temprange_ss",  "BO2_salinityrange_ss"),
               ignore.order = TRUE)
  expect_named(import_biooracle_env(env_vars = c("temp", "salinity", "depth"), env_modes = c("mean", "range"), data_dir = "data"),
               c("BO2_tempmean_ss", "BO2_salinitymean_ss", "BO2_temprange_ss",  "BO2_salinityrange_ss", "MS_bathy_5m"),
               ignore.order = TRUE)


})
