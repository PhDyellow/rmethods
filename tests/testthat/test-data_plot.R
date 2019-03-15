context("test-data_plot")

#Libraries here
library(archivist)

dir.create("test_plots", showWarnings = FALSE)

bio_data <- archivist::loadFromLocalRepo("a38d13707ac343e705a54297842d7157",
                                         value = TRUE,
                                         repoDir = "/vmshare/phd/projects/aus_bioregions/experiments/2019-03-08-1035/archivist/")




test_that("Test cor plot", {
  expect_s3_class(plot_cor_heatmap(bio_data[-(1:3)]), "ggplot")
  set.seed(1000)
  expect_known_hash(plot_cor_heatmap(bio_data[-(1:3)]), hash = "3f27a53d0c")

})

test_that("Test marginal plotting", {
  expect_length(plot_marginals_species(bio_data[-(1:3)]), 237)
  expect_type(plot_marginals_species(bio_data[-(1:3)]), "list")
  expect_s3_class(plot_marginals_species(bio_data[-(1:3)])[[1]], "ggplot")
  set.seed(1000)
  expect_known_hash(plot_marginals_species(bio_data[-(1:3)]), hash = "f9cfda9fc1")


  expect_length(plot_marginals(bio_data[-(1:3)]), 237)
  expect_type(plot_marginals(bio_data[-(1:3)]), "list")
  expect_s3_class(plot_marginals(bio_data[-(1:3)])[[1]], "ggplot")
  set.seed(1000)
  expect_known_hash(plot_marginals(bio_data[-(1:3)]), hash = "543b6a2152")

})

test_that("Test pair plotting", {

  expect_length(plot_pairs(bio_data[4:8]), 5*(5-1)/2)
  expect_type(plot_pairs(bio_data[4:8]), "list")
  expect_s3_class(plot_pairs(bio_data[4:8])[[1]], "ggplot")
  set.seed(1000)
  expect_known_hash(plot_pairs(bio_data[4:8]), hash = "7e98a7ca99")
})


global_map <- sf::st_read("/vmshare/phd/data/World_EEZ_v8_20140228", layer = "World_EEZ_v8_2014_HR")
country <- "Australia"
bounding_poly  <- global_map[global_map$Country == country,]


test_that("Test map plotting", {
  expect_length(plot_maps_raster(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly), 237)
  expect_type(plot_maps_raster(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly), "list")
  expect_s3_class(plot_maps_raster(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly)[[1]], "ggplot")
  set.seed(1000)

  #Buffer size failure
  expect_known_hash(plot_maps_raster(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = NULL)[[1]], hash = "bb449420aa")

  #expect_error(plot_maps_raster(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly)[[1]],
  #             "cannot allocate vector of size")

  expect_length(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly), 237)
  expect_type(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly), "list")
  expect_s3_class(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly)[[1]], "ggplot")
  set.seed(1000)
  #expect_known_hash(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly)[[1]], hash = "b1f7e6da98")
  expect_known_hash(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = NULL)[[1]], hash = "7ddeeb367e")
})

