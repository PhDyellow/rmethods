context("test-data_plot")

#Libraries here
library(archivist)
test_plot_dir <- "test_plots"
dir.create(test_plot_dir, showWarnings = FALSE)

bio_data <- archivist::loadFromLocalRepo("a38d13707ac343e705a54297842d7157",
                                         value = TRUE,
                                         repoDir = "/vmshare/phd/projects/aus_bioregions/experiments/2019-03-08-1035_data_csiro_copepod/archivist/")




test_that("Test cor plot", {
  expect_s3_class(plot_cor_heatmap(bio_data[-(1:3)]), "ggplot")
  set.seed(1000)

  cor_heatmap_test <-plot_cor_heatmap(bio_data[-(1:3)])
  ggplot2::ggsave(filename = "tmp.png",
                  plot = cor_heatmap_test,
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"
  )
  unlink(paste0(test_plot_dir, "/tmp.png"))
  current_hash <- substring(digest::digest(cor_heatmap_test), 1, 10)
  ggplot2::ggsave(filename = paste0("cor_heatmap_testB_",
                                    format(Sys.time(), "%Y-%m-%d"),
                                    "_",
                                    current_hash,
                                    ".png"),
                  plot = cor_heatmap_test,
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"
                  )
  expect_known_hash(cor_heatmap_test, hash = "4a1bb71c31")

})

test_that("Test marginal plotting", {
  expect_length(plot_marginals_species(bio_data[-(1:3)]), 237)
  expect_type(plot_marginals_species(bio_data[-(1:3)]), "list")
  expect_s3_class(plot_marginals_species(bio_data[-(1:3)])[[1]], "ggplot")
  set.seed(1000)

  marginals_species_test <-plot_marginals_species(bio_data[-(1:3)])
  ggplot2::ggsave(filename = "tmp.png",
                  plot = marginals_species_test[[1]],
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"
  )
  unlink(paste0(test_plot_dir, "/tmp.png"))
  current_hash <- substring(digest::digest(marginals_species_test), 1, 10)
  ggplot2::ggsave(filename = paste0("marginals_species_test_",
                                    format(Sys.time(), "%Y-%m-%d"),
                                    "_",
                                    current_hash,
                                    ".png"),
                  plot = marginals_species_test[[1]],
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )
  expect_known_hash(marginals_species_test, hash = "ffd0658f87")


  expect_length(plot_marginals(bio_data[-(1:3)]), 237)
  expect_type(plot_marginals(bio_data[-(1:3)]), "list")
  expect_s3_class(plot_marginals(bio_data[-(1:3)])[[1]], "ggplot")
  set.seed(1000)
  marginals_test <-plot_marginals(bio_data[-(1:3)])
  current_hash <- substring(digest::digest(marginals_test), 1, 10)
  ggplot2::ggsave(filename = paste0("marginals_test_",
                                    format(Sys.time(), "%Y-%m-%d"),
                                    "_",
                                    current_hash,
                                    ".png"),
                  plot = marginals_test[[1]],
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )
  expect_known_hash(marginals_test, hash = "1ad9c4d788")

})

test_that("Test pair plotting", {

  expect_length(plot_pairs(bio_data[4:8]), 5*(5-1)/2)
  expect_type(plot_pairs(bio_data[4:8]), "list")
  expect_s3_class(plot_pairs(bio_data[4:8])[[1]], "ggplot")
  set.seed(1000)
  pairs_test <- plot_pairs(bio_data[4:8])
  current_hash <- substring(digest::digest(pairs_test), 1, 10)
  ggplot2::ggsave(filename = paste0("pairs_test_",
                                    format(Sys.time(), "%Y-%m-%d"),
                                    "_",
                                    current_hash,
                                    ".png"),
                  plot = pairs_test[[1]],
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )
  expect_known_hash(pairs_test, hash = "85e8327afe")
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
  expect_known_hash(plot_maps_raster(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = NULL)[[1]], hash = "859d4e38a8")

  #expect_error(plot_maps_raster(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly)[[1]],
  #             "cannot allocate vector of size")

  expect_length(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly), 237)
  expect_type(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly), "list")
  expect_s3_class(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly)[[1]], "ggplot")
  set.seed(1000)
  #expect_known_hash(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = bounding_poly)[[1]], hash = "b1f7e6da98")
  expect_known_hash(plot_maps_points(bio_data[-3], lat_col = "lat", lon_col = "lon", sf_poly = NULL)[[1]], hash = "d4e7cc96cc")
})

gf_model <- archivist::loadFromLocalRepo("6c257ac6fde03138579831b3b2abc17f",
                                         value = TRUE,
                                         repoDir = "/vmshare/phd/projects/aus_bioregions/experiments/2019-03-21-1330_model_gf_sau_2010_aus_eez/archivist/")



test_that("Gradient Importance plots work", {


  gf_imp_test <- plot_gf_importance(gf_model)
  current_hash <- substring(digest::digest(gf_imp_test), 1, 10)
  ggplot2::ggsave(filename = paste0("gf_imp_test_",
                                    format(Sys.time(), "%Y-%m-%d"),
                                    "_",
                                    current_hash,
                                    ".png"),
                  plot = gf_imp_test,
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )
  expect_known_hash(gf_imp_test, hash = "716e09f226")
})

