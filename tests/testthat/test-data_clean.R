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

test_data_dir <- "test_import_sau"
dir.create(test_data_dir, showWarnings = FALSE)
#min,min
#min, max,
#max,max,
#max,min
x_lims <- c(-114.45833,  -94.89167 )
y_lims <- c(-13.97500,  18.25000)
#list, because a polygon can contain holes
test_area <- sf::st_polygon(list(rbind(c(x_lims[1], y_lims[1]),
                                  c(x_lims[1], y_lims[2]),
                                  c(x_lims[2], y_lims[2]),
                                  c(x_lims[2], y_lims[1]),
                                  c(x_lims[1], y_lims[1]))))
#Visually validate the test area, covers South America
#mapWorld <- ggplot2::borders("world", colour="gray50", fill="gray50") # create a layer of borders
#ggplot2::ggplot() + mapWorld + ggplot2::geom_sf(data = test_area, fill = NA, color = "red")

test_that("import_seaaroundus_bio returns correctly", {
  expect_named(import_seaaroundus_bio(data_dir = test_data_dir,
                                      sau_groups = c("pelagiclg", "pelagicmd", "pelagicsm"),
                                      extra_sau_sp = c("Prionace glauca", "Carcharhinus longimanus",
                                                       "Sphyrna lewini", "Sphyrna zygaena", "Alopias vulpinus",
                                                       "Alopias superciliosus", "Lamna nasus",
                                                       "Emmelichthys nitidus nitidus"),
                                      bounding_wkt = sf::st_as_text(test_area),
                                      chunk_size = 100,
                                      years = 2010:2012)[["data"]],
               c("lon", "lat", "Acanthocybium.solandri", "Carcharhinus.longimanus",
                 "Cetengraulis.mysticetus", "Coryphaena.hippurus", "Elops.saurus",
                 "Euthynnus.lineatus", "Harengula.jaguana", "Istiompax.indica",
                 "Istiophorus.albicans", "Istiophorus.platypterus", "Joturus.pichardi",
                 "Kajikia.albida", "Kajikia.audax", "Katsuwonus.pelamis", "Makaira.mazara",
                 "Megalops.atlanticus", "Opisthonema.libertate", "Prionace.glauca",
                 "Rachycentron.canadum", "Sarda.sarda", "Sardinops.sagax", "Scomber.japonicus",
                 "Scomberomorus.maculatus", "Selar.crumenophthalmus", "Tetrapturus.angustirostris",
                 "Thunnus.alalunga", "Thunnus.albacares", "Thunnus.obesus", "Thunnus.orientalis",
                 "Thunnus.thynnus", "Xiphias.gladius"))
  expect_type(import_seaaroundus_bio(data_dir = test_data_dir,
                                     sau_groups = c("pelagiclg", "pelagicmd", "pelagicsm"),
                                     extra_sau_sp = c("Prionace glauca", "Carcharhinus longimanus",
                                                      "Sphyrna lewini", "Sphyrna zygaena", "Alopias vulpinus",
                                                      "Alopias superciliosus", "Lamna nasus",
                                                      "Emmelichthys nitidus nitidus"),
                                     bounding_wkt = sf::st_as_text(test_area),
                                     chunk_size = 100,
                                     years = 2010:2012), "list")
  expect_s3_class(import_seaaroundus_bio(data_dir = test_data_dir,
                                     sau_groups = c("pelagiclg", "pelagicmd", "pelagicsm"),
                                     extra_sau_sp = c("Prionace glauca", "Carcharhinus longimanus",
                                                      "Sphyrna lewini", "Sphyrna zygaena", "Alopias vulpinus",
                                                      "Alopias superciliosus", "Lamna nasus",
                                                      "Emmelichthys nitidus nitidus"),
                                     bounding_wkt = sf::st_as_text(test_area),
                                     chunk_size = 100,
                                     years = 2010)[["data"]], "data.frame")
})

test_data <- data.frame(a = c(-1:1), b = c(0.2:3), c = (-5:-3))
test_that("Taking the log of Env data works reliably", {
  expect_equal(log_env_data(dataset = test_data,
                            exclude_cols = "a"),
               data.frame(a = c(-1:1), b = log(0.2:3), c = log(abs(-5:-3))))
  expect_equal(log_env_data(dataset = test_data,
                            exclude_cols = NULL),
               data.frame(a = log((-1:1)--1+.Machine$double.eps), b = log(0.2:3), c = log(abs(-5:-3))))


})


set.seed(201903)
test_data_length <- 1000
test_data <- data.frame(a = rnorm(n = test_data_length),
                        b = rnorm(n = test_data_length),
                        c = rnorm(n = test_data_length))
#test_data <- c(-10:10, -3:3, -3:3)
test_mean <- apply(X = as.matrix(test_data), MARGIN = 2, FUN = mean)
test_std <- apply(X = as.matrix(test_data), MARGIN = 2, FUN = sd)

test_that("Remvoing env outliers behaves as expected", {
  expect_length(outlier_rows_env(dataset = test_data, range = 1), test_data_length)
  #At range = 1 for Gaussian data, 1/3 of samples are outliers
  #Probability of a whole row (3 variables) with no outliers is (2/3)^3 = 0.29
  # Probability of a row being excluded is 1 - 0.29 = 0.71 or 710 out of 1000
  expect_equal(sum(outlier_rows_env(dataset = test_data, range = 1)), 685) #With seed 201903
  #At range = 2 for Gaussian data, 1/20 of samples are outliers
  #Probability of a whole row (3 variables) with no outliers is (19/20)^3 = 0.85
  # Probability of a row being excluded is 1 - 0.85 = 0.15 or 150 out of 1000
  expect_equal(sum(outlier_rows_env(dataset = test_data, range = 2)), 127) # With seed 201903
  expect_error(outlier_rows_env(dataset = as.matrix(test_data), range = 1), "outlier_rows_env requires a data.frame")

})


set.seed(201903)
test_data_length <- 10000
prob_zero <- 0.5
lambda <- 16
test_data_raw <- data.frame(a = rpois(n = test_data_length, lambda = lambda),
                        b = rpois(n = test_data_length, lambda = lambda),
                        c = rpois(n = test_data_length, lambda = lambda))
zeros <- data.frame(a = rbinom(test_data_length, 1, prob = prob_zero),
                    b = rbinom(test_data_length, 1, prob = prob_zero),
                    c = rbinom(test_data_length, 1, prob = prob_zero))
test_data <- test_data_raw*zeros
#test_data <- c(-10:10, -3:3, -3:3)
test_mean <- apply(X = as.matrix(test_data), MARGIN = 2, FUN = mean)
test_std <- apply(X = as.matrix(test_data), MARGIN = 2, FUN = sd)

test_that("Remvoing species outliers behaves as expected", {
  expect_length(outlier_rows_sp(dataset = test_data, range = 1), test_data_length)
  # The data is zero inflated, but 0's are ignored for finding mean and variance
  # Also, small sample counts are never outliers for species data, but may be noise in the signal.
  # Assuming std deviation is approximately similar for Poisson as for Gaussian at lambda = 16
  # At range = 1 for Poisson data (lambda = 16), 1/6 of samples are upper outliers, and
  # 1/2 of all samples are discarded as 0. So 1/12 samples will be outliers.
  # Probability of a whole row (3 variables) with no outliers is (11/12)^3 = 0.77
  # Probability of a row being excluded is 1 - 0.77 = 0.23 or 2300 out of 10000
  expect_equal(sum(outlier_rows_sp(dataset = test_data, range = 1)), 2320) #With seed 201903
  # At range = 2 for Poisson data (lambda = 16), 1/40 of samples are upper outliers, and
  # 1/2 of all samples are discarded as 0. So 1/80 samples will be outliers.
  # Probability of a whole row (3 variables) with no outliers is (79/80)^3 = 0.96
  # Probability of a row being excluded is 1 - 0.96 = 0.04 or 400 out of 10000
  expect_equal(sum(outlier_rows_sp(dataset = test_data, range = 2)), 469) # With seed 201903
  expect_error(outlier_rows_sp(dataset = as.matrix(test_data), range = 1), "outlier_rows_sp requires a data.frame")
})


test_data <- data.frame(a = 1:10, b= rep(0, length.out = length(1:10)), c= rep(c(0,1), length.out = length(1:10), each= 6))
test_that("Excluding rare species columns works as expected",{
  expect_equal(rare_sp_cols(test_data, n = 5, exclude_cols = NULL), c(FALSE, TRUE, TRUE))
  expect_equal(rare_sp_cols(test_data, n = 2, exclude_cols = NULL), c(FALSE, TRUE, FALSE))
  expect_equal(rare_sp_cols(test_data, n = 2, exclude_cols = c("b")), c(FALSE, FALSE, FALSE))
  expect_equal(rare_sp_cols(test_data, n = 20, exclude_cols = c("b")), c(TRUE, FALSE, TRUE))

})
