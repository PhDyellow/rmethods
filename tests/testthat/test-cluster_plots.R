context("test-cluster_plots")


library(mclust)

# 4 clusters, spherical with different variances and means
test_plot_dir <- "test_plots"
dir.create(test_plot_dir, showWarnings = FALSE)

set.seed(20190322)
test_size <- 100

test_data <- rbind(MASS::mvrnorm(n = test_size, mu = c(1,1), Sigma = diag(2)*.2),
                   MASS::mvrnorm(n = test_size, mu = c(1,-1), Sigma = diag(2)*.6),
                   MASS::mvrnorm(n = test_size, mu = c(-1,-1), Sigma = diag(2)*.2),
                   MASS::mvrnorm(n = test_size, mu = c(-1,1), Sigma = diag(2)*.1))
test_data <- as.data.frame(test_data)
names(test_data) <- c("x", "y")

clust <- mclust::Mclust(data = test_data, G = 4, modelNames = "VVI")

#plot(clust, what = "classification")
test_that("ellipse plots work", {
  expect_s3_class(plot_clusters(dataset = test_data, col_x = "x", col_y ="y",
                                cluster_model = clust, level = 0.683, legend_thres = 10, alpha=0.3),
                  "ggplot")
  set.seed(1000)

  plot_clusters_test <-plot_clusters(dataset = test_data, col_x = "x", col_y ="y",
                                     cluster_model = clust, level = 0.683, legend_thres = 10, alpha=0.3)
  ggplot2::ggsave(filename = "tmp.png",
                  plot = plot_clusters_test,
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )
  unlink(paste0(test_plot_dir, "/tmp.png"))
  current_hash <- substring(digest::digest(plot_clusters_test), 1, 10)
  ggplot2::ggsave(filename = paste0("plot_clusters_test_",
                                    format(Sys.time(), "%Y-%m-%d"),
                                    "_",
                                    current_hash,
                                    ".png"),
                  plot = plot_clusters_test,
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )
  expect_known_hash(plot_clusters_test, hash = "ff3ac3f61b")
})

set.seed(20190322)
test_size <- 100

test_data <- rbind(MASS::mvrnorm(n = test_size, mu = c(1,1,1,1), Sigma = diag(4)*.2),
                   MASS::mvrnorm(n = test_size, mu = c(1,1,1,-1), Sigma = diag(4)*.6),
                   MASS::mvrnorm(n = test_size, mu = c(1,1,-1,-1), Sigma = diag(4)*.2),
                   MASS::mvrnorm(n = test_size, mu = c(1,1,-1,1), Sigma = diag(4)*.1))
test_data <- as.data.frame(test_data)
names(test_data) <- c("x", "y", "z", "t")



clust <- mclust::Mclust(data = test_data, G = 4, modelNames = "VVI")

test_that("PCA projection of both data and ellipses is working", {
  expect_s3_class(plot_clusters_project(dataset = test_data, transform = "pca",
                                cluster_model = clust, level = 0.683, legend_thres = 10, alpha=0.3),
                  "ggplot")
  set.seed(1000)

  plot_clusters_test <-plot_clusters_project(dataset = test_data, transform = "pca",
                                     cluster_model = clust, level = 0.683, legend_thres = 10, alpha=0.3)
  ggplot2::ggsave(filename = "tmp.png",
                  plot = plot_clusters_test,
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )
  unlink(paste0(test_plot_dir, "/tmp.png"))
  current_hash <- substring(digest::digest(plot_clusters_test), 1, 10)
  ggplot2::ggsave(filename = paste0("plot_clusters_project_test_",
                                    format(Sys.time(), "%Y-%m-%d"),
                                    "_",
                                    current_hash,
                                    ".png"),
                  plot = plot_clusters_test,
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )
  expect_known_hash(plot_clusters_test, hash = "026d6681b8")


})

global_map <- sf::st_read("/vmshare/phd/data/World_EEZ_v8_20140228", layer = "World_EEZ_v8_2014_HR")
country <- "Australia"
bounding_poly  <- global_map[global_map$Country == country,]

library(archivist)
env_trans_hash <- "/vmshare/phd/projects/aus_bioregions/experiments/2019-03-21-1330/archivist/119074835a617f36e24e02d5022dea6b" # in "path/hash" format
env_trans <- archivist::loadFromLocalRepo(repoDir = dirname(env_trans_hash), md5hash = basename(env_trans_hash), value = TRUE)

set.seed(20190322)

test_size <- nrow(env_trans)
n_groups <- 4

group_counts <- split(1:test_size, sort(1:test_size%%n_groups))
group_counts <- sapply(group_counts, length)

test_data <- foreach::foreach(g = group_counts, i = 1:n_groups, .combine = rbind) %do% {
  ret <- as.data.frame(MASS::mvrnorm(n = g, mu = c(i,i)*10, Sigma = diag(2)*i))
  names(ret) <- c("a", "b")
  return(ret)
}

test_data_spatial <- cbind(env_trans[, c("x", "y")], test_data)

clust <- mclust::Mclust(data = test_data, G = n_groups, modelNames = "VVI")
plot(clust, what = "classification")
test_that("Clusters maps plot properly", {
  expect_s3_class(plot_cluster_map(dataset = test_data_spatial, col_x = "x", col_y ="y",
                                   cluster_model = clust, legend_thres = 10, sf_poly= bounding_poly), "ggplot")
  set.seed(1000)
  plot_clusters_map_test <- plot_cluster_map(dataset = test_data_spatial, col_x = "x", col_y ="y",
                                             cluster_model = clust, legend_thres = 10, sf_poly= bounding_poly)
  ggplot2::ggsave(filename = "tmp.png",
                  plot = plot_clusters_map_test,
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )
  unlink(paste0(test_plot_dir, "/tmp.png"))
  current_hash <- substring(digest::digest(plot_clusters_map_test), 1, 10)
  ggplot2::ggsave(filename = paste0("plot_clusters_map_test_",
                                    format(Sys.time(), "%Y-%m-%d"),
                                    "_",
                                    current_hash,
                                    ".png"),
                  plot = plot_clusters_map_test,
                  path = test_plot_dir,
                  device = "png",
                  width = 21,
                  height = 15,
                  units = "cm"

  )

  #Buffer size failure
  expect_known_hash(plot_cluster_map(dataset = test_data_spatial, col_x = "x", col_y ="y",
                                     cluster_model = clust, legend_thres = 10, sf_poly= bounding_poly), hash = "6b5e8dc891")

})


set.seed(20190322)

test_size <- nrow(env_trans)
n_groups <- 4
n_vars <- 6
group_counts <- split(1:test_size, sort(1:test_size%%n_groups))
group_counts <- sapply(group_counts, length)

test_data <- foreach::foreach(g = group_counts, i = 1:n_groups, .combine = rbind) %do% {
  ret <- as.data.frame(MASS::mvrnorm(n = g, mu = rep(i, times = n_vars)*10, Sigma = diag(n_vars)*i))
  names(ret) <- letters[1:n_vars]
  return(ret)
}

test_data_spatial <- cbind(env_trans[, c("x", "y")], test_data)

clust <- mclust::Mclust(data = test_data, G = n_groups, modelNames = "VVI")

test_that("Test pair plotting", {

  expect_length(plot_cluster_pairs(test_data, clust,
                                   level = 0.683, legend_thres = 10, alpha=0.3),
                n_vars*(n_vars-1)/2)
  expect_type(plot_cluster_pairs(test_data, clust,
                                 level = 0.683, legend_thres = 10, alpha=0.3),
              "list")
  expect_s3_class(plot_cluster_pairs(test_data, clust,
                                     level = 0.683, legend_thres = 10, alpha=0.3)[[1]],
                  "ggplot")
  set.seed(1000)
  pairs_test <- plot_cluster_pairs(test_data, clust, level = 0.683, legend_thres = 10, alpha=0.3)
  current_hash <- substring(digest::digest(pairs_test), 1, 10)
  ggplot2::ggsave(filename = paste0("cluster_pairs_test_",
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
  expect_known_hash(pairs_test, hash = "ae9f7f3196")
  sub_vars <- 3
  expect_length(plot_cluster_pairs(test_data, clust, plot_vars = letters[1:sub_vars],
                                   level = 0.683, legend_thres = 10, alpha=0.3),
                sub_vars*(sub_vars-1)/2)
  expect_type(plot_cluster_pairs(test_data, clust, plot_vars = letters[1:sub_vars],
                                 level = 0.683, legend_thres = 10, alpha=0.3),
              "list")
  expect_s3_class(plot_cluster_pairs(test_data, clust, plot_vars = letters[1:sub_vars],
                                     level = 0.683, legend_thres = 10, alpha=0.3)[[1]],
                  "ggplot")
})
