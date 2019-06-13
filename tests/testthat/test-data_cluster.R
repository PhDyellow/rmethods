context("test-data_cluster")

#4 distinct clusters

set.seed(20190612)
test_size <- 100
thres <- 0.1
test_distinct <- rbind(MASS::mvrnorm(n = test_size/2, mu = c(5,5), Sigma = diag(2)*.2),
                   MASS::mvrnorm(n = test_size*2, mu = c(5,-5), Sigma = diag(2)*.6),
                   MASS::mvrnorm(n = test_size, mu = c(-5,-5), Sigma = diag(2)*.2),
                   MASS::mvrnorm(n = test_size, mu = c(-5,5), Sigma = diag(2)*.1))
test_distinct <- as.data.frame(test_distinct)
names(test_distinct) <- c("x", "y")
test_distinct_dist <- fields::rdist(test_distinct)
test_distinct_sim <- 1/(1+test_distinct_dist)



# 4 overlapping clusters

  test_overlap <- rbind(MASS::mvrnorm(n = test_size/2, mu = c(1,1), Sigma = diag(2)*.2),
                     MASS::mvrnorm(n = test_size*2, mu = c(1,-1), Sigma = diag(2)*.6),
                     MASS::mvrnorm(n = test_size, mu = c(-1,-1), Sigma = diag(2)*.2),
                     MASS::mvrnorm(n = test_size, mu = c(-1,1), Sigma = diag(2)*.1))
  test_overlap <- as.data.frame(test_overlap)
names(test_overlap) <- c("x", "y")
test_overlap_dist <- fields::rdist(test_overlap)
test_overlap_sim <- 1/(1+test_overlap_dist)


# 4 Gaussians, nearly indistinguishable

test_blur <- rbind(MASS::mvrnorm(n = test_size/2, mu = c(0.1,0.1), Sigma = diag(2)*.2),
                   MASS::mvrnorm(n = test_size*2, mu = c(0.1,-0.1), Sigma = diag(2)*.6),
                   MASS::mvrnorm(n = test_size, mu = c(-0.1,-0.1), Sigma = diag(2)*.2),
                   MASS::mvrnorm(n = test_size, mu = c(-0.1,0.1), Sigma = diag(2)*.1))
test_blur <- as.data.frame(test_blur)
names(test_blur) <- c("x", "y")
test_blur_dist <- fields::rdist(test_blur)
test_blur_sim <- 1/(1+test_blur_dist)


clust_distinct <- conn_clust_rec(test_distinct_sim > thres, stop_all = TRUE, break_ties = "random", min_sites = 1, max_primary = Inf)
clust_overlap <- conn_clust_rec(test_overlap_sim > thres, stop_all = TRUE, break_ties = "random", min_sites = 1, max_primary = Inf)
clust_blur <- conn_clust_rec(test_blur_sim > thres, stop_all = TRUE, break_ties = "random", min_sites = 1, max_primary = Inf)

membership_distinct <- (test_distinct_sim > thres)[, clust_distinct ]
membership_overlap <- (test_overlap_sim > thres)[, clust_overlap ]
membership_blur <- (test_blur_sim > thres)[, clust_distinct ]

test_that("conn_clust_rec", {
  expect_equal(clust_distinct, c(173, 410, 348,   37))
  expect_equal(clust_overlap, c(265))
  expect_equal(clust_blur, c(396))
})

#I have not tested all the permutations of input options here.
