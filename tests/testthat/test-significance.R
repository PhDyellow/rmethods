context("test-significance")

#Testing R.block.multi

test_env <- data.frame(x = c(0, 5, 10), y = c(0,5,10))

test_sample <- data.frame(x = c(5, 7.5), y =c(5,7.5))

test_res <-  data.frame(x = c(1,2), y =c(2,1))

test_that("r.block.multi", {
  #intermediate neighbourhood, residuals mixed
  expect_equal(R.block.multi(test_env, test_sample, test_res, prox = 1),
               as.matrix(data.frame(x = c(1.247664, 1.444672, 1.660756),
                          y = c(1.752336, 1.555328, 1.339244))),
                tolerance = 1e-6
  )
  #small neighbourhood, so weight given to nearest neighbour
  expect_equal(R.block.multi(test_env, test_sample, test_res, prox = 0.1),
               as.matrix(data.frame(x = c(1, 1, 2),
                                    y = c(2, 2, 1))),
               tolerance = 1e-6
  )
  #very large neighbourhood, average of all points
  expect_equal(R.block.multi(test_env, test_sample, test_res, prox = 1000),
               as.matrix(data.frame(x = c(1.5, 1.5, 1.5),
                                    y = c(1.5, 1.5, 1.5))),
               tolerance = 1e-6
  )
})

#test hotellings

interpolated_res <- R.block.multi(test_env, test_sample, test_res, prox = 1)

#Calculate each box with hotellings_summary()
test_dist <- matrix(1, nrow = length(test_env$x), ncol = length(test_env$x)) #1 so diagonals are set properly

for( i in 1:nrow(test_env)) {
  for(j in 1:nrow(test_env)) {
    if(i >= j){
      next
    }

    ret <- hotellings_summary(x_bar = test_env[i,],
                              y_bar = test_env[j,],
                              x_var_vec =  interpolated_res[i,],
                              y_var_vec =  interpolated_res[j,],
                              x_n = 100000,
                              y_n = 100000)
    #symmetric matrix, but I want to operate on full rows
    test_dist[i, j] <- ret$pval
    test_dist[j, i] <- ret$pval

  }

}

test_that("hotellings T^2", {
  expect_equivalent(test_dist, hotellings_bulk(test_env, interpolated_res^2),
    tolerance = 1e-6)
})


