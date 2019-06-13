#' Interpolate residuals
#'
#' Given a set of survey sites with residuals and a set of environmental sites,
#' estimate the residual at each environmental site with a
#' Gaussian density kernel.
#'
#' Distances are in transformed predictor space, not latitude or longitude.
#' Sites with similar environmental conditions will have similar residuals
#' regardless of geographic distance.
#'
#'  @param Xblock data frame of site by predictors that need residual estimates
#'  @param Xsite data.frame of site by predictors that have known residuals
#'  @param res data.frame of site by predictor_residuals, must match dimensions of Xsite
#'  @param prox Gaussian kernel scale term. The kernel width is standardised by the median distance between Xsites, prox controls further scaleing of the kernel width.
#'
#'  @return data.frame of residuals for Xblock, rows of the returned data.frame align with rows of Xblock

R.block.multi <- function(Xblock,Xsite,res,prox = 1) {
  #At prox = 1, information in the distances is maximised
  #At larger prox values, sites are closer and more sites strongly influence the interpolated residuals.
  #At smaller prox values, sites are further apart and it is harder for sites to stongly influence the interpoloated residuals.
  #At extreme values of prox, all sites are equally "close" or "distant" and the interpolation just becomes the global mean of the residuals
  if(prox < 0){
    stop("R.block.multi: prox must be >=0")
  }

  if(any(res < 0)){
    stop("R.block.multi: all residuals (res) must be positive")
  }

  d <- fields::rdist(Xblock, Xsite) #distance matrix between Xblock and Xsite, which are in GF space

  #Gaussian kernel
  #Maximal information is when distances are at maximal gradient point, which
  #(by double derivative for exp(-(x)^2) ) is at sqrt(0.5)
  #So scale the data to put the median at sqrt(0.5)
  w <- prox * (median(d) / sqrt(0.5) )
  #then calculate the weighted distances from XBlock to Xsite
  d1 <- exp(-(d/w)^2)

  ##Now I need to take the weighted average sum(w_i * x_i)
  d2 <- (d1 %*% as.matrix(res))

  #now to divide by sum of distances
  ret <- d2 / matrix(data  = rowSums(d1), nrow = nrow(d2), ncol = ncol(d2), byrow = FALSE)

  return(ret)
}



#' Hotellings T^2 test for a pair of sites
#'
#' Calculates the Hotellings T^2 test statistic and resturns a list of relevant values.
#'
#' hotellings_summary() assumes that the covariance matricies are both diagonals and
#' represent standard deviation, not variance.
#'
#' @param x_bar mean of x sample
#' @param y_bar mean of y sample
#' @param x_var_vec standard deviation of x sample
#' @param y_var_vec standard deviation of y sample
#' @param x_n number of samples in x
#' @param y_n number of samples in y
#'
#' @return list of  \describe{
#'   \item{\code{df}}{degrees of freedom}
#'   \item{\code{t2}}{T^2 test statistic}
#'   \item{\code{f}}{F statistic}
#'   \item{\code{pval}}{p-value}
#'  }
hotellings_summary <- function(x_bar, y_bar, x_var_vec, y_var_vec, x_n = 1, y_n = 1){
  #The hotellings T^2 test:
  #T^2 = (x_bar - y_bar)' * (x_var/x_n + y_var/y_n)^-1 * (x_bar - y_bar)
  #Hotellings requires the inverse of the summed variance matrices
  #For speed of calculations, this function assumes:
  #residuals are standard deviations, which have already been squared before calling this function
  #Residuals are diagonals, there are no off diagonal terms in the variance matricies
  #Residuals represent populations, so there is no scaling for sample size when adding the variance matricies
  #so solve(x_cov/x_n + y_cov/y_n)


  mean_diff <- x_bar - y_bar
  k <- length(x_bar)
  df <- x_n+y_n-k-1
  inv <- 1/(x_var_vec^2 + y_var_vec^2) #The residuals come from distance to mean, but variance is distance to mean squared, so I must square the residuals
  #== solve(diag(x_cov_diag)/x_n + diag(y_cov_diag)/y_n) #yes, identical

  #t2 <-  (t(mean_diff)%*%solve(x_cov/x_n + y_cov/y_n)%*%mean_diff )
  #t2 <-  (t(mean_diff)%*%inv%*%mean_diff ) #thanks to diagonal assumption
  t2<- sum(mean_diff*inv*mean_diff)

  f  <- (df/((x_n + y_n - 2)*k))*t2

  pvalue = 1-pf(f, k, df)

  return(list(df = df, t2 = t2, f = f, pval = pvalue))
}

#' Hotellings T^2 statistic for many points
#'
#' Calculates the Hotellings T^2 test statistic for all pairs of points at once.
#' Vectorized for speed.
#'
#' Unlike hotellings_summary(), hotellings_bulk requires Procrustus residuals to be squared before
#' being passed in.
#'
#' @param means data.frame of sites by predictor, each row is the mean of a sample
#' @param res_sq data.frame of sites by diagonal variance, each row is the diagonal variance of a sample
#'
#' @return site by site p-value matrix of all pairs of sites
#'
#' @import data.table
hotellings_bulk <- function(means, res_sq){
  ##Find the pairwise distances between all pairs of points
  ##the sample sizes are assumed to approach infinity
  ##Assumes res have already been squared into variances
  if(any(dim(means) != dim(res_sq))){
    stop("means and res must have the same dimensions")
  }

  means <- data.table::as.data.table(means)
  res_sq <- data.table::as.data.table(res_sq)




  k <- ncol(means)
  ind <- data.table::CJ(x = seq.int(1,nrow(means)), y = seq.int(1, nrow(means)))
  ind <- ind[x < y,]

  mean_diffs <- means[ind$x,] - means[ind$y,] #may hit memory issues, this should create 3 matricies, each of size nrow(means)*(1 - nrow(means))/2 and throw away two.
  pair_vars <- res_sq[ind$x,] + res_sq[ind$y,] #as above, 3 matrices
  pair_inv <- 1/pair_vars

  t2_tmp <-  mean_diffs * pair_inv * mean_diffs
  t2 <- rowSums(t2_tmp)
  f <-  t2/k #(x_n + y_n -k - 1)/((x_n + y_n -2)*k) as x_n, y_n approach inf, becomes 1/k

  pvalue <- 1 - pf(f, k, Inf) #am I using the degrees of freedom properly?


  #fill in p matrix
  longform <- data.table(ind, p = pvalue)
  #add corners so dcast creates correct dimesions
  longform <- rbind(longform, data.table(x = c(1, nrow(means)), y = c(1, nrow(means)), p = c(0,0)))
  #I assume data.table can do this better than a for,for loop
  p_mat_raw <- data.table::dcast(longform, y ~ x, fill = 0, value.var = "p" )
  p_mat_raw[ , y := NULL]
  p_mat_raw <- as.matrix(p_mat_raw)


  p_mat <- as.dist(p_mat_raw, upper = TRUE, diag= TRUE)
  return(as.matrix(p_mat) + diag(nrow(means)))

}

