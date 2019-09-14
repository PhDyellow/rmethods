#' Plot correlation heatmap from a data.frame
#'
#'
#' @param dataset A data.frame to plot
#'
#' @return a ggplot object. The ggplot object is not printed
#' @export
plot_cor_heatmap <- function(dataset){
  corr <- cor(dataset)
  pl <- ggcorrplot::ggcorrplot(corr, hc.order = TRUE,
                               type = "lower",
                               outline.color = "white" )
  return(pl)
}


#' Plot marginals of a dataset
#'
#' Returns a list of ggplot objects
#' showing the marginal histograms of each variable
#'
#' @param dataset A data.frame of variables to plot
#' @param binwidth Width of bins for the histogram plots. Defaults to 1.
#'
#' @return List of ggplot objects. The ggplot objects are not printed
#' @export
plot_marginals <- function(dataset, binwidth = 1){
  dNames <- names(dataset)

  var_plots <- lapply(X = dNames, FUN = function(x, dataset){
    var_plot <- ggplot2::ggplot(data = dataset[,x, drop=FALSE], mapping = ggplot2::aes_string(x = x)) +
      ggplot2::geom_histogram(binwidth = binwidth, boundary = 0)

    return(var_plot)

  }, dataset=dataset)
  return(var_plots)
}

#' Plots marginals of species
#'
#'
#' Species observations usually have very large numbers
#' of 0's, so plot_marginals_species excludes 0 from the
#' histogram for easier visualisation
#'
#' @param dataset A data.frame of variables to plot
#' @param binwidth Width of bins for the histogram plots. Defaults to 1.
#'
#' @return List of ggplot objects. The ggplot objects are not printed
#' @export
plot_marginals_species <- function(dataset, binwidth = 1){
  dNames <- names(dataset)

  sp_plots <- lapply(X = dNames, FUN = function(x, dataset){
    sp_plot <- ggplot2::ggplot(data = dataset[dataset[,x] !=0, x , drop=FALSE], mapping = ggplot2::aes_string(x = x)) +
      ggplot2::geom_histogram(binwidth = binwidth, boundary = 0) +
      ggplot2::xlab("Abundance") +
      ggplot2::ylab("Frequency") +
      ggplot2::ggtitle(sprintf("%s\n%g zero abundance observations excluded out of %g total observations", x, sum(dataset[,x] == 0), nrow(dataset)))
    return(sp_plot)
  }, dataset=dataset)

  return(sp_plots)
}


#' Plots pairs of variables
#'
#'
#' Generates a gglot scatterplot for each pair of variables in \code{dataset}
#'
#' @param dataset A data.frame of variables to plot
#'
#' @return List of ggplot objects. The ggplot objects are not printed
#' @export
plot_pairs <- function(dataset, coord_equal = FALSE){
  dNames <- names(dataset)
  pair_plots <- mapply(
    FUN = function(v1, v2, dataset){
      dNames <- names(dataset)
      if (match(v1, dNames) < match(v2, dNames)){
        pair_plot <- ggplot2::ggplot(data = dataset[, c(v1, v2)], mapping = ggplot2::aes_string(x = v1, y = v2)) +
          ggplot2::geom_point(shape = ".", alpha = 0.1)
        if(coord_equal){
          pair_plot <- pair_plot + ggplot2::coord_equal()
        }
        return(pair_plot)
      } else {
        return(NULL)
      }

    },
    rep(dNames, times = length(dNames)),
    rep(dNames, each = length(dNames)),
    MoreArgs = list(dataset=dataset),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  pair_plots <- pair_plots[-which(sapply(pair_plots, is.null))]
  return(pair_plots)
}


#' Plots gridded variables projected onto maps
#'
#'
#' Generates a gglot map for each variable in \code{dataset}.
#'
#' An sf_polygon can be included to provide context for the data.
#'
#' @param dataset A data.frame of variables to plot.
#' The data.frame is expected to have either "lat" and "lon" columns, or "x" and "y" columns.
#' @param lat_col Name of column containing latitude or y coordinates
#' @param lon_col Name of column containing longitude or x coordinates
#' @param sf_poly A polygon of class sf for plotting on the map.
#'
#' @return List of ggplot objects. The ggplot objects are not printed
#' @export
plot_maps_raster <- function(dataset, lat_col, lon_col, sf_poly = NULL){
  dNames <- names(dataset[, !(names(dataset) %in% c(lat_col, lon_col))])

  map_plots <- lapply(X = dNames, FUN = function(x, dataset){
    map_plot <- ggplot2::ggplot(data= dataset[,c(x, lon_col, lat_col)], ggplot2::aes_string(x = lon_col, y = lat_col, fill = x)) +
      ggplot2::geom_raster()  +
      ggplot2::labs(fill = x) +
      ggplot2::coord_sf()+
      ggthemes::theme_tufte()

    if(!is.null(sf_poly)){
      map_plot <- map_plot +  ggplot2::geom_sf(data = sf_poly, inherit.aes = FALSE, color = "black", fill= NA)
    }
    return(map_plot)
  }, dataset=dataset)
  return(map_plots)
}



#' Plots point data, projected onto maps
#'
#'
#' Generates a gglot map for each variable in \code{dataset}.
#'
#' An sf_polygon can be included to provide context for the data.
#'
#' @param dataset A data.frame of variables to plot.
#' The data.frame is expected to have either "lat" and "lon" columns, or "x" and "y" columns.
#' @param lat_col Name of column containing latitude or y coordinates
#' @param lon_col Name of column containing longitude or x coordinates
#' @param sf_poly A polygon of class sf for plotting on the map.
#'
#' @return List of ggplot objects. The ggplot objects are not printed
#' @export
plot_maps_points <- function(dataset, lat_col, lon_col, sf_poly = NULL){
  dNames <- names(dataset[, !(names(dataset) %in% c(lat_col, lon_col))])

  point_plots <- lapply(X = dNames, FUN = function(x, dataset){
    point_plot <- ggplot2::ggplot(data= dataset[dataset[x] > 0, c(x, lon_col, lat_col)],
                                  ggplot2::aes_string(x = lon_col, y = lat_col, size = x )) +
      ggplot2::geom_point(data= dataset[dataset[x] == 0, c(x, lon_col, lat_col)],
                          ggplot2::aes_string(x = lon_col, y = lat_col), size = .25, colour = "grey")  +
      ggplot2::geom_point(colour = "blue")  +
      ggplot2::scale_size_area(max_size = 5) +
      ggplot2::labs(fill = x) +
      ggplot2::coord_sf()+
      ggthemes::theme_tufte()

    if(!is.null(sf_poly)){
      point_plot <- point_plot + ggplot2::geom_sf(data = sf_poly,inherit.aes = FALSE, color = "black", fill= NA)
    }
    return(point_plot)
  }, dataset=dataset)
  return(point_plots)
}

#' Plot cluster pairs
#'
#'
#' Generates a list of plots, of each pair of variables,
#' showing clustered data with ellipses.
#'
#' @param dataset A data.frame of variables to plot.
#' @param cluster_model Cluster model (currently mclust) fitted to \code{dataset}
#' @param plots_vars Character vector of variables to plot. Defaults to all variables fitted in cluster_model
#' @param level The heigth of the contour line for the ellipses. Defaults to 68.3\% or 1 standard deviation from the mean
#' @param legend_thres When the number of clusters exceeds legend_thres, the legend is not plotted
#' @param alpha Alpha value of scatterplot points
#'
#' @return A list of ggplot object. The ggplot objects are not printed
#' @export
plot_cluster_pairs <- function(dataset, cluster_model, plot_vars = NULL,  level = 0.683, legend_thres = 10, alpha=0.3){
  if(is.null(plot_vars)){
    dNames <- colnames(cluster_model$data)
  } else {
    dNames <- plot_vars
  }

  cluster_pair_plots <- mapply(
    FUN = function(v1, v2, dataset, cluster_model, level, legend_thres, alpha){
      dNames <- colnames(cluster_model$data)
      if (match(v1, dNames) < match(v2, dNames)){
        pair_plot <- rmethods::plot_clusters(dataset = dataset, col_x = v1, col_y = v2,
                                             cluster_model = cluster_model,
                                             level = level, legend_thres = legend_thres, alpha = alpha)
        return(pair_plot)
      } else {
        return(NULL)
      }

    },
    rep(dNames, times = length(dNames)),
    rep(dNames, each = length(dNames)),
    MoreArgs = list(dataset=dataset, cluster_model = cluster_model,
                    level = level, legend_thres = legend_thres, alpha = alpha),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  cluster_pair_plots <- cluster_pair_plots[-which(sapply(cluster_pair_plots, is.null))]
  return(cluster_pair_plots)
}



#' Plot clustered data with cluster ellipses
#'
#'
#' Produces a scatter plot of the clustered data,
#' coloured by cluster. Clusters are shown
#' as ellipses at the 1 standard deviation contour line.
#' If not all columns in dataset were used for fitting the
#' cluster model, the relevant columns will be read from
#' the cluster_model
#'
#' Because each Gaussian is scaled by the unconditioned prior, pi, two curves are plotted.
#'
#' The solid filled curve is scaled by pi, and the dashed curve shows the ellipse before scaling by pi.
#'
#' @param dataset A dataframe of observations
#' @param col_x Column name for x axis
#' @param col_y Column name for y axis
#' @param cluster_model Cluster model (currently mclust) fitted to \code{dataset}
#' @param level The heigth of the contour line for the ellipses. Defaults to 68.3\%, or 1 standard deviation from the mean
#' @param legend_thres When the number of clusters exceeds legend_thres, the legend is not plotted
#' @param alpha Alpha value of scatterplot points
#'
#' @return A ggplot object. The ggplot object is not printed
#'
#' @export
#' @importFrom foreach foreach %do%
plot_clusters <- function(dataset, col_x, col_y, cluster_model, level = 0.683, legend_thres = 10, alpha=0.3){

  level_scaled <- level/cluster_model$G
  #generate ellipses
  ellipse_points_independent <- foreach(i = 1:cluster_model$G, .combine = rbind) %do% {
    ellipse_points <- ellipse::ellipse(x = cluster_model$parameters$variance$sigma[c(col_x, col_y),c(col_x, col_y),i],
                                       centre = cluster_model$parameters$mean[c(col_x, col_y),i],
                                       level = level_scaled)
    ellipse_points <- as.data.frame(ellipse_points)
    ellipse_points <- cbind(ellipse_points, data.frame(g=i))
  }

  #generate ellipses
  ellipse_points_scaled_pi <- foreach(i = 1:cluster_model$G, .combine = rbind) %do% {
    ellipse_points <- ellipse::ellipse(x = cluster_model$parameters$variance$sigma[c(col_x, col_y),c(col_x, col_y),i],
                                       centre = cluster_model$parameters$mean[c(col_x, col_y),i],
                                       level = (cluster_model$parameters$pro[i]*level_scaled*cluster_model$G))
    ellipse_points <- as.data.frame(ellipse_points)
    ellipse_points <- cbind(ellipse_points, data.frame(g=i))
  }

  classification <- predict(cluster_model, dataset[, colnames(cluster_model$data)])

  pl <- ggplot2::ggplot(data = dataset, mapping = ggplot2::aes_string(x=col_x, y=col_y,
                                                                      color=as.factor(classification$classification))) +
    ggplot2::scale_colour_manual(values = rainbow(cluster_model$G))  +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_path(data = ellipse_points_scaled_pi, mapping = ggplot2::aes_string(x=col_x, y=col_y, color = as.factor(ellipse_points_scaled_pi$g)), inherit.aes = FALSE) +
    ggplot2::geom_polygon(data = ellipse_points_scaled_pi, mapping = ggplot2::aes_string(x=col_x, y=col_y, fill = as.factor(ellipse_points_scaled_pi$g)), alpha = alpha,  inherit.aes = FALSE) +
    ggplot2::geom_path(data = ellipse_points_independent, mapping = ggplot2::aes_string(x=col_x, y=col_y, color = as.factor(ellipse_points_independent$g)), inherit.aes = FALSE, linetype = 2) +
    ggthemes::theme_tufte()+
    ggplot2::coord_fixed()

  if (cluster_model$G <= legend_thres){
    pl <- pl + ggplot2::guides(color=FALSE ) +  ggplot2::labs(fill = "Cluster")
  } else {
    pl <- pl + ggplot2::guides(color=FALSE, fill=FALSE )
  }

  return(pl)
}

#' Project and Plot clustered data with cluster ellipses
#'
#'
#' Produces a scatter plot of the clustered data,
#' coloured by cluster and projected onto the first two PCA axes.
#' Clusters are shown
#' as ellipses at the 1 standard deviation contour line.
#' If not all columns in dataset were used for fitting the
#' cluster model, the relevant columns will be read from
#' the cluster_model
#'
#' @param dataset A dataframe of observations
#' @param transform A string indicating transform to use from c("pca"). Defaults to "pca"
#' @param cluster_model Cluster model (currently mclust) fitted to \code{dataset}
#' @param level The heigth of the contour line for the ellipses. Defaults to 68.3\%, or 1 standard deviation from the mean
#' @param legend_thres When the number of clusters exceeds legend_thres, the legend is not plotted
#' @param alpha Alpha value of scatterplot points
#'
#' @return A ggplot object. The ggplot object is not printed
#'
#' @export
#' @importFrom foreach foreach %do%
plot_clusters_project <- function(dataset, transform = c("pca"),
                                  cluster_model, level = 0.683, legend_thres = 10, alpha=0.3){

  level_scaled <- level/cluster_model$G
  #Get projection
  trans <- prcomp(dataset[, colnames(cluster_model$data)])

  #Transform data and ellipses

  #https://math.stackexchange.com/questions/332441/affine-transformation-applied-to-a-multivariate-gaussian-random-variable-what
  #https://stats.stackexchange.com/questions/2592/how-to-project-a-new-vector-onto-pca-space
  trans_clusters <- list()

  combine_func <- function(...){abind::abind(..., along = 3)}
  trans_clusters$sigma <- foreach(i = 1:dim(cluster_model$parameters$variance$sigma)[3], .combine = combine_func) %do% {
    trans$rotation %*% cluster_model$parameters$variance$sigma[,,i] %*% t(trans$rotation)

  }
  trans_clusters$mean = predict(trans, t(cluster_model$parameters$mean))
  trans_clusters$pro = cluster_model$parameters$pro


  #generate ellipses
  ellipse_points_independent <- foreach(i = 1:cluster_model$G, .combine = rbind) %do% {
    ellipse_points <- ellipse::ellipse(x = trans_clusters$sigma[,,i],
                                       centre = trans_clusters$mean[i,],
                                       level = level_scaled)
    ellipse_points <- as.data.frame(ellipse_points)
    ellipse_points <- cbind(ellipse_points, data.frame(g=i))
  }

  #generate ellipses
  ellipse_points_scaled_pi <- foreach(i = 1:cluster_model$G, .combine = rbind) %do% {
    ellipse_points <- ellipse::ellipse(x = trans_clusters$sigma[,,i],
                                       centre = trans_clusters$mean[i,],
                                       level = (trans_clusters$pro[i]*level_scaled*cluster_model$G))
    ellipse_points <- as.data.frame(ellipse_points)
    ellipse_points <- cbind(ellipse_points, data.frame(g=i))
  }

  classification <- predict(cluster_model, dataset[, colnames(cluster_model$data)])
  col_x <- colnames(cluster_model$data)[1]
  col_y <- colnames(cluster_model$data)[2]
  pl <- ggplot2::ggplot(data = as.data.frame(trans$x[,1:2]), mapping = ggplot2::aes_string(x="PC1", y="PC2",
                                                                      color=as.factor(classification$classification))) +
    ggplot2::scale_colour_manual(values = rainbow(cluster_model$G))  +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_path(data = ellipse_points_scaled_pi, mapping = ggplot2::aes_string(x=col_x, y=col_y, color = as.factor(ellipse_points_scaled_pi$g)), inherit.aes = FALSE) +
    ggplot2::geom_polygon(data = ellipse_points_scaled_pi, mapping = ggplot2::aes_string(x=col_x, y=col_y, fill = as.factor(ellipse_points_scaled_pi$g)), alpha = alpha,  inherit.aes = FALSE) +
    ggplot2::geom_path(data = ellipse_points_independent, mapping = ggplot2::aes_string(x=col_x, y=col_y, color = as.factor(ellipse_points_independent$g)), inherit.aes = FALSE, linetype = 2) +
    ggthemes::theme_tufte()+
    ggplot2::coord_fixed()

  if (cluster_model$G <= legend_thres){
    pl <- pl + ggplot2::guides(color=FALSE ) +  ggplot2::labs(fill = "Cluster")
  } else {
    pl <- pl + ggplot2::guides(color=FALSE, fill=FALSE )
  }

  return(pl)
}



#' Plot map coloured by assigned cluster
#'
#'
#' Produces a map plot coloured by cluster.
#' A bounding polygon may be provided for context.
#' If not all columns in dataset were used for fitting the
#' cluster model, the relevant columns will be read from
#' the cluster_model
#'
#' @param dataset A dataframe of observations with an associated fitted model
#' @param col_x Column name for x axis
#' @param col_y Column name for y axis
#' @param depth_col Column name for bathymetry
#' @param depth_contour Vector of depths used for contour lines.
#' @param cluster_model Cluster model (currently mclust) fitted to \code{dataset}
#' @param legend_thres When the number of clusters exceeds legend_thres, the legend is not plotted
#' @param sf_poly A polygon of class sf for plotting on the map.
#'
#' @return A ggplot object. The ggplot object is not printed
#'
#' @export
plot_cluster_map <- function(dataset, col_x = "lon", col_y = "lat",
                             depth_col = NULL, depth_contour = NULL,
                             cluster_model, legend_thres = 10,
                             sf_poly = NULL){

  classification <- predict(cluster_model, dataset[, colnames(cluster_model$data)])

  pl <- ggplot2::ggplot(data = dataset,
               mapping = ggplot2::aes_string(x = col_x, y = col_y, fill = as.factor(classification$classification))) +
    ggplot2::scale_fill_manual(values = rainbow(cluster_model$G))  +
    ggplot2::geom_raster(hjust = 0, vjust = 1) +
    ggplot2::labs(fill = "cluster") +
    ggplot2::coord_sf()+
    ggthemes::theme_tufte()

  if(!is.null(sf_poly)){
    pl <- pl +  ggplot2::geom_sf(data = sf_poly, inherit.aes = FALSE, color = "black", fill= NA)
  }

  if(!is.null(depth_col) & !is.null(depth_contour)){
    pl <- pl + ggplot2::geom_contour(data = dataset[, c(col_x, col_y, depth_col)],
                          mapping = ggplot2::aes_string(x=col_x, y=col_y, z=depth_col), inherit.aes = FALSE, breaks = depth_contour)
  }


  if (cluster_model$G <= legend_thres){
    pl <- pl + ggplot2::labs(fill = "Cluster")
  } else {
    pl <- pl + ggplot2::guides(color=FALSE, fill=FALSE )
  }

  return(pl)

}

#' Plot clusters individually, coloured by confidence
#'
#'
#' Generates a list of ggplot objects, one per cluster.
#'
#' Each cluster is plotted alone, showing the sites each cluster is associated with, and how strongly the
#' site associates with each cluster.
#'
#' Across all clusters, each site membership sums to 1
#'
#' @param dataset A dataframe of observations with an associated fitted model
#' @param col_x Column name for x axis
#' @param col_y Column name for y axis
#' @param depth_col Column name for bathymetry
#' @param depth_contour Vector of depths used for contour lines.
#' @param cluster_model Cluster model (currently mclust) fitted to \code{dataset}
#' @param sf_poly A polygon of class sf for plotting on the map.
#'
#' @return A ggplot object. The ggplot object is not printed
#'
#' @export
#' @importFrom foreach foreach %do%
plot_cluster_own <- function(dataset, col_x = "lon", col_y = "lat",
                             depth_col = NULL, depth_contour = NULL,
                             cluster_model, sf_poly = NULL){
  classification <- predict(cluster_model, dataset[, colnames(cluster_model$data)])
  map_plots <- foreach(i = 1:cluster_model$G) %do% {
    pl <- ggplot2::ggplot(data= data.frame(dataset[,c(col_x, col_y)], z =  classification$z[,i]),
                          ggplot2::aes_string(x = col_x, y = col_y, fill = "z")) +
      ggplot2::geom_raster()  +
      ggplot2::scale_fill_viridis_c(limits=c(0,1))  +
      ggplot2::labs(fill = "z") +
      ggplot2::coord_sf()+
      ggthemes::theme_tufte()

    if(!is.null(sf_poly)){
      pl <- pl +  ggplot2::geom_sf(data = sf_poly, inherit.aes = FALSE, color = "black", fill= NA)
    }

    if(!is.null(depth_col) & !is.null(depth_contour)){
      pl <- pl + ggplot2::geom_contour(data = dataset[, c(col_x, col_y, depth_col)],
                                       mapping = ggplot2::aes_string(x=col_x, y=col_y, z=depth_col), inherit.aes = FALSE, breaks = depth_contour)
    }

    return(pl)
  }

  return(map_plots)

}

#' Plot cluster map, but hide sites with low confidence
#'
#'
#' @export
plot_cluster_map_ci <- function(x, rawx, cluster_model, legend_thres = 10, landPoly, confThres){
  stop("not implemented yet")
  keepSites <- cluster_model$uncertainty < (1-confThres)


  pl <- ggplot(data = data.frame(lon= x[["lon"]][keepSites], lat= x[["lat"]][keepSites], c = cluster_model$classification[keepSites]), mapping = aes(x = lon, y = lat, fill = as.factor(c))) +
    scale_fill_manual(values = rainbow(cluster_model$G))  +
    geom_raster(hjust = 0, vjust = 1) +
    labs(fill = "cluster") +
    geom_sf(data = landPoly, inherit.aes = FALSE, color = "black", fill= NA) +
    coord_sf()+
    theme_tufte() +
    geom_contour(data = rawx[keepSites, .(x,y,MS_bathy_5m)], mapping = aes(x=x, y=y, z=MS_bathy_5m), inherit.aes = FALSE, breaks = c(-200))

  if (cluster_model$G <= legend_thres){
    pl <- pl + labs(fill = "Cluster")
  } else {
    pl <- pl + guides(color=FALSE, fill=FALSE )
  }

  return(pl)

}




#' Gradient Forest Variable importance weighted by R^2
#'
#'
#' Plots variable importance weighted by R^2.
#'
#' Uses S3 dispatch to handle both gradientForest and combinedGradientForest
#'
#' @param x A gradient forest model
#'
#' @return ggplot2 plot object showing R^2 weighted importance
#' @export
plot_gf_importance <- function(x) UseMethod("plot_gf_importance")

plot_gf_importance.default <- function(x, ...){

  warning(paste("PLOT_GF_IMPORTANCE does not know how to handle object of class ",
                class(x),
                "and can only be used on classes gradientForest and combinedGradientForest"))

}



plot_gf_importance.gradientForest <- function(x){
  imp_vars <- gradientForest::importance.gradientForest(x, type = "Weighted")
  imp_vars_df <- data.frame(env_vars = names(imp_vars), imp = imp_vars)
  imp_vars_df <- imp_vars_df[order(imp_vars_df$imp) ,]

  pl <- ggplot2::ggplot(data = imp_vars_df, mapping = ggplot2::aes(x = env_vars, y = imp)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(limits = imp_vars_df$env_vars) +
    ggthemes::theme_tufte() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "", y = "Variable Importance") +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank())

  return(pl)
}

plot_gf_importance.combinedGradientForest <- function(x){
  imp_vars <- gradientForest::importance.combinedGradientForest(x, type = "Weighted")

  imp_vars_df <- data.frame(env_vars = names(imp_vars), imp = imp_vars)
  imp_vars_df <- imp_vars_df[order(imp_vars_df$imp) ,]

  pl <- ggplot2::ggplot(data = imp_vars_df, mapping = ggplot2::aes(x = env_vars, y = imp)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(limits = imp_vars_df$env_vars) +
    ggthemes::theme_tufte() +
    ggplot2::coord_flip()+
    ggplot2::labs(x = "", y = "Variable Importance") +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank())

return(pl)
}



#' GradientForest summary plot collection
#'
#'
#' Generates a set of gradientForest plots
#'
#' These are base R plots, apart from the variable importance plot
#' @param gf_fit A gradientForest or combinedGradientForest object
#' @param vars integer or logical array specifying which variables to plot. 1 is the most important variable. If logical, must match length of names(importance(gf_fit))
#' @export
gf_plots <- function(gf_fit, vars){

  plot_gf_importance(gf_fit)

  if(class(gf_fit)[1] == "gradientForest"){

    plot(gf_fit, plot.type = "S", imp.vars = names(importance(gf_fit))[vars],
         leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
         cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5,
                                                                 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))

    plot(gf_fit, plot.type = "C", imp.vars = names(importance(gf_fit))[vars],
         show.overall = T, show.species = F, legend = T, leg.posn = "topleft",
         common.scale = T,
         leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
         cex.axis = 0.6, line.ylab = 0.9,
         par.args =
           list(mgp = c(1.5,0.5, 0),mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))

  } else if (class(gf_fit)[1] == "combinedGradientForest") {
    plot(gf_fit,plot.type="Predictor.Ranges",imp.vars = names(importance(gf_fit))[vars])
    plot(gf_fit,plot.type="Predictor.Density",imp.vars = names(importance(gf_fit))[vars])
    plot(gf_fit,plot.type="Cumulative.Importance",imp.vars = names(importance(gf_fit))[vars])
    plot(gf_fit,plot.type="Performance")



  }
}
