#' Plot correlation heatmap from a data.frame
#'
#' @param dataset A data.frame to plot
#'
#' @return a ggplot object. The ggplot object is not printed
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
#' Species observations usually have very large numbers
#' of 0's, so plot_marginals_species excludes 0 from the
#' histogram for easier visualisation
#'
#' @param dataset A data.frame of variables to plot
#' @param binwidth Width of bins for the histogram plots. Defaults to 1.
#'
#' @return List of ggplot objects. The ggplot objects are not printed
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
#' Generates a gglot scatterplot for each pair of variables in \code{dataset}
#'
#' @param dataset A data.frame of variables to plot
#'
#' @return List of ggplot objects. The ggplot objects are not printed
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
plot_maps_raster <- function(dataset, lat_col, lon_col, sf_poly = NULL){
  dNames <- names(dataset[, !(names(dataset) %in% c(lat_col, lon_col))])

  map_plots <- lapply(X = dNames, FUN = function(x, dataset){
    map_plot <- ggplot2::ggplot(data= dataset[,c(x, lon_col, lat_col)], ggplot2::aes_string(x = lon_col, y = lat_col, fill = x)) +
      ggplot2::geom_raster()  +
      ggplot2::geom_sf(data = sf_poly, inherit.aes = FALSE, color = "black", fill= NA) +
      ggplot2::labs(fill = x) +
      ggplot2::coord_sf()+
      ggthemes::theme_tufte()
    return(map_plot)
  }, dataset=dataset)
  return(map_plots)
}



#' Plots point data, projected onto maps
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
plot_maps_points <- function(dataset, lat_col, lon_col, sf_poly = NULL){
  dNames <- names(dataset[, !(names(dataset) %in% c(lat_col, lon_col))])

  point_plots <- lapply(X = dNames, FUN = function(x, dataset){
    point_plot <- ggplot2::ggplot(data= dataset[dataset[x] > 0, c(x, lon_col, lat_col)],
                                  ggplot2::aes_string(x = lon_col, y = lat_col, size = x )) +
      ggplot2::geom_point(data= dataset[dataset[x] == 0, c(x, lon_col, lat_col)],
                          ggplot2::aes_string(x = lon_col, y = lat_col), size = .25, colour = "grey")  +
      ggplot2::geom_point(colour = "blue")  +
      ggplot2::scale_size_area(max_size = 5) +
      ggplot2::geom_sf(data = sf_poly,inherit.aes = FALSE, color = "black", fill= NA) +
      ggplot2::labs(fill = x) +
      ggplot2::coord_sf()+
      ggthemes::theme_tufte()
    return(point_plot)
  }, dataset=dataset)
  return(point_plots)
}

#' Plot cluster pairs
#'
#' Generates a list of plots, of each pair of variables,
#' showing clustered data with ellipses.
#'
#' @param dataset A data.frame of variables to plot.
#' @param cluster_model Cluster model (currently mclust) fitted to \code{dataset}
#' @param level The heigth of the contour line for the ellipses. Defaults to 68.3%, or 1 standard deviation from the mean
#' @param legend_thres When the number of clusters exceeds legend_thres, the legend is not plotted
#' @param alpha Alpha value of scatterplot points
#'
#' @return A list of ggplot object. The ggplot objects are not printed
plot_cluster_pairs <- function(dataset, cluster_model,  level = 0.683, legend_thres = 10, alpha=0.3){
  dNames <- colnames(cluster_model$data)
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
#' Produces a scatter plot of the clustered data,
#' coloured by cluster. Clusters are shown
#' as ellipses at the 1 standard deviation contour line.
#' If not all columns in dataset were used for fitting the
#' cluster model, the relevant columns will be read from
#' the cluster_model
#'
#' @param dataset A dataframe of observations
#' @param col_x Column name for x axis
#' @param col_y Column name for y axis
#' @param cluster_model Cluster model (currently mclust) fitted to \code{dataset}
#' @param level The heigth of the contour line for the ellipses. Defaults to 68.3%, or 1 standard deviation from the mean
#' @param legend_thres When the number of clusters exceeds legend_thres, the legend is not plotted
#' @param alpha Alpha value of scatterplot points
#'
#' @return A ggplot object. The ggplot object is not printed
#'
#' @importFrom foreach foreach %do%
plot_clusters <- function(dataset, col_x, col_y, cluster_model, level = 0.683, legend_thres = 10, alpha=0.3){

  #generate ellipses
  ellipse_points <- foreach(i = 1:cluster_model$G, .combine = rbind) %do% {
    ellipse_points <- ellipse::ellipse(x = cluster_model$parameters$variance$sigma[c(col_x, col_y),c(col_x, col_y),i], centre = cluster_model$parameters$mean[c(col_x, col_y),i], level = level^2)
    ellipse_points <- as.data.frame(ellipse_points)
    ellipse_points <- cbind(ellipse_points, data.frame(g=i))
  }

  classification <- predict(cluster_model, dataset[, colnames(cluster_model$data)])

  pl <- ggplot2::ggplot(data = dataset, mapping = ggplot2::aes_string(x=col_x, y=col_y,
                                                                      color=as.factor(classification$classification))) +
    ggplot2::scale_colour_manual(values = rainbow(cluster_model$G))  +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_path(data = ellipse_points, mapping = aes_string(x=col_x, y=col_y, color = as.factor(ellipse_points$g)), inherit.aes = FALSE) +
    ggplot2::geom_polygon(data = ellipse_points, mapping = aes_string(x=col_x, y=col_y, fill = as.factor(ellipse_points$g)), alpha = alpha,  inherit.aes = FALSE) +
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
    ggplot2::geom_sf(data = sf_poly, inherit.aes = FALSE, color = "black", fill= NA) +
    ggplot2::coord_sf()+
    ggthemes::theme_tufte()

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

#' Plot cluster map, but hide sites with low confidence
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
