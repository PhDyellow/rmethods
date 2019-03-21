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


plot_clusters <- function(x, colX, colY, clustMod, level = 0.683, legendThres = 10, alpha=0.3){

  #generate ellipses
  ellipsePoints <- foreach(i = 1:clustMod$G, .combine = rbind) %do% {
    ellipsePoints <- ellipse(x = clustMod$parameters$variance$sigma[,,i], centre = clustMod$parameters$mean[,i], level = level^2)
    ellipsePoints <- as.data.frame(ellipsePoints)
    ellipsePoints <- cbind(ellipsePoints, data.frame(g=i))
  }

  pl <- ggplot(data = data.frame(x = x[[colX]], y = x[[colY]], c = clustMod$classification), mapping = aes(x=x, y=y, color=as.factor(c))) +
    scale_colour_manual(values = rainbow(clustMod$G))  +
    geom_point(shape = ".", alpha = alpha) +
    geom_path(data = ellipsePoints, mapping = aes_string(x=colX, y=colY, color = as.factor(ellipsePoints$g)), inherit.aes = FALSE) +
    geom_polygon(data = ellipsePoints, mapping = aes_string(x=colX, y=colY, fill = as.factor(ellipsePoints$g)), alpha = alpha,  inherit.aes = FALSE) +
    theme_tufte()+
    coord_fixed()

  if (clustMod$G <= legendThres){
    pl <- pl + guides(color=FALSE ) +  labs(fill = "Cluster")
  } else {
    pl <- pl + guides(color=FALSE, fill=FALSE )
  }

  return(pl)


}

plot_cluster_map <- function(x, rawx, clustMod, legendThres = 10, landPoly){
  pl <- ggplot(data = data.frame(lon= x[["lon"]], lat= x[["lat"]], c = clustMod$classification), mapping = aes(x = lon, y = lat, fill = as.factor(c))) +
    scale_fill_manual(values = rainbow(clustMod$G))  +
    geom_raster(hjust = 0, vjust = 1) +
    labs(fill = "cluster") +
    geom_sf(data = landPoly, inherit.aes = FALSE, color = "black", fill= NA) +
    coord_sf()+
    theme_tufte() +
    geom_contour(data = rawx[, .(x,y,MS_bathy_5m)], mapping = aes(x=x, y=y, z=MS_bathy_5m), inherit.aes = FALSE, breaks = c(-200))


  if (clustMod$G <= legendThres){
    pl <- pl + labs(fill = "Cluster")
  } else {
    pl <- pl + guides(color=FALSE, fill=FALSE )
  }

  return(pl)

}

plot_cluster_map_ci <- function(x, rawx, clustMod, legendThres = 10, landPoly, confThres){
  keepSites <- clustMod$uncertainty < (1-confThres)


  pl <- ggplot(data = data.frame(lon= x[["lon"]][keepSites], lat= x[["lat"]][keepSites], c = clustMod$classification[keepSites]), mapping = aes(x = lon, y = lat, fill = as.factor(c))) +
    scale_fill_manual(values = rainbow(clustMod$G))  +
    geom_raster(hjust = 0, vjust = 1) +
    labs(fill = "cluster") +
    geom_sf(data = landPoly, inherit.aes = FALSE, color = "black", fill= NA) +
    coord_sf()+
    theme_tufte() +
    geom_contour(data = rawx[keepSites, .(x,y,MS_bathy_5m)], mapping = aes(x=x, y=y, z=MS_bathy_5m), inherit.aes = FALSE, breaks = c(-200))

  if (clustMod$G <= legendThres){
    pl <- pl + labs(fill = "Cluster")
  } else {
    pl <- pl + guides(color=FALSE, fill=FALSE )
  }

  return(pl)

}
