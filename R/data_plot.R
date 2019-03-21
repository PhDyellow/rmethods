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
