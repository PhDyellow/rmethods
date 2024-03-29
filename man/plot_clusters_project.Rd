% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_plot.R
\name{plot_clusters_project}
\alias{plot_clusters_project}
\title{Project and Plot clustered data with cluster ellipses}
\usage{
plot_clusters_project(dataset, transform = c("pca"), cluster_model,
  level = 0.683, legend_thres = 10, alpha = 0.3)
}
\arguments{
\item{dataset}{A dataframe of observations}

\item{transform}{A string indicating transform to use from c("pca"). Defaults to "pca"}

\item{cluster_model}{Cluster model (currently mclust) fitted to \code{dataset}}

\item{level}{The heigth of the contour line for the ellipses. Defaults to 68.3\%, or 1 standard deviation from the mean}

\item{legend_thres}{When the number of clusters exceeds legend_thres, the legend is not plotted}

\item{alpha}{Alpha value of scatterplot points}
}
\value{
A ggplot object. The ggplot object is not printed
}
\description{
Produces a scatter plot of the clustered data,
coloured by cluster and projected onto the first two PCA axes.
Clusters are shown
as ellipses at the 1 standard deviation contour line.
If not all columns in dataset were used for fitting the
cluster model, the relevant columns will be read from
the cluster_model
}
