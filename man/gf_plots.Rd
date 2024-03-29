% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_plot.R
\name{gf_plots}
\alias{gf_plots}
\title{GradientForest summary plot collection}
\usage{
gf_plots(gf_fit, vars)
}
\arguments{
\item{gf_fit}{A gradientForest or combinedGradientForest object}

\item{vars}{integer or logical array specifying which variables to plot. 1 is the most important variable. If logical, must match length of names(importance(gf_fit))}
}
\description{
Generates a set of gradientForest plots
}
\details{
These are base R plots, apart from the variable importance plot
}
