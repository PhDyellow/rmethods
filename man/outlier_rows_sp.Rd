% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_clean.R
\name{outlier_rows_sp}
\alias{outlier_rows_sp}
\title{Remove outliers from a species data.frame}
\usage{
outlier_rows_sp(dataset, range = 3, exclude_cols = NULL)
}
\arguments{
\item{dataset}{A data.frame of site by environmental obervations}

\item{range}{Keep points that lie within \code{range}*\eqn{\sigma} of the mean}

\item{exclude_cols}{Columns to ignore for identifying rows with outliers}
}
\description{
Returns a boolean vector of length \code{nrow(dataset)} showing which
observations to keep
}
\details{
Outliers are specified as being more than \code{range}*\eqn{\sigma} from the mean
within a column. Any row with at least one outlier is excluded.

Species observations are zero-inflated, so zero's are ignored for calculating
the mean and standard deviation. Also, only very large counts above the mean
are excluded as outliers, small counts below the mean are always kept.
}
