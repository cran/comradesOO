% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genericMethods.R
\docType{methods}
\name{clusterTableFolded}
\alias{clusterTableFolded}
\alias{clusterTableFolded,comradesDataSet-method}
\title{clusterTableFolded}
\usage{
clusterTableFolded(x)
}
\arguments{
\item{x}{A comradesDataSet object}
}
\value{
A table showing the vienna structures of each cluster
}
\description{
Extract the cluster coordinates with fold prediciton in data frame format
}
\examples{
cds = makeExampleComradesDataSet()

clusterTableFolded(cds)
}
