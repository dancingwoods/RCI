\name{PlotClustering}
\alias{PlotClustering}
\title{Plots a given clustering}
\usage{
  PlotClustering(db, ids, clust, chan = 2,
    cols = list(c(0, 0, 1), c(0, 1, 0)))
}
\arguments{
  \item{db}{the database with the segmentation}

  \item{ids}{a vector of ids of the masks involved in the
  clustering}

  \item{clust}{a vector with integers indicating the
  clusters}

  \item{chan}{the channel to use as the background image}

  \item{cols}{a list of vectors of length 3 giving the rgb
  values for each cluster}
}
\value{
  NULL
}
\description{
  Plots a given clustering
}

