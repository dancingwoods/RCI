\name{ClusterCells}
\alias{ClusterCells}
\title{Cluster segmented ROIs based on correlation or phase distance using k-means}
\usage{
  ClusterCells(calexp, mask, k, criteria = "cor",
    freq = c(0.78, 0.81), dt = 0.1247232)
}
\arguments{
  \item{calexp}{the calexp object with the data}

  \item{mask}{a mask identifying the cells to be clustered.
  Each unique non-zero/NA value in the mask indicates a
  cell to be clustered.}

  \item{k}{the number of clusters to find}

  \item{criteria}{the criteria to use for clustering --
  'cor' (correlation) 'phase' (phase of frequency specified
  in freq)}

  \item{freq}{the frequency band to use to extract the
  phase for phase-clustering}
}
\description{
  Cluster segmented ROIs based on correlation or phase
  distance using k-means
}

