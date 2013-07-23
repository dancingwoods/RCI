\name{PhaseDistMat}
\alias{PhaseDistMat}
\title{Get phase distance matrix}
\usage{
  PhaseDistMat(mat, low = 0.789, high = 0.791,
    dt = 0.1247232)
}
\arguments{
  \item{mat}{a matrix with the calcium traces on the
  columns}

  \item{low}{the lower bound of the frequency to consider}

  \item{high}{the upper bound of the frequency to consider}

  \item{dt}{the sampling rate of the calcium traces}
}
\value{
  a matrix of phase distances between the calcium traces
}
\description{
  Get phase distance matrix
}

