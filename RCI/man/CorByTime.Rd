\name{CorByTime}
\alias{CorByTime}
\title{Return array specifying the correlation matrix for a sliding
window of the data}
\usage{
  CorByTime(seriesmat, window = 500)
}
\arguments{
  \item{seriesmat}{a matrix with the calcium traces on the
  columns}

  \item{window}{the size of the sliding window}
}
\value{
  an array with the first two dimensions giving the
  correlation matrices and the third dimension indicating
  the start time of the window
}
\description{
  Return array specifying the correlation matrix for a
  sliding window of the data
}

