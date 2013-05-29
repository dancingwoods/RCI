\name{GetSeries}
\alias{GetSeries}
\title{Return average time series for each cell in a mask}
\usage{
  GetSeries(mask, calexp, channel = 2)
}
\arguments{
  \item{calexp}{the calexp with the data}

  \item{mask}{the mask identifying cells.  Each unique
  non-zero/NA value in the mask indicates a cell to be
  clustered.}

  \item{channel}{the channel to get the cell traces from}
}
\description{
  Return average time series for each cell in a mask
}

