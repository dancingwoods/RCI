\name{FilterVector}
\alias{FilterVector}
\title{INTERNAL
Filters a vector by frequency using a butterworth filter}
\usage{
  FilterVector(vec, low, high, order = 8, dt = 0.1247232,
    type = "BP")
}
\arguments{
  \item{vec}{the vector to filter}

  \item{low}{the lower value of the filter}

  \item{high}{the higher value of the filter}

  \item{order}{the order of the butterworth filter}

  \item{dt}{the time (in seconds) of one datapoint.
  1/frequency in hz}

  \item{type}{the type of filter, defaults to "BP" bandpass
  filter. Can also choose other filters offered by the
  butfilt function}
}
\value{
  the filtered vector
}
\description{
  INTERNAL Filters a vector by frequency using a
  butterworth filter
}

