\name{GetAllSeries}
\alias{GetAllSeries}
\title{Get average time series of the given class in the segmentation stored in a mask database}
\usage{
  GetAllSeries(db, calexp, classids, chan = 2)
}
\arguments{
  \item{db}{the database}

  \item{calexp}{the calcium experiment with the data to use
  to extract the series}

  \item{classids}{a vector of ids specifying which types of
  ROI to extract traces for}

  \item{chan}{the channel to use for activity traces
  (defaults to 2)}
}
\description{
  Get average time series of the given class in the
  segmentation stored in a mask database
}

