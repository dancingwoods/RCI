\name{GetDataFeatures}
\alias{GetDataFeatures}
\title{INTERNAL
Computes the features related to the data under a mask, adding them to the database}
\usage{
  GetDataFeatures(db, data, cormat)
}
\arguments{
  \item{db}{a database connection}

  \item{data}{the data array for this experiment}

  \item{cormat}{the pixel-pixel corrlations for channel 2
  for this data}
}
\value{
  NULL
}
\description{
  INTERNAL Computes the features related to the data under
  a mask, adding them to the database
}
\details{
  Computes features of all masks in the database or a list
  of masks specified by id.  The features computed are
  currently:\cr var1 - the variance of the pixel means for
  channel 1 \cr var2 - the variance of the pixel means for
  channel 2\cr var1eq - the variance of the pixel means for
  the equalized version of channel 1\cr var2eq - the
  variance of the pixel means for the equalized version of
  channel 2\cr mean1eq - the mean of the pixel means for
  the equalized version of channel 1\cr mean2eq - the mean
  of the pixel means for the equalized version of channel
  2\cr cor2 - the mean pixel-pixel correlation between the
  map pixels in channel 2\cr cor2min - the min pixel-pixel
  correlation between the map pixels in channel 2\cr
  cor2max - the max pixel-pixel correlation between the map
  pixels in channel 2\cr
}

