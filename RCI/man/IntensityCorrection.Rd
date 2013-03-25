\name{IntensityCorrection}
\alias{IntensityCorrection}
\title{Performs intensity correction on the given calcium experiment}
\usage{
  IntensityCorrection(calexp, cortype = "ar", order = 25,
    naclip = T)
}
\arguments{
  \item{calexp}{the data to be corrected is in the $data
  element of this calexp object}

  \item{cortype}{the type of correction to perform.  `ar`
  for autoregressive filter}

  \item{order}{the order of the model to fit (for ar type)}

  \item{naclip}{should NAs produced at the beginning of the
  experiment be clipped off (by AR model, for instance)}
}
\description{
  Performs intensity correction on the given calcium
  experiment
}

