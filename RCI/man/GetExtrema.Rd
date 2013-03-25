\name{GetExtrema}
\alias{GetExtrema}
\title{Finds the extrema in an image.}
\usage{
  GetExtrema(image, maxima = T)
}
\arguments{
  \item{image}{the image matrix}

  \item{maxima}{boolean, should this function find maxima
  (default).  If false, finds minima}
}
\value{
  a matrix with 1 at maxima (or minima) and 0 elsewhere
}
\description{
  Finds the extrema in an image.
}

