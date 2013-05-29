\name{AddMasks}
\alias{AddMasks}
\title{Generate masks according to the given method and add them to the database}
\usage{
  AddMasks(db, calexp, method, channel = 2, scales = NULL,
    invert = F)
}
\arguments{
  \item{db}{the mask database object}

  \item{calexp}{the calcium experiment data object}

  \item{method}{what method should be used to generate
  masks to add. 'LoG' Laplacian of Gaussian, EqThresh
  thresholding of equalized image}

  \item{channel}{which data channel to use for computing
  the masks}

  \item{scales}{which smoothing scales to use for method
  (LoG)}

  \item{invert}{boolean, should the mean image be inverted
  to find dark regions instead of bright?}
}
\description{
  Generate masks according to the given method and add them
  to the database
}

