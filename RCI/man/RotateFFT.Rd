\name{RotateFFT}
\alias{RotateFFT}
\title{Rotates an image by the given angle using a sequence of Fourier domain shears as
described in Eddy 1996.}
\usage{
  RotateFFT(img, theta, fdomain = FALSE)
}
\arguments{
  \item{img}{the image to rotate}

  \item{theta}{the angle to rotate the image}

  \item{fdomain}{is the image given already in the Fourier
  domain?  It will be returned in the same domain as given
  (passing in the Fourier domain is helpful to reduce
  superfluous transforms if performing additional
  operations in the Fourier domain).}
}
\description{
  Rotates an image by the given angle using a sequence of
  Fourier domain shears as described in Eddy 1996.
}

