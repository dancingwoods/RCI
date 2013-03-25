\name{TranslateFFT}
\alias{TranslateFFT}
\title{Shifts an image by the given (fractional pixel) amounts}
\usage{
  TranslateFFT(img, xshift, yshift, fdomain = FALSE)
}
\arguments{
  \item{img}{the image (matrix) to shift}

  \item{xshift}{the amount to shift the in x dimension
  (columns)}

  \item{yshift}{the amount to shift in the y dimension
  (rows)}

  \item{fdomain}{is the image given in the Fourier domain?
  It will be returned in the same domain as given (passing
  in the Fourier domain is helpful to reduce superfluous
  transforms if performing additional operations in the
  Fourier domain).}
}
\value{
  the shifted image (matrix)
}
\description{
  Shifts an image by the given (fractional pixel) amounts
}
\details{
  Uses the shift theorem to shift the given image by
  transforming to the Fourier domain.  The shift can be
  sub-pixel, resulting in Fourier interpolation.
}

