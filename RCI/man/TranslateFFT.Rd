\name{TranslateFFT}
\alias{TranslateFFT}
\title{INTERNAL
Shifts an image by the given (fractional pixel) amounts}
\usage{
  TranslateFFT(img, xshift, yshift)
}
\arguments{
  \item{img}{the image (matrix) to shift}

  \item{xshift}{the amount to shift the in x dimension
  (columns)}

  \item{yshift}{the amount to shift in the y dimension
  (rows)}
}
\value{
  the shifted image (matrix)
}
\description{
  INTERNAL Shifts an image by the given (fractional pixel)
  amounts
}
\details{
  Uses the shift theorem to shift the given image by
  transforming to the Fourier domain.  The shift can be
  sub-pixel, resulting in Fourier interpolation.
}

