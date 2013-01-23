\name{FFTXCor}
\alias{FFTXCor}
\title{INTERNAL
Computes sub-pixel shifts values using FFT}
\usage{
  FFTXCor(img1, img2, upsamp = 1, taper = 0)
}
\arguments{
  \item{img1}{matrix giving the first image (the
  reference)}

  \item{img2}{matrix giving the second image (to be
  shifted)}

  \item{upsamp}{the factor by which the fft matrix should
  be expanded}

  \item{taper}{number of pixels to taper the data on the
  edges of the image}
}
\value{
  a vector of length 2 giving the magnitude of the estimted
  x and y shift returns NA in the case of improper input
}
\description{
  INTERNAL Computes sub-pixel shifts values using FFT
}
\details{
  Computes the sub-pixel shifts by computing the upsampled
  cross-correlation between the two images and finding the
  maximum. Computes the upsampled cross-correlation by
  embedding the product of FT(img1)* and FFT(img2) in a
  larger matrix of 0's determined by the upsampling factor.
}

