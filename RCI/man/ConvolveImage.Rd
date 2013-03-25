\name{ConvolveImage}
\alias{ConvolveImage}
\title{Convolves an image with the given kernel matrix}
\usage{
  ConvolveImage(image, kernel, circular = T)
}
\arguments{
  \item{image}{a matrix with the image}

  \item{kernel}{a matrix with the kernel (should be smaller
  than the image)}

  \item{circular}{boolean, should the convolution be
  circular (default) or should the image be padded with
  zeros to prevent circular convolution}
}
\value{
  a matrix of the same size as image with the convolved
  image
}
\description{
  Convolves an image with the given kernel matrix
}
\details{
  Uses Fourier methods to convolve the given image with the
  given kernel
}

