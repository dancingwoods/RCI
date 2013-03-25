\name{ClipImage}
\alias{ClipImage}
\title{Clips a border from around an image matrix}
\usage{
  ClipImage(image, border = NULL, size = NULL)
}
\arguments{
  \item{image}{the image matrix to clip}

  \item{border}{the size of the border to clip.  Must be
  less than half the image size}

  \item{size}{the resulting size of the image.}
}
\value{
  a matrix with the center (nrow-2*border) by
  (ncol-2*border) pixels of the image
}
\description{
  Clips a border from around an image matrix
}
\details{
  uses the size argument if given, else uses the border
  argument, else returns the original image
}

