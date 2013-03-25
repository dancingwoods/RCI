\name{EmbedImage}
\alias{EmbedImage}
\title{Embeds an image matrix in a larger matrix with a border of 0's}
\usage{
  EmbedImage(image, border = NULL, size = NULL)
}
\arguments{
  \item{image}{the image to embed}

  \item{border}{the width of the border to add around the
  edges}

  \item{size}{the resulting size of the image - this must
  be bigger than the dimensions of image}
}
\value{
  a matrix of size (nrow+2*border) by (ncol+ 2*border)
}
\description{
  Embeds an image matrix in a larger matrix with a border
  of 0's
}
\details{
  uses the size argument if given, else uses the border
  argument, else returns the original image
}

