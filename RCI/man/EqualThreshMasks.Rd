\name{EqualThreshMasks}
\alias{EqualThreshMasks}
\title{Generates a set of masks using thresholding of the sliding histogram equalized
version of an image}
\usage{
  EqualThreshMasks(image, thresh, radius = 8,
    fullmax = 4096, sparse = T)
}
\arguments{
  \item{image}{the image to use to generate masks}

  \item{thresh}{the threshold to use (pixels above thresh
  in equalized image are found)}

  \item{radius}{the radius for the window used for
  equalization}

  \item{fullmax}{the maximum value possible (for
  equalization)}

  \item{sparse}{boolean, should the function return sparse
  masks instead of a matrix for plotting}
}
\value{
  a matrix with unique integers at mask locations and 0 in
  the background, or a a sparse masks object
}
\description{
  Generates a set of masks using thresholding of the
  sliding histogram equalized version of an image
}

