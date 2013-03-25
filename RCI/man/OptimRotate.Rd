\name{OptimRotate}
\alias{OptimRotate}
\title{INTERNAL
Uses optimization of an objective function to compute the best alignment rotation
between two images}
\usage{
  OptimRotate(img1, img2, taper = TRUE, error = "mse",
    searchrange = c(-0.1, 0.1))
}
\arguments{
  \item{img1}{the reference image}

  \item{img2}{the image to align}

  \item{taper}{should the images be tapered before the
  rotation is computed (hanning window)}

  \item{error}{objective function to be used - "mse" mean
  squared error, "mae" mean absolute error, "cor"
  correlation}

  \item{searchrange}{the range of rotations to search over
  in the optimization}
}
\value{
  a real valued estimate of the optimal alignment rotation
}
\description{
  INTERNAL Uses optimization of an objective function to
  compute the best alignment rotation between two images
}

