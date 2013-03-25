\name{OptimShift}
\alias{OptimShift}
\title{INTERNAL
Computes the rigid body motion alignment parameters by optimizing some error function
comparing the two images. (uses optimization routines in the neldermead package)}
\usage{
  OptimShift(img1, img2, taper = TRUE, error = "mse",
    startval = c(0.1, 0.1, 0), pocstart = TRUE,
    bigsize = NULL)
}
\arguments{
  \item{img1}{the reference image}

  \item{img2}{the image to align}

  \item{taper}{boolean, should the images be tapered before
  aligning}

  \item{bigsize}{the size of the array in which to embed
  the tapered images (defaults to next power of 2)}

  \item{error}{the error function to use. Can be "mse" for
  mean squared error, mae" for mean absolute error, or
  "cor" for correlation.}

  \item{startval}{a length 3 vector giving the inital
  values for the optimization (xshift, yshift, theta)}

  \item{pocstart}{should the POC method be used to
  initialize the start values}
}
\value{
  a vector of length 3 giving the translation and rotation
  estimates
}
\description{
  INTERNAL Computes the rigid body motion alignment
  parameters by optimizing some error function comparing
  the two images. (uses optimization routines in the
  neldermead package)
}

