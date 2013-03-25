\name{OptimTranslate}
\alias{OptimTranslate}
\title{Uses optimization of an objective function to compute the best alignment translation
between two images (uses optimization routines in the neldermead package)}
\usage{
  OptimTranslate(img1, img2, taper = TRUE, error = "mse",
    startval = c(0.1, 0.1), bigsize = NULL)
}
\arguments{
  \item{img1}{the reference image}

  \item{img2}{the image to align}

  \item{taper}{should the images be tapered before the
  rotation is computed (hanning window)}

  \item{error}{objective function to be used - "mse" mean
  squared error, "mae" mean absolute error, "cor"
  correlation}

  \item{startval}{the inital estimate of the shift
  parameters}
}
\value{
  a real valued vector of length 2, giving estimates of x
  and y translation
}
\description{
  Uses optimization of an objective function to compute the
  best alignment translation between two images (uses
  optimization routines in the neldermead package)
}

