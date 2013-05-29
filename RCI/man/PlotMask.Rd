\name{PlotMask}
\alias{PlotMask}
\title{Plots a mask over an already plotted image}
\usage{
  PlotMask(mask, rgb = runif(3), alpha = 0.5, ...)
}
\arguments{
  \item{mask}{the specification of the mask}

  \item{rgb}{a vector of length 3 giving the color of the
  mask in RGB (defaults to random)}

  \item{alpha}{the alpha transparency value of the mask
  (between 0 and 1)}

  \item{...}{additional graphing parameters}
}
\value{
  NULL
}
\description{
  Plots a mask over an already plotted image
}
\details{
  Given a mask as either a matrix of logicals or a matrix
  with 1's on the mask, over-plot a semi-transparent
  colored region on an already plotted image.
}

