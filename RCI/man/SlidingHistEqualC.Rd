\name{SlidingHistEqualC}
\alias{SlidingHistEqualC}
\title{INTERNAL
Computes the sliding window histogram equalization of a matrix}
\usage{
  SlidingHistEqualC(mat, radius, fullmax = 4096)
}
\arguments{
  \item{mat}{the matrix to equalize}

  \item{radius}{the radius of the sliding window (total
  window size is a square window with sides 2*radius+1)}

  \item{fullmax}{the maximum value in the equalized image}
}
\value{
  The equalized matrix
}
\description{
  INTERNAL Computes the sliding window histogram
  equalization of a matrix
}
\details{
  Uses C code in slidinghistequalC.c
}

