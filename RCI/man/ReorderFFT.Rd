\name{ReorderFFT}
\alias{ReorderFFT}
\title{INTERNAL
Reorders the matrix returned by fft}
\usage{
  ReorderFFT(mat, inverse = F)
}
\arguments{
  \item{mat}{a matrix of values to reorder}

  \item{inverse}{if true, takes reordered matrix and
  returns to order expected by fft. if false, takes matrix
  from fft and reorders it}
}
\value{
  the reordered matrix
}
\description{
  INTERNAL Reorders the matrix returned by fft
}
\details{
  Reorders the matrix returned by the R function fft.  The
  R function returns the coefficients from
  low-to-high-to-low frequencies in both dimensions.  The
  reordering puts the low frequencies in the center of the
  matrix so that the coefficients go from
  high-to-low-to-high in each dimension
}

