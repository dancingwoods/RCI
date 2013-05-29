\name{SparseToMatrix}
\alias{SparseToMatrix}
\title{INTERNAL
Converts a sparse mask to a matrix mask}
\usage{
  SparseToMatrix(sm, ny = 126, nx = 126, background = NA)
}
\arguments{
  \item{sm}{the sparse representation of the mask (a vector
  whose positive values are the indices of the mask
  pixels)}

  \item{ny}{the number of rows of the matrix mask}

  \item{nx}{the number of columbs of the matrix mask}

  \item{background}{the value to put in the non-mask pixels
  of the matrix}
}
\value{
  A matrix of dimension (ny, nx) with 1's in the mask
  pixels and background elsewhere
}
\description{
  INTERNAL Converts a sparse mask to a matrix mask
}

