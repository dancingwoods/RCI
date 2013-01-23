\name{MatrixToSparse}
\alias{MatrixToSparse}
\title{INTERNAL
Converts a matrix mask into a sparse mask.  Assumes that the non-mask pixels
of the matrix are NA.}
\usage{
  MatrixToSparse(mat)
}
\arguments{
  \item{mat}{The mask as a matrix with NA in non-mask
  pixels}
}
\value{
  a vector of indices of the mask pixels
}
\description{
  INTERNAL Converts a matrix mask into a sparse mask.
  Assumes that the non-mask pixels of the matrix are NA.
}

