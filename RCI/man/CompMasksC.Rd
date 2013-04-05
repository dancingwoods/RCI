\name{CompMasksC}
\alias{CompMasksC}
\title{INTERNAL
Compute the overlap matrix between a set of masks}
\usage{
  CompMasksC(masklist)
}
\arguments{
  \item{masklist}{A list of sparse masks.  Each element of
  this list is a vector whose first element is the negative
  id of the mask and whose other elements are the sorted
  indices of the mask pixels.}
}
\value{
  a matrix whose elements are 0 or 1 giving the overlap
  relationships between the masks.  The masks are sorted in
  the matrix in the same order as they are given in
  masklist.
}
\description{
  INTERNAL Compute the overlap matrix between a set of
  masks
}
\details{
  Computes the overlap matrix between a set of masks using
  C code for efficiency.
}

