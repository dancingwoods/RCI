\name{CompMaskC}
\alias{CompMaskC}
\title{INTERNAL
Compute the overlap of a single mask with a list of masks}
\usage{
  CompMaskC(mask, masklist)
}
\arguments{
  \item{mask}{A single mask formatted as a vector of mask
  indices with or without the negative id as the first
  element of the vector.}

  \item{masklist}{A list of sparse masks as returned by
  GetMasks.  Each element of this list is a vector whose
  first element is the negative id of the mask and whose
  other elements are the sorted indices of the mask
  pixels.}
}
\value{
  a vector whose elements are 0 or 1 giving the overlap
  relationships between the masks.  The values are sorted
  in the vector in the same order as they are given in
  masklist.
}
\description{
  INTERNAL Compute the overlap of a single mask with a list
  of masks
}
\details{
  Computes the overlap matrix between a mask and a list of
  other masks using C code for efficiency
}

