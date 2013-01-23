\name{GetInnerMasks}
\alias{GetInnerMasks}
\title{INTERNAL
Selects the masks from the given list that are contained in a region}
\usage{
  GetInnerMasks(framemat, masklist)
}
\arguments{
  \item{framemat}{a matrix of the same size as the masks in
  masklist with non-NA pixels specifying the region in
  which to find masks}

  \item{masklist}{a list, as returned by GetMasks, of
  sparse masks (vectors where the first element is the
  negative id of the mask and the other elements are the
  mask indices)}
}
\description{
  INTERNAL Selects the masks from the given list that are
  contained in a region
}
\details{
  Given a list of masks as returned by GetMasks and a
  matrix with a mask specifying a region, returns the masks
  in the masklist that are completely contained in the
  given region.
}

