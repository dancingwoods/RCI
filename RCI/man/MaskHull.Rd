\name{MaskHull}
\alias{MaskHull}
\title{INTERNAL
Computes the convex hull of a mask}
\usage{
  MaskHull(sparsemask, nr, nc)
}
\arguments{
  \item{sparsemask}{the sparse mask (vector) for which to
  find the convex hull.}

  \item{nr}{the number of rows in the image}

  \item{nc}{the number of columns in the image}
}
\value{
  a sparse representation of the convex hull of the given
  mask
}
\description{
  INTERNAL Computes the convex hull of a mask
}
\details{
  FIXME: there's the issue that
  maphull(maphull(x))!=maphull(x), but using this anyway
}

