\name{MaskHull}
\alias{MaskHull}
\title{INTERNAL
Computes the convex hull of a mask}
\usage{
  MaskHull(mask)
}
\arguments{
  \item{mask}{the mask for which to find the convex hull.
  Background pixels should be NA}
}
\value{
  a matrix with 1's on the convex hull of the mask and NA
  in the background
}
\description{
  INTERNAL Computes the convex hull of a mask
}
\details{
  FIXME: there's the issue that
  maphull(maphull(x))!=maphull(x), but using this anyway
}

