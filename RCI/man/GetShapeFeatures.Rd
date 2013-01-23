\name{GetShapeFeatures}
\alias{GetShapeFeatures}
\title{INTERNAL
Computes the features related to just the shape of masks, adding them to the database}
\usage{
  GetShapeFeatures(db, mids = NULL)
}
\arguments{
  \item{db}{a database connection}

  \item{mids}{an optional vector of the mask ids for which
  to extract features}
}
\value{
  NULL
}
\description{
  INTERNAL Computes the features related to just the shape
  of masks, adding them to the database
}
\details{
  Computes features of all masks in the database or a list
  of masks specified by id.  The features computed are
  currently:\cr npixels - the number of pixels in a mask
  \cr nholes - the number of non-mask pixels that are
  surrounded by at least 3 mask pixels\cr bboxratio - the
  ratio of the area of the mask's bounding box and the
  number of pixels\cr in the mask hullratio - the ratio of
  the area of the mask's convex hull and the number of
  pixels in the mask\cr
}

