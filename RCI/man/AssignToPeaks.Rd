\name{AssignToPeaks}
\alias{AssignToPeaks}
\title{Assigns the non-zero pixels of 'region' to one of the maxima of the iamge by hillclimbing
on image}
\usage{
  AssignToPeaks(region, image, restrict = T)
}
\arguments{
  \item{region}{a matrix with 1 in the regions to be
  assigned and 0 elsewhere}

  \item{image}{the image matrix}

  \item{restrict}{boolean.  should the hill-climbing be
  restricted to a path entirely within region}
}
\value{
  a matrix with unique integers in the pixels of region
  corresponding to each local maxima
}
\description{
  Assigns the non-zero pixels of 'region' to one of the
  maxima of the iamge by hillclimbing on image
}

