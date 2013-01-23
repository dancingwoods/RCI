\name{ImageToCoordMat}
\alias{ImageToCoordMat}
\title{INTERNAL
Converts an image matrix to a matrix with coordinates and values in the columns}
\usage{
  ImageToCoordMat(img)
}
\arguments{
  \item{img}{the matrix to convert}
}
\value{
  A matrix of size npixels-by-3. The first coordinate is
  the row, the the column and the third the intensity.
}
\description{
  INTERNAL Converts an image matrix to a matrix with
  coordinates and values in the columns
}

