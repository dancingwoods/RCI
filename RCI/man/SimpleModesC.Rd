\name{SimpleModesC}
\alias{SimpleModesC}
\title{INTERNAL
Finds the local maxima in an image}
\usage{
  SimpleModesC(img, min = 0)
}
\arguments{
  \item{img}{the image in which to find the local maxima}

  \item{min}{if this is set to 1, find local minima
  instead}
}
\value{
  matrix with 1 at the maxima and NA elsewhere
}
\description{
  INTERNAL Finds the local maxima in an image
}
\details{
  Uses C code in localmaxC.c
}

