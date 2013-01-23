\name{HillClimbC}
\alias{HillClimbC}
\title{INTERNAL
Perform hill climbing on a matrix starting from a given point and returning the
local maxima that is reached.}
\usage{
  HillClimbC(y, x, img)
}
\arguments{
  \item{y}{Starting row}

  \item{x}{Starting column}

  \item{img}{The matrix on which to performt the
  hillclimbing}
}
\value{
  a vector of 2 numbers giving the coordinates of the peak
  found by hillclimbing
}
\description{
  INTERNAL Perform hill climbing on a matrix starting from
  a given point and returning the local maxima that is
  reached.
}
\details{
  Uses C code in hillclimbC.c
}

