\name{LogSeq}
\alias{LogSeq}
\title{Generate a sequence with log scale jumps}
\usage{
  LogSeq(minv, maxv, length.out)
}
\arguments{
  \item{minv}{the value at the low end of the scale.  This
  is set to 0.01 if it is <=0}

  \item{maxv}{the value at the high end of the scale}

  \item{length.out}{the number of values to return}
}
\value{
  a vector with the requested sequence
}
\description{
  Generate a sequence with log scale jumps
}

