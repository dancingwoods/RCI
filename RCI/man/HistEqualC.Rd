\name{HistEqualC}
\alias{HistEqualC}
\title{INTERNAL
Computed the histogram equalization of a matrix.}
\usage{
  HistEqualC(mat, fullmax = 4096)
}
\arguments{
  \item{mat}{the matrix to equalize}

  \item{fullmax}{the range to equalize to}
}
\value{
  the equalized matrix
}
\description{
  INTERNAL Computed the histogram equalization of a matrix.
}
\details{
  Uses C code in histequalC.c
}

