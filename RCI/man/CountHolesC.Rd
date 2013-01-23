\name{CountHolesC}
\alias{CountHolesC}
\title{INTERNAL
Counts the number of pixels not in a mask that are surrounded by at least 3 mask pixels}
\usage{
  CountHolesC(mask)
}
\arguments{
  \item{mask}{the mask in which to count holes.  NA or 0 in
  the background.}
}
\value{
  an integer giving the number of holes in the mask
}
\description{
  INTERNAL Counts the number of pixels not in a mask that
  are surrounded by at least 3 mask pixels
}
\details{
  Uses C code from the file countholesC.c
}

