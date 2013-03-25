\name{ShiftFFTVector}
\alias{ShiftFFTVector}
\title{INTERNAL
Shifts a vector by the specified amount using FFT phase shift,
but assuming the Fourier transform has already been performed.}
\usage{
  ShiftFFTVector(vec, amt)
}
\arguments{
  \item{vec}{the vector to shift}

  \item{amt}{the amount to shift}
}
\value{
  the circularly shifted vector
}
\description{
  INTERNAL Shifts a vector by the specified amount using
  FFT phase shift, but assuming the Fourier transform has
  already been performed.
}

