\name{ShiftFFT}
\alias{ShiftFFT}
\title{Shifts an image by the given amount, both translation and rotation}
\usage{
  ShiftFFT(img, pars, fdomain = FALSE, rotatefirst = FALSE)
}
\arguments{
  \item{img}{the image to shift}

  \item{pars}{a length-3 vector giving (x-translation,
  y-translation, rotation angle)}

  \item{fdomain}{is the image given in the Fourier domain?
  It will be returned in the same domain as given (passing
  in the Fourier domain is helpful to reduce superfluous
  transforms if performing additional operations in the
  Fourier domain).}

  \item{rotatefirst}{should rotation be performed before
  translation}
}
\value{
  the shifted image
}
\description{
  Shifts an image by the given amount, both translation and
  rotation
}
\details{
  Uses RotateFFT and TranslateFFT to compute result
}

