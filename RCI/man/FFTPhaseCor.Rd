\name{FFTPhaseCor}
\alias{FFTPhaseCor}
\title{INTERNAL
Computes sub-pixel shifts values using phase correlation (FFT implementation)}
\usage{
  FFTPhaseCor(img1, img2, upsamp = 2, taper = TRUE,
    cortaper = TRUE, subpixel = "gauss", subrad = 3)
}
\arguments{
  \item{img1}{matrix giving the first image (the
  reference)}

  \item{img2}{matrix giving the second image (to be
  shifted)}

  \item{upsamp}{the factor by which the fft matrix should
  be expanded}

  \item{taper}{boolean, should the images be tapered before
  alignment}

  \item{cortaper}{boolean, should the normalized
  cross-spectrum be tapered before being (inverse)
  transformed}

  \item{subpixel}{'none' for no additional subpixel
  fitting, 'gauss' for Gaussian fit, 'poc' for poc function
  fitting}

  \item{subrad}{the radius of the submatrix used to compute
  the subpixel fits}
}
\value{
  a vector of length 2 giving the magnitude of the estimted
  x and y shift returns NA in the case of improper input
}
\description{
  INTERNAL Computes sub-pixel shifts values using phase
  correlation (FFT implementation)
}
\details{
  Computes the sub-pixel shifts by computing the upsampled
  phase correlation between the two images and finding the
  maximum.  If the parameter gausfit is TRUE, then a
  gaussian is fit around the peak of the phase correlation
  function to get additional sub-pixel shift information.
  This is on top of any upsampling
}

