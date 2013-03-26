\name{LoGMasks}
\alias{LoGMasks}
\title{Generates a set of masks using the Laplacian of Gaussian technique
for the given scale and kernel size}
\usage{
  LoGMasks(image, scale, ksize = 15, sparse = T)
}
\arguments{
  \item{image}{the image to use to generate masks}

  \item{scale}{the scale of the kernel to use}

  \item{kside}{the size of the kernel}

  \item{sparse}{boolean, should the function return sparse
  masks instead of a matrix for plotting}
}
\value{
  a matrix with unique integers at mask locations and 0 in
  the background, or a a sparse masks object
}
\description{
  Generates a set of masks using the Laplacian of Gaussian
  technique for the given scale and kernel size
}

