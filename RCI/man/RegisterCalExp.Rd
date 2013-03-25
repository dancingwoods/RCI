\name{RegisterCalExp}
\alias{RegisterCalExp}
\title{Removes in-plane motion effects using rigid body alignment of the image frames}
\usage{
  RegisterCalExp(calexp, refimg, channel = 1,
    bigsize = c(256, 256))
}
\arguments{
  \item{calexp}{a calexp object with a \$data field}

  \item{refimg}{a reference image to use for alignment.
  Should be the same size as the images in calexp\$data}

  \item{channel}{the channel to use for alignment
  (typically the structual channel)}

  \item{upsamp}{the upsampling factor (this gives the
  sup-pixel precision of 1/upsamp)}
}
\value{
  a calexp object with a \$registration field.  The \$data
  in the returned object has been registered.  The
  \$registration field records the details of the estimated
  shifts.  \item{refimg}{the reference image used}
  \item{mpars}{the estimated shifts. This is a matrix of
  size nframes-by-2}
}
\description{
  Removes in-plane motion effects using rigid body
  alignment of the image frames
}
\details{
  Registers the images in a calexp object by rigid body
  image alignment of the images in a particular channel to
  the reference image given.  Initial translation
  parameters are estimated using Phase-Only correlation.
  The parameters are then optimized using Nelder-Mead
  optimization of the mean squared error between the
  images.
}

