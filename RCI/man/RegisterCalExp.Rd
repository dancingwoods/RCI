\name{RegisterCalExp}
\alias{RegisterCalExp}
\title{Removes in-plane motion effects using rigid body alignment of the image frames}
\usage{
  RegisterCalExp(calexp, refimg, channel = 1, upsamp = 2)
}
\arguments{
  \item{calexp}{a calexp object with a \$data field}

  \item{refimg}{a reference image to use for alignment.
  Should be the same size as the images in calexp\$data}

  \item{channel}{the channel to use for alignment}

  \item{upsamp}{the upsampling factor (this gives the
  sup-pixel precision of 1/upsamp)}
}
\value{
  a calexp object with a \$registration field.  The \$data
  in the returned object has been registered.  The
  \$registration field records the details of the estimated
  shifts.  \item{upsamp}{the usampling factor used}
  \item{refimg}{the reference image used} \item{mpars}{the
  estimated shifts. This is a matrix of size nframes-by-2}
}
\description{
  Removes in-plane motion effects using rigid body
  alignment of the image frames
}
\details{
  Registers the images in a calexp object by rigid body
  image alignment of the images in a particular channel to
  the reference image given.  The shifts for each image are
  estimated by maximum cross-correlation. Sub-pixel shifts
  are achieved using upsampling by the factor given
  (providing accuracy of 1/upsamp pixels).  Only one
  channel is used to register the data.  This should be the
  channel with the cleanest spatial information for the
  best results
}

