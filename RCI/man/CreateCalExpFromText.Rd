\name{CreateCalExpFromText}
\alias{CreateCalExpFromText}
\title{Convert a folder of text images to a calexp data object}
\usage{
  CreateCalExpFromText(name, imgdir, nchans = 2)
}
\arguments{
  \item{name}{a short name to identify this experiment}

  \item{imgdir}{a string giving the directory path for the
  directory containing the csv images}

  \item{nchans}{the number of channels that exist in the
  data}
}
\value{
  an object of class calexp \item{name}{the name passed in
  as an argument to this function} \item{data}{an array
  containing the image data, with dimensions
  nchans-nrows-ncols}
}
\description{
  Convert a folder of text images to a calexp data object
}
\details{
  This function Converts a directory of text files into a
  calexp data object in R. Assumes that the images are
  individual text files and that they are alphabetically in
  order by channel and then by time index.  The directory
  must contain only these image files. Each image must have
  the same dimensions, and there must be the same number of
  images for each channel.
}

