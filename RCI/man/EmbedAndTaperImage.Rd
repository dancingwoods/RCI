\name{EmbedAndTaperImage}
\alias{EmbedAndTaperImage}
\title{Embeds an image in a larger matrix of 0's and tapers the image edges using a Hanning window}
\usage{
  EmbedAndTaperImage(img, taperamt, size = NULL,
    border = NULL)
}
\arguments{
  \item{img}{the image to embed and taper}

  \item{taperamt}{the width of the taper on the edges of
  the image.  Must be less than or equal to half the image
  width}

  \item{border}{the width of the border of 0's to add}
}
\value{
  an image that has been embedded and tapered
}
\description{
  Embeds an image in a larger matrix of 0's and tapers the
  image edges using a Hanning window
}
\details{
  uses size if given, else uses border, else doesn't embed
}

