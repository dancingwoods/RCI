\name{RestrictMaskSize}
\alias{RestrictMaskSize}
\title{Removes masks from the database whose size is less than min pixels or greater
than max pixels}
\usage{
  RestrictMaskSize(db, minsize = NA, maxsize = NA)
}
\arguments{
  \item{db}{the mask database object}

  \item{min}{the minimum mask size (pixels) to retain}

  \item{max}{the maximum mask size (pixels) to retain}
}
\description{
  Removes masks from the database whose size is less than
  min pixels or greater than max pixels
}

