\name{GetSegmentation}
\alias{GetSegmentation}
\title{Get a matrix giving the segmentation stored in a mask database}
\usage{
  GetSegmentation(db, classids = c(2, 3), val = "id")
}
\arguments{
  \item{db}{the mask database}

  \item{classids}{a vector giving the ids of the classes to
  include in the returned segmentation}

  \item{val}{if "id", puts the mask ids in the ROI in the
  returned matrix, otherwise uses the class id}
}
\value{
  a data.frame with the extracted data
}
\description{
  Get a matrix giving the segmentation stored in a mask
  database
}

