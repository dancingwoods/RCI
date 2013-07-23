\name{EvaluateConfidence}
\alias{EvaluateConfidence}
\title{Computes several confidence measures on the segmentation in a database these include
- number of masks returned by the classifier
- min, max, and mean probability assigned by the classifier}
\usage{
  EvaluateConfidence(db, class)
}
\arguments{
  \item{db}{the database for which to evaluate the
  segmentation confidence.}

  \item{class}{the classifier used to generate the
  segmentation}
}
\value{
  a data.frame with a row for each ROI and columns
  specifying confidence measures
}
\description{
  Computes several confidence measures on the segmentation
  in a database these include - number of masks returned by
  the classifier - min, max, and mean probability assigned
  by the classifier
}

