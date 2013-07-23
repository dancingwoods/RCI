\name{EvaluateSegmentation}
\alias{EvaluateSegmentation}
\title{Takes the current labels and segmentation in the given database and evaluates perfomance
of the segmenter against the hand labels.
For each group of labels cells, a false negative if no mask overlaps, false positive if segmented
mask overlaps with no group, marginal result if segmented cell overlaps with but is not of of the
labeled masks.}
\usage{
  EvaluateSegmentation(db)
}
\arguments{
  \item{db}{the database for which to evaluate the
  segmentation.  should have both a segmentation and hand
  labels.}
}
\value{
  a data.frame with a row for each class and columns
  specifying performance
}
\description{
  Takes the current labels and segmentation in the given
  database and evaluates perfomance of the segmenter
  against the hand labels. For each group of labels cells,
  a false negative if no mask overlaps, false positive if
  segmented mask overlaps with no group, marginal result if
  segmented cell overlaps with but is not of of the labeled
  masks.
}

