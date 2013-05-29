\name{GrpRFCreate}
\alias{GrpRFCreate}
\title{Creates and trains a custom random forest classifier}
\usage{
  GrpRFCreate(labels, groups, data)
}
\arguments{
  \item{labels}{a vector giving the labels of the training
  data.  this should be a factor}

  \item{groups}{a vector giving the grouping of the
  training data.  Any data points the same group value are
  treated as members of the same group.}

  \item{data}{a data.frame giving the training data
  features}
}
\description{
  Creates and trains a custom random forest classifier
}

