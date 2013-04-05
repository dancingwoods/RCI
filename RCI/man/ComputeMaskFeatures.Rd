\name{ComputeMaskFeatures}
\alias{ComputeMaskFeatures}
\title{Computed features of masks that are currently in the mask database.
This is sort of hack - should be more general and allow users to
specify feature functions}
\usage{
  ComputeMaskFeatures(db, calexp, feature)
}
\arguments{
  \item{db}{the mask database object}

  \item{calexp}{the calcium experiment object}

  \item{feature}{string indicating which feature to
  compute}
}
\description{
  Computed features of masks that are currently in the mask
  database. This is sort of hack - should be more general
  and allow users to specify feature functions
}
\details{
  Current feature options: "size" - mask size in pixels
}

