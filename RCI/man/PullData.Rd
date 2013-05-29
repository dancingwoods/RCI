\name{PullData}
\alias{PullData}
\title{Pulls data from a mask database into a data.frame that has one row for
each mask with the label, id, and feature values of that mask in the columns}
\usage{
  PullData(db, labeled = TRUE, group = TRUE)
}
\arguments{
  \item{db}{the mask database object}

  \item{labeled}{boolean, should the results be restricted
  to masks that have been labeled}
}
\value{
  a data.frame with the extracted data
}
\description{
  Pulls data from a mask database into a data.frame that
  has one row for each mask with the label, id, and feature
  values of that mask in the columns
}

