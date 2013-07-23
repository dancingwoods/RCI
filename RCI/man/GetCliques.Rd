\name{GetCliques}
\alias{GetCliques}
\title{Assigns each of the specified masks to a clique to use when
solving the MWIS segmentation}
\usage{
  GetCliques(db, ids, minid = 1)
}
\arguments{
  \item{db}{the mask database}

  \item{ids}{the vector of ids of the masks to assign to
  cliques}

  \item{minid}{the minimum value of the clique ids to
  return}
}
\value{
  a vector giving the clique id for each of the specified
  masks
}
\description{
  Assigns each of the specified masks to a clique to use
  when solving the MWIS segmentation
}

