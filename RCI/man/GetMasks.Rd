\name{GetMasks}
\alias{GetMasks}
\title{Returns a list of the masks in a database}
\usage{
  GetMasks(db)
}
\arguments{
  \item{db}{a database connection}
}
\value{
  a list of vectors, each vector specifying a mask.  The
  first element of each mask vector is the negative index
  of the mask.  The remaining elements of each vector are
  the indices of the mask pixels.
}
\description{
  Returns a list of the masks in a database
}

