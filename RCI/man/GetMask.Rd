\name{GetMask}
\alias{GetMask}
\title{Return the requested mask from the specified database}
\usage{
  GetMask(db, id, sparse = T)
}
\arguments{
  \item{db}{a database connection}

  \item{id}{the id of the mask to return}

  \item{format}{"sparse" for a sparse mask in vector form,
  "matrix" for a matrix mask}
}
\value{
  either a vector giving the indices of the requested mask
  or a matrix version of the mask
}
\description{
  Return the requested mask from the specified database
}

