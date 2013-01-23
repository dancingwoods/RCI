\name{DbAddMask}
\alias{DbAddMask}
\title{Add a mask to a database}
\usage{
  DbAddMask(db, mask, source)
}
\arguments{
  \item{db}{a database connection object}

  \item{mask}{a matrix giving the mask to add to the
  database (T/F, 0/1, or NA/1)}

  \item{source}{a string giving the tag for the source of
  the mask}
}
\value{
  NULL
}
\description{
  Add a mask to a database
}
\details{
  Adds the given mask to the database.  If the mask is
  already in the database, increments the count for the
  source of the mask (or adds a new count for a new source)
}

