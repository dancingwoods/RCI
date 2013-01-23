\name{ConMaskDb}
\alias{ConMaskDb}
\title{Connects to an experiment's mask database}
\usage{
  ConMaskDb(path)
}
\arguments{
  \item{path}{the path to the SQLite database to connect
  to}
}
\value{
  a connection object as returned by dbConnect in the DBI
  package
}
\description{
  Connects to an experiment's mask database
}

