\name{SetMaskLabel}
\alias{SetMaskLabel}
\title{Sets the label field for a particular mask in a mask database}
\usage{
  SetMaskLabel(db, id, label)
}
\arguments{
  \item{db}{a database connection}

  \item{id}{the id of the mask to label}

  \item{label}{the label to assign to the mask (0=unknown,
  1=cell, 2=not cell)}
}
\value{
  NULL
}
\description{
  Sets the label field for a particular mask in a mask
  database
}

