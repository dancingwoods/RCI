\name{ViewCI}
\alias{ViewCI}
\title{Opens the GUI viewer to manipulate the segmentation process.}
\usage{
  ViewCI(db, cf = NULL)
}
\arguments{
  \item{db}{if specified, the viewer opens with the given
  dbController (looking in the directories stored in that
  object)}

  \item{cf}{a classifier object.  Needed to allow redoing
  segmentation after correcting labels}
}
\value{
  NULL
}
\description{
  Opens the GUI viewer to manipulate the segmentation
  process.
}

