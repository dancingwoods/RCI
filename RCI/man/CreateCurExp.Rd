\name{CreateCurExp}
\alias{CreateCurExp}
\title{INTERNAL
Create an object to store information about the currently selected experiment.}
\usage{
  CreateCurExp()
}
\value{
  A list with fields to store information about the
  experiment (initially empty) \item{name}{the experiment
  name} \item{data}{if loaded, the data object for this
  object} \item{db}{the database connection for this
  experiment} \item{nmasks}{the number of candidate masks
  in the database for this experiment} \item{features}{the
  tags of the features that exist in this database}
  \item{sources}{the tage for the mask sources present in
  this database} \item{selmat}{a matrix where the first
  column is the ID of the mask and the second column gives
  the annotation for the mask} \item{sms}{the list of
  sparse masks for the experiment retreived from the
  database with GetMasks()} \item{mimg1}{the matrix giving
  the mean image for channel 1} \item{mimg2}{the matrix
  giving the mean image for channel 2} \item{nx}{the number
  of columns in the images for this experiment}
  \item{ny}{the number of rows in the images}
}
\description{
  INTERNAL Create an object to store information about the
  currently selected experiment.
}
\details{
  Creates a list that stores information about the
  currently selected experiment.
}

