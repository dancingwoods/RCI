\name{PredictExperiment}
\alias{PredictExperiment}
\title{Predicts an experiment using the classifier and MWIS
TODO - this is currently a hack with a heuristic to find cliques.  should really find connected
components and solve the MWIS}
\usage{
  PredictExperiment(classifier, db, thresh = NULL)
}
\arguments{
  \item{classifier}{the classifier}

  \item{dn}{the database to predict}

  \item{enforcelabels}{boolen, should the final
  segmentation be forced to correspond to the hand labels}
}
\value{
  null, modifies the database
}
\description{
  Predicts an experiment using the classifier and MWIS TODO
  - this is currently a hack with a heuristic to find
  cliques.  should really find connected components and
  solve the MWIS
}

