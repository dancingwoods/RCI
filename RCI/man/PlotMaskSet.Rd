\name{PlotMaskSet}
\alias{PlotMaskSet}
\title{Plots sets of masks over an already plotted image}
\usage{
  PlotMaskSet(mask, alpha = 0.5, rgb = NULL, ...)
}
\arguments{
  \item{mask}{the specification of the mask, unique values
  for each mask set, and 0 or NA in background}

  \item{alpha}{the alpha transparency value of the mask
  (between 0 and 1)}

  \item{...}{additional graphing parameters}
}
\value{
  NULL
}
\description{
  Plots sets of masks over an already plotted image
}
\details{
  Given a matrix with unique integers for each mask set,
  overplot each mask set in a different color (randomly
  chosen)
}

