\name{plot.MTSpectrum}
\alias{plot.MTSpectrum}
\title{Plots a multitaper spectral estimate created by MultiTaperSpectrum}
\usage{
  plot.MTSpectrum(spect, maglog = TRUE, minfreq = 0,
    maxfreq = NULL, ...)
}
\arguments{
  \item{spect}{the multitaper spectrum object}

  \item{maglog}{should the magnitide be plotted on the log
  scale}

  \item{minfreq}{the minimum frequency to plot}

  \item{maxfreq}{the maximum frequency to plot (NULL for
  Nyquist frequency)}

  \item{...}{other graphical parameters}
}
\value{
  NULL
}
\description{
  Plots a multitaper spectral estimate created by
  MultiTaperSpectrum
}

