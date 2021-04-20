\name{tsmetrics.bsts.estimate}
\alias{tsmetrics}
\alias{tsmetrics.bsts.estimate}
\alias{tsmetrics.arima.estimate}
\alias{tsmetrics.bsts.predict}
\alias{tsmetrics.arima.predict}
\title{
Performance Metrics
}
\description{
Performance metrics from an estimated model.
}
\usage{
\method{tsmetrics}{bsts.estimate}(object, ...)
\method{tsmetrics}{arima.estimate}(object, ...)
\method{tsmetrics}{bsts.predict}(object, actual, alpha = 0.1, ...)
\method{tsmetrics}{arima.predict}(object, actual, alpha = 0.1, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{arima.estimate} or \dQuote{bsts.estimate}.
}
\item{actual}{
The actual data matched to the dates of the forecasts.
}
\item{alpha}{
The coverage level for distributional forecast metrics.
}
\item{\dots}{
Optional arguments passed to the MASE function which includes 
\dQuote{frequency}, else will be read from the object (act as an
ovveride for instance when using fourier seasonality in xreg for arima).
}
}
\author{
Alexios Galanos
}
