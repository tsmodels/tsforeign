\name{summary}
\alias{summary}
\alias{summary.bsts.estimate}
\alias{summary.arima.estimate}
\title{
Model Estimation Summary
}
\description{
Summary method.
}
\usage{
\method{summary}{bsts.estimate}(object, 
quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)
\method{summary}{arima.estimate}(object, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{bsts.estimate} or \dQuote{arima.estimate}.
}
\item{quantiles}{
A vector of quantiles to evaluate for each variable.
}
\item{\dots}{
Not currently used
}
}
\value{
A printout of the parameter summary, model type and some model metrics. For the
bsts model, the posterior distribution of parameters is fist converted into an 
\link{mcmc} object (coda package).
}
\author{
Alexios Galanos
}
