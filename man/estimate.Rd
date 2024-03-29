% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arima-model.R, R/bsts-model.R
\name{estimate.arima.spec}
\alias{estimate.arima.spec}
\alias{estimate}
\alias{estimate.bsts.spec}
\title{Model Estimation}
\usage{
\method{estimate}{arima.spec}(object, ...)

\method{estimate}{bsts.spec}(
  object,
  n_iter = 5000,
  timeout.seconds = Inf,
  bma.method = "SSVS",
  trace = TRUE,
  ...
)
}
\arguments{
\item{object}{an object of class \dQuote{arima.spec} or \dQuote{bsts.spec}.}

\item{...}{for the bsts model, additional arguments passed to the underlying 
estimation functions in the BSTS package.}

\item{n_iter}{mcmc draws for the bsts model.}

\item{timeout.seconds}{timeout passed to the BstsOptions function.}

\item{bma.method}{Bayesian Model Averaging method in the presence of 
regressors, passed to the BstsOptions function.}

\item{trace}{whether to show the MCMC iterations.}
}
\value{
An object of class \dQuote{arima.estimate} or \dQuote{bsts.estimate}.
}
\description{
Estimates a model given a specification object using
maximum likelihood.
}
