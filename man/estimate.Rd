\name{estimate}
\alias{estimate}
\alias{estimate.arima.spec}
\alias{estimate.bsts.spec}
\title{Estimate Model}
\description{
Estimate a model given a specification object using either 
maximum likelihood or MCMC (for BSTS).
}
\usage{
\method{estimate}{arima.spec}(object, ...)
\method{estimate}{bsts.spec}(object, n_iter = 5000, 
timeout.seconds = Inf, bma.method = "SSVS", trace = TRUE,...)
}
\arguments{
\item{object}{
An object of the appropriate class dispatched to the model's estimation
function.
}
\item{n_iter}{
The MCMC draws for the bsts model.
}
\item{timeout.seconds}{
The time out passed to the BstsOptions function.
}
\item{bma.method}{
The Bayesian Model Averaging method in the presence of regressors, 
passed to the BstsOptions function.
}
\item{trace}{
Whether to show the MCMC iterations.
}
\item{\dots}{
For the bsts model additional arguments passed to the
underlying estimation functions in the BSTS package.
}
}
\author{
Alexios Galanos
}
