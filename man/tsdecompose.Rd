\name{tsdecompose}
\alias{tsdecompose}
\alias{tsdecompose.bsts.estimate}
\alias{tsdecompose.bsts.predict}
\title{
Model Decomposition
}
\description{
Decomposes the estimated model or prediction into its component parts.
}
\usage{
\method{tsdecompose}{bsts.estimate}(object, ...)
\method{tsdecompose}{bsts.predict}(object, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{bsts.estimate} or \dQuote{bsts.predict}.
}
\item{\dots}{
Not currently used
}
}
\value{
A list of class \dQuote{tsmodel.distribution} representing the distribution of 
the state components of the model.
}
\author{
Alexios Galanos
}
