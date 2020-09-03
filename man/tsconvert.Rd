\name{tsconvert}
\alias{tsconvert}
\alias{tsconvert.bsts.estimate}
\title{Convert an estimated model to DLM}
\description{
Converts an estimated BSTS object to a DLM object.
}
\usage{
\method{tsconvert}{bsts.estimate}(object, to = "dlm", draw = "mean", 
burn = SuggestBurn(0.1, object$model), ...)
}
\arguments{
\item{object}{
An object class \dQuote{bsts.estimate}.
}
\item{to}{
Conversion to. Currently only DLM supported, but KFAS may be added in the
future.
}
\item{draw}{
Either an integer or \dQuote{mean}, after eliminating any burn draws.
}
\item{burn}{
The draws to burn.
}
\item{\dots}{
Nothing here.
}
}
\author{
Alexios Galanos
}
