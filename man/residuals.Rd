\name{residuals}
\alias{residuals}
\alias{residuals.bsts.estimate}
\alias{residuals.arima.estimate}
\title{
Model Residuals
}
\description{
Extract the residual values from an estimated model.
}
\usage{
\method{residuals}{bsts.estimate}(object, distribution = FALSE, 
invdiff = TRUE, standardize = FALSE, raw = FALSE, 
type = c("filtered", "smoothed"), ...)
\method{residuals}{arima.estimate}(object, raw = FALSE, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{bsts.estimate} or \dQuote{arima.estimate}.
}
\item{distribution}{
Whether to return the full distribution for the bsts model.
}
\item{invdiff}{
If a model in differences was estimated, whether to inverse the
differencing when calculating the residual values.
}
\item{raw}{
If a Box Cox transformation was used, whether to return the
untransformed errors (raw).
}
\item{standardize}{
Whether to scale the errors by the model standard deviation.
}
\item{type}{
For the bsts model, the filtered are the one step ahead 
forecast errors, whilst the smoothed represent those 
values using the whole information set.
}
\item{\dots}{
Not currently used
}
}
\value{
For the distribution choice, an object of class \dQuote{tsmodel.distribution},
else an xts vector.
}
\author{
Alexios Galanos
}
