\name{fitted}
\alias{fitted}
\alias{fitted.bsts.estimate}
\alias{fitted.arima.estimate}
\title{
Model Fitted Values
}
\description{
Extract the fitted values from an estimated model.
}
\usage{
\method{fitted}{bsts.estimate}(object, distribution = FALSE, 
invdiff = TRUE, raw = FALSE, type = c("filtered", "smoothed"),...)
\method{fitted}{arima.estimate}(object, raw = FALSE,  ...)
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
differencing when calculating the fitted values.
}
\item{raw}{
If a Box Cox transformation was used, whether to return the
untransformed fitted values (raw).
}
\item{type}{
For the bsts model, the filtered are the one step ahead 
forecast values, whilst the smoothed represent those 
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
