% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arima-model.R
\name{arima_modelspec}
\alias{arima_modelspec}
\title{Auto ARIMA specification}
\usage{
arima_modelspec(
  y,
  xreg = NULL,
  frequency = NULL,
  seasonal = FALSE,
  seasonal_type = "regular",
  seasonal_harmonics = NULL,
  transformation = "box-cox",
  lambda = NULL,
  lower = 0,
  upper = 1.5,
  ...
)
}
\arguments{
\item{y}{an xts vector.}

\item{xreg}{an xts matrix of external regressors.}

\item{frequency}{frequency of y (if using a seasonal model).}

\item{seasonal}{whether to include a seasonal component.}

\item{seasonal_type}{type of seasonality (regular or trigonometric).}

\item{seasonal_harmonics}{number of harmonics to include in the seasonal 
component when seasonal_type is trigonometric.}

\item{transformation}{a valid transformation for y from the \dQuote{tstransform} 
function in the \dQuote{tsaux} package (currently box-cox or logit are available).}

\item{lambda}{the Box Cox lambda. If not NULL, then either a numeric value or NA 
denoting automatic calculation.}

\item{lower}{lower bound for the transformation.}

\item{upper}{upper bound for the transformation.}

\item{...}{additional terms passed to the \code{\link{auto.arima}} function.}
}
\value{
An object of class \dQuote{arima.spec}.
}
\description{
Sets up an automatic ARIMA specification ready for estimation.
}
\note{
This is a wrapper to the auto.arima function from the forecast package.
}
