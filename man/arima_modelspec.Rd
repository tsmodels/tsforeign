\name{arima_modelspec}
\alias{arima_modelspec}
\title{
Auto ARIMA specification
}
\description{
Setups an auto.arima specification.
}
\usage{
arima_modelspec(y, xreg = NULL, frequency = NULL, seasonal = FALSE, 
seasonal_type = "regular", lambda = NULL, seasonal_harmonics = NULL, 
lambda_lower = 0, lambda_upper = 1.5, ...)
}
\arguments{
\item{y}{
An xts vector.
}
\item{xreg}{
An optional xts matrix of regressors (pre-lagged).
}
\item{frequency}{
The frequency of y (if using a seasonal model).
}
\item{seasonal}{
Whether to include a seasonal component.
}
\item{seasonal_type}{
The type of seasonality (regular or trigonometric).
}
\item{lambda}{
The Box Cox power parameter (lambda). If NA will 
estimate this using the method of Guerrero.
}
\item{lambda_lower}{
Lower bound on Box Cox lambda transformation.
}
\item{lambda_upper}{
Upper bound on Box Cox lambda transformation.
}
\item{seasonal_harmonics}{
The number of harmonics to include in the seasonal
component when seasonal_type is trigonometric.
}
\item{\dots}{
Additional terms passed to the auto.arima function.
}
}
\value{
An object of class \dQuote{arima.spec}.
}
\author{
Alexios Galanos
}
\note{
This is a wrapper to the auto.arima function from the
forecast package.
}

