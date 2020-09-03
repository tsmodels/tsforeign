\name{bsts_modelspec}
\alias{bsts_modelspec}
\title{
BSTS Specification
}
\description{
Specifies an BSTS model prior to estimation.
}
\usage{
bsts_modelspec(y, xreg = NULL, frequency = NULL, differences = 0, 
level = TRUE, slope = TRUE, damped = FALSE, seasonal = FALSE, 
seasonal_frequency = 4, ar = FALSE, ar_max = 1, cycle = FALSE, 
cycle_frequency = NULL, cycle_names = NULL, seasonal_type = "regular", 
lambda = NULL, lambda_lower = 0, lambda_upper = 1, seasonal_harmonics = NULL, 
distribution = "gaussian", ...)
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
\item{differences}{
The number of differences to apply to the outcome
variable y (max of 2).
}
\item{level}{
Whether to include a level component (Local Level Model).
}
\item{slope}{
Whether to include a slope component (Local Linear Model).
}
\item{damped}{
Whether to include a damped trend (damped Local Linear
Model).
}
\item{seasonal}{
Whether to include a seasonal component.
}
\item{seasonal_frequency}{
A vector of seasonal frequencies.
}
\item{seasonal_type}{
The type of seasonality (regular or trigonometric).
}
\item{seasonal_harmonics}{
A vector of length equal to the length of the seasonal
frequency denoting the harmonics to use when seasonal 
type is trigonometric.
}
\item{ar}{
Whether to include a sparse AR component.
}
\item{ar_max}{
The number of lags for the AR component.
}
\item{cycle}{
Whether to include a cyclical component.
}
\item{cycle_frequency}{
The number of periods in a cycle. This can
be a vector in which case multiple cycles
are included.
}
\item{cycle_names}{
An optional vector of cycle names.
}
\item{lambda}{
The Box Cox power parameter (lambda). If NA will estimate this using 
the method of Guerrero.
}
\item{lambda_lower}{
Lower bound on Box Cox lambda transformation.
}
\item{lambda_upper}{
Upper bound on Box Cox lambda transformation.
}
\item{distribution}{
Valid choices are currently only \dQuote{gaussian}.
}
\item{\dots}{
Not currently used/
}
}
\value{
An object of class \dQuote{bsts.spec}
}
\author{
Alexios Galanos
}
\note{
This is a wrapper to part of the functionality of the bsts package.
Once an object is estimated, all other methods are implemented 
locally (including prediction).
}