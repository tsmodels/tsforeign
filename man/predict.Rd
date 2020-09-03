\name{predict}
\alias{predict}
\alias{predict.arima.estimate}
\alias{predict.bsts.estimate}
\title{
Prediction
}
\description{
Prediction function for class \dQuote{arima.estimate} and \dQuote{bsts.estimate}.
}
\usage{
\method{predict}{arima.estimate}(object, h = NULL, newxreg = NULL, nsim = 5000, 
forc_dates = NULL, bootstrap = FALSE, innov = NULL, ...)
\method{predict}{bsts.estimate}(object, h = 1, newxreg = NULL, 
forc_dates = NULL, last_state_means = colMeans(bsts_final_state(object)), 
posterior_means = NULL, innov = NULL, burn = NULL, ...)
}
\arguments{
\item{object}{
An object of the appropriate class which is then dispatched to the model
specific predict function.
}
\item{h}{
The forecast horizon.
}
\item{newxreg}{
The matrix of external regressors in the forecast horizon.
}
\item{nsim}{
The number of simulations to use for generating the simulated
predictive distribution. For the bsts model, this is equal 
to the number of MCMC samples generated during estimation, less
any burn-in draws.
}
\item{forc_dates}{
An optional vector of forecast dates equal to h. If NULL will
use the implied periodicity of the data to generate a regular
sequence of dates after the last available date in the data.
}
\item{last_state_means}{
This is an optional named vector which will re-center the last_state
means. It is best to call \link{bsts_final_state} function in order
to get the full matrix and make any changes to the state_means prior
to submitting to the function. It is required that the full vector
is provided and minimal checks other than length are performed.
}
\item{posterior_means}{
An optional vector of posterior parameter means which are then
use to re-center the posterior parameters. In the case of 
AR and regressors coefficients, these should represent the
means of the non-zero values (since both of these types of
parameters are based on a spike and slab prior). The posterior
distribution of the parameters can be obtain by calling
\link{bsts_posterior}.
}
\item{bootstrap}{
For the arima model whether to bootstrap the residuals.
}
\item{innov}{
For the bsts and arima models this is an optional vector
of uniform innovations which will be translated to regular
innovations using the appropriate distribution quantile
function and model standard deviation. The length of this 
vector should be equal to nsim x horizon for the arima
model and MCMC draws (less burn) x horizon for the
bsts model. Burn is equal to SuggestBurn(0.1, object$model)
for the bsts model.
}
\item{burn}{
An optional scalar denoting the numbers of draws to burn
from the posterior prior to prediction.
}
\item{\dots}{
Not currently used
}
}
\value{
An object which inherits class \dQuote{tsmodel.predict} with slots
for the simulated or posterior predictive distribution, the original 
series (as a zoo object), the original specification object and the mean forecast. 
The predictive distribution is inversed difference (if differencing > 0) and 
back transformed if lambda was not NULL in the original specification.
The innov argument for bsts is unlikely to be useful for ensembling, as
has been done for other models with a single source of error as a result
of multiple source of errors which disables our ability to infuse the
required dependence structure via this approach. 
}
\author{
Alexios Galanos
}
