\name{bsts_final_state}
\alias{bsts_final_state}
\title{
Extract the BSTS Final State Distribution
}
\description{
Extracts the final state distribution from the BSTS object.
}
\usage{
bsts_final_state(object)
}
\arguments{
\item{object}{
An object of class \dQuote{bsts.estimate}
}
}
\details{
This is the final smoothed state used as the initial state 
in forecasting.
}
\value{
A named matrix of draws x parameters.
}
\author{
Alexios Galanos
}