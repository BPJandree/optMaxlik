% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optMaxlik.R
\name{aicc}
\alias{aicc}
\title{Corrected AIC}
\usage{
aicc(x, n = Inf, penalty = 0)
}
\arguments{
\item{x}{maxLik object}

\item{n}{number of observations for correction. Defaults to Inf, which coincides with the AIC.}

\item{penalty}{numeric that should be added to the log likelihood to adjust when using penalized likelihood or regularized estimators.}
}
\value{
numeric
}
\description{
Compute the corrected AIC from a maxLik object.
}
\keyword{AICc}

