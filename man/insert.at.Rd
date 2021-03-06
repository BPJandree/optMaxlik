% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optMaxlik.R
\name{insert.at}
\alias{insert.at}
\title{Mainly auxiliary. Inserts the content of an object into the supplied vector.}
\usage{
insert.at(a, pos, ...)
}
\arguments{
\item{a}{vector in which to insert a value}

\item{pos}{numeric indicating the position where to insert a value.}

\item{...}{objects to be inserted.}
}
\value{
vector
}
\description{
Function to insert values into a vector.
}
\examples{
a <- 1:10
insert.at(a, 5, 0)
}
\keyword{fixed}
\keyword{parameters}

