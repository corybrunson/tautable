% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vec.inv.R
\name{vec.inv}
\alias{vec.inv}
\title{Count the inversions in a vector.}
\usage{
vec.inv(vec)
}
\arguments{
\item{vec}{A numerical vector (theoretically, with no repeated values).}
}
\value{
The number of pairs of entries *i* and *j* of \code{vec} for which *i* > *j*.
}
\description{
Count the inversions in a vector.
}
\examples{
vec.inv(1:5)
vec.inv(5:1)
}

