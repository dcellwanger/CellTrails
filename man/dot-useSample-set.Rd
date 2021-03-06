% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessor-methods.R
\docType{methods}
\name{.useSample<-}
\alias{.useSample<-}
\alias{useSample<-,SingleCellExperiment-method}
\title{SET trajectory samples indicator}
\usage{
.useSample(object) <- value
}
\arguments{
\item{object}{An object of class \code{SingleCellExperiment}}

\item{value}{A logical vector}
}
\value{
An updated object of class \code{SingleCellExperiment}
}
\description{
Sets indicator if sample was used for trajectory reconstruction.
}
\author{
Daniel C. Ellwanger
}
\keyword{internal}
