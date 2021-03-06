% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster-methods.R
\name{.findStates_def}
\alias{.findStates_def}
\title{DEF: Find states}
\usage{
.findStates_def(
  X,
  ordi,
  link.method = "ward.D2",
  min.size,
  max.pval = 1e-04,
  min.fc = 2,
  min.g = 5,
  reverse = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{X}{Expression matrix}

\item{ordi}{Ordination of samples}

\item{link.method}{Linkage criteria}

\item{min.size}{Min size of inital clusters}

\item{max.pval}{Pval threshold}

\item{min.fc}{Fold-change threshold}

\item{reverse}{Reverse order}

\item{verbose}{For debug}

\item{ming.g}{Feature cound threshold}
}
\description{
For details see \code{findStates}
}
\author{
Daniel C. Ellwanger
}
\keyword{internal}
