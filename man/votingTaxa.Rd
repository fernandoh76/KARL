% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/votingTaxa.R
\name{votingTaxa}
\alias{votingTaxa}
\title{votingTaxa}
\usage{
votingTaxa(ids, nproc, plim = 3)
}
\arguments{
\item{ids}{Is a vector containing identifiers for misclassified genomes.}

\item{nproc}{Is the number of processors to perform calculations.}

\item{plim}{Is the maximum number of Principal Components taken for distance calculation.}
}
\description{
This function helps to re-assign misclassified instances in default classification models by calculating distances to neighboring points in the PCA space.
}
\examples{
load('evalActinobacteria.df.i3.Rdata')
actinoTprF<-evalActinobacteria.df.i3@perfectErrors$TprF
actino.vt7<-votingTaxa(ids=actinoTprF,nproc=10,plim=7)
}
\keyword{KARL}
\keyword{PCA}
\keyword{classification}
\keyword{enzymes}

