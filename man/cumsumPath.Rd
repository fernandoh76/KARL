% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cumsumPath.R
\name{cumsumPath}
\alias{cumsumPath}
\title{cumsumPath}
\usage{
cumsumPath(rank, taxa, pid, ecol)
}
\arguments{
\item{rank}{Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.}

\item{taxa}{A vector of taxa names.}

\item{pid}{The pathway id.}

\item{ecol}{A vector of colors of same length as taxa.}
}
\description{
This function plots a cumulative curve for pathway completeness.
}
\examples{
taxcum<-c('Actinobacteria','Tenericutes','Firmicutes','Proteobacteria','Bacteroidetes')
taxcol<-c('red','green','blue','violet','orange')
cumsumPath(rank='phylum',taxa=taxcum,pid='00550',ecol=taxcol)
}
\keyword{KARL}
\keyword{enzymes}
\keyword{pathways}

