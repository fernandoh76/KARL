% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluateSVM.R
\name{evaluateSVM}
\alias{evaluateSVM}
\title{evaluateSVM}
\usage{
evaluateSVM(idata = alldata, rank, tax, folds = 10, iters = 3,
  parallel = T)
}
\arguments{
\item{rank}{Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.}

\item{tax}{The taxon name.}

\item{folds}{Is the number of folds in cross-validation.}

\item{iters}{Is the number of iterations in the repeated cross-validation.}

\item{parallel}{Is a logical indication if doing parallelization or sequential computing. Number of cores will be equal to number of folds. calculation.}
}
\description{
This function performs a m-repeated n-fold cross-validation and stores detailed results for each genome in a class called evalModel.
}
\examples{
evalCyanobacteria.df.i3<-evaluateSVM(rank='phylum',tax='Cyanobacteria',folds=10,parallel=T,iters=3)
}
\keyword{KARL}
\keyword{SVM}
\keyword{cross-validation}
\keyword{evaluate}
\keyword{model}

