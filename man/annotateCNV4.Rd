% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnv_caller_functions_packaged.R
\name{annotateCNV4}
\alias{annotateCNV4}
\title{annotateCNV4}
\usage{
annotateCNV4(
  cnvResults,
  saveOutput = TRUE,
  maxClust2 = 4,
  outputSuffix = "_1",
  sdCNV = 0.6,
  filterResults = TRUE,
  filterRange = 0.8,
  minAlteredCells = 40
)
}
\arguments{
\item{cnvResults}{output of filtered CNVs (set of lists)}

\item{saveOutput}{save output to a file on disk}

\item{maxClust2}{maximum # of clusters}

\item{outputSuffix}{output suffix to save files}

\item{sdCNV}{sd for CNV copy number estimation}

\item{filterResults}{filter results to remove likely unaltered chromosomes}

\item{filterRange}{copy number range minimum}

\item{minAlteredCells}{filtering parameter - remove all alterations with a smallest group less than X cells (default is 40)}
}
\description{
Annotates the filtered CNV calls and saves the results, and estimates absolute copy numbers
}
\examples{
annotateCNV4(cleanCNV,saveOutput=TRUE,outputSuffix="_nolog",sdCNV=0.5)
}
\keyword{CNV}
\keyword{output}
