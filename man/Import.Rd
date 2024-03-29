% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Import.R
\name{Import}
\alias{Import}
\title{This function imports data that will be analyzed in downstream functions.}
\usage{
Import(NControl, NCondition, Directory)
}
\arguments{
\item{NControl}{the number of Control replicate experiments that are to be analyzed}

\item{NCondition}{the number of Condition replicate experiments that are to be analyzed}

\item{Directory}{the directory where the source data files to be analyzed are saved. This is also the location where the results will be saved.}
}
\value{
Imported data from all experiments
}
\description{
This function imports data that will be analyzed in downstream functions.
}
\examples{
\dontrun{
Data_Imported<-Import(NControl,NCondition,Directory)
}
}
