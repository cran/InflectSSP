% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Quantify.R
\name{Quantify}
\alias{Quantify}
\title{This function determines the median abundance value across the proteome for all experiments together}
\usage{
Quantify(Data_Normalized, NReps)
}
\arguments{
\item{Data_Normalized}{the normalized abundance data calculated in the Normalize function}

\item{NReps}{the number of replicates to be analyzed}
}
\value{
The median abundance data for all experiments at the proteome level
}
\description{
This function determines the median abundance value across the proteome for all experiments together
}
\examples{
\dontrun{
     Data_Quantified<-Quantify(Data_Normalized)}
}
