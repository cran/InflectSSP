% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Correction.R
\name{Correction}
\alias{Correction}
\title{This function corrects the normalized abundance of each protein using a correction constant that is calculated in this function. The correction constant is determined using the difference between actual and predicted fit at the proteome level.}
\usage{
Correction(PSM, UP, Data_CurveFit1Parameters, Data_Normalized, Data_Quantified)
}
\arguments{
\item{PSM}{the number of peptide spectrum matches that are deemed acceptable for reporting}

\item{UP}{the number of unique peptides for a protein that are deemed acceptable for reporting}

\item{Data_CurveFit1Parameters}{the parameters determined from Curve Fit 1 operation for proteome melts}

\item{Data_Normalized}{the normalized abundance data for each protein determined in the Normalize function.}

\item{Data_Quantified}{the median normalized abundance data at the proteome level}
}
\value{
the corrected and normalized abundance data for each protein
}
\description{
This function corrects the normalized abundance of each protein using a correction constant that is calculated in this function. The correction constant is determined using the difference between actual and predicted fit at the proteome level.
}
\examples{
\dontrun{
Data_Corrected<-Correction(PSM,UP,Data_CurveFit1Parameters,
Data_Normalized,Data_Quantified)
}
}
