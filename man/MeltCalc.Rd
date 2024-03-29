% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MeltCalc.R
\name{MeltCalc}
\alias{MeltCalc}
\title{This function determines melt shifts for all proteins that meet quality criteria and also determines the melt shift p-values}
\usage{
MeltCalc(
  Directory,
  Data_CurveFit2_Complete_Unique,
  CurveRsq,
  PValMelt,
  MeltLimit,
  PValMeltFDR
)
}
\arguments{
\item{Directory}{the directory data is saved to}

\item{Data_CurveFit2_Complete_Unique}{the curve fit data from the CurveFit2 function}

\item{CurveRsq}{the criteria for melt curve p-values}

\item{PValMelt}{the criteria for the melt shift p-values}

\item{MeltLimit}{the melt shift temperature limit used for determining which proteins are significant}

\item{PValMeltFDR}{Whether or not the FDR correction for pvalue is used in designation of melts of interest}
}
\value{
Proteins melt shifts
}
\description{
This function determines melt shifts for all proteins that meet quality criteria and also determines the melt shift p-values
}
\examples{
\dontrun{
     Data_Melts<-MeltCalc(Directory,Data_CurveFit2_Complete_Unique,
     CurveRsq,PValMelt,MeltLimit,PValMeltFDR)}
}
