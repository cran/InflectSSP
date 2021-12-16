#' This function is the primary function that calls other functions in the program.
#' @param Directory the directory where the source data files to be analyzed are saved. This is also the location where the results will be saved.
#' @param NControl the number of Control replicate experiments that are to be analyzed
#' @param NCondition the number of Condition replicate experiments that are to be analyzed
#' @param PSM the number of peptide spectrum matches that are deemed acceptable for reporting
#' @param UP the number of unique peptides for a protein that are deemed acceptable for reporting
#' @param CurveRsq Coefficient of determination criteria for melt curves
#' @param PValMelt p-value criteria for melt shifts
#' @param PValMeltFDR Whether or not the FDR correction for pvalue is used in designation of melts of interest
#' @param RunSTRING whether or not the STRING function will be run or not in the analysis
#' @param STRINGScore the score to be used in the STRING analysis
#' @param RunPANTHER whether or not the PANTHER function will be run or not in the analysis
#' @param Species species number for bioinformatics search
#' @param PANTHERpvalue p-value for PANTHER analysis
#' @param MeltLimit the melt shift temperature limit used for determining which proteins to report as significant
#' @importFrom readxl read_xlsx
#' @importFrom xlsx write.xlsx
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom tidyr spread
#' @importFrom grDevices dev.off
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics axis
#' @importFrom graphics legend
#' @importFrom plotrix addtable2plot
#' @importFrom data.table data.table
#' @importFrom grDevices pdf
#' @importFrom network network
#' @export
#' @return the proteins that have significant melt shifts from an experiment
#' @examples
#' \dontrun{
#'      Directory<-'/Users/Einstein'
#'      NControl<-2
#'      NCondition<-3
#'      PSM<-2
#'      UP<-3
#'      CurveRsq<-.95
#'      PValMelt<-0.05
#'      PValMeltFDR<-"No"
#'      MeltLimit<-3
#'      RunSTRING<-"Yes"
#'      STRINGScore<-0.99
#'      RunPANTHER<-"Yes"
#'      Species<-9606
#'      PANTHERpvalue<-0.05
#'      InflectSSP(Directory,NControl,
#'      NCondition,PSM,UP,CurveRsq,PValMelt,PValMeltFDR,
#'      MeltLimit,RunSTRING,STRINGScore,RunPANTHER,
#'      Species,PANTHERpvalue)
#'      }

InflectSSP<-function(Directory,NControl,NCondition,PSM,UP,CurveRsq,PValMelt,PValMeltFDR,MeltLimit,RunSTRING,STRINGScore,RunPANTHER,Species,PANTHERpvalue){
  Data_Imported<-Import(NControl,NCondition,Directory)
  message("Import Complete, Data Normalization Pending")
  Data_Normalized<-Normalize(Data_Imported)
  message("Data Normalization Complete, Proteome Quantification Pending")
  Data_Quantified<-Quantify(Data_Normalized)
  message("Proteome Quantification Complete, Proteome Curve Fitting Pending")
  Data_CurveFit1Parameters<-CurveFit1(Data_Quantified)
  message("Proteome Curve Fitting Complete, Protein Abundance Correction Pending")
  Data_Corrected<-Correction(PSM,UP,Data_CurveFit1Parameters,Data_Normalized,Data_Quantified)
  message("Protein Abundance Correction Complete, Protein Melt Curve Fitting Pending")
  Data_Corrected_Control<-subset(Data_Corrected,Data_Corrected$`Condition/Control`=="Control")
  Data_Corrected_Condition<-subset(Data_Corrected,Data_Corrected$`Condition/Control`=="Condition")
  Data_CurveFit2_Control<-CurveFit2(Data_Corrected_Control)
  Data_CurveFit2_Condition<-CurveFit2(Data_Corrected_Condition)
  Data_CurveFit2_Control_Unique<-Data_CurveFit2_Control[!duplicated(Data_CurveFit2_Control$Accession),]
  Data_CurveFit2_Condition_Unique<-Data_CurveFit2_Condition[!duplicated(Data_CurveFit2_Condition$Accession),]
  Data_CurveFit2_Complete_Unique<-rbind(Data_CurveFit2_Control_Unique,Data_CurveFit2_Condition_Unique)
  message("Protein Melt Curve Fitting Complete, Protein Melt Shift Calculation Pending")
  Data_Melts<-MeltCalc(Directory,Data_CurveFit2_Complete_Unique,CurveRsq,PValMelt,MeltLimit,PValMeltFDR)
  message("Protein Melt Shift Calculation Complete, Melt Shift Results Reporting Pending")
  ReportDataMelts(Data_Melts,Data_CurveFit2_Control,Data_CurveFit2_Condition,Directory)
  message("Melt Shift Results Reporting Complete")
  if(RunSTRING=="Yes"){
  message("STRING Reporting Pending")
  ReportSTRING(Data_Melts,STRINGScore,Directory,Species)
  message("STRING Reporting Complete")
  }
  if(RunPANTHER=="Yes"){
    message("PANTHER Reporting Pending")
    ReportPANTHER(Data_Melts,Directory,Species,PANTHERpvalue,PValMeltFDR)
    message("PANTHER Reporting Complete")
  }

}
