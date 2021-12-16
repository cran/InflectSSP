#' This function corrects the normalized abundance of each protein using a correction constant that is calculated in this function. The correction constant is determined using the difference between actual and predicted fit at the proteome level.
#' @param Data_CurveFit1Parameters the parameters determined from Curve Fit 1 operation for proteome melts
#' @param Data_Normalized the normalized abundance data for each protein determined in the Normalize function.
#' @param Data_Quantified the median normalized abundance data at the proteome level
#' @param PSM the number of peptide spectrum matches that are deemed acceptable for reporting
#' @param UP the number of unique peptides for a protein that are deemed acceptable for reporting
#' @importFrom stats complete.cases
#' @export
#' @return the corrected and normalized abundance data for each protein
#' @examples
#' \dontrun{
#' Data_Corrected<-Correction(PSM,UP,Data_CurveFit1Parameters,
#' Data_Normalized,Data_Quantified)
#' }

Correction<-function(PSM,UP,Data_CurveFit1Parameters,Data_Normalized,Data_Quantified){
  #Prepares data from previous step for analysis. Includes determining number of temperatures and number of proteins.
  NumTemps<-as.data.frame(Data_Normalized$NumberTemps)[1,]
  NumProteins<-as.numeric(nrow(Data_Normalized))
  Data_Normalized_Exclude_1<-Data_Normalized[(Data_Normalized[,2]>PSM),]
  Data_Normalized_Exclude<-Data_Normalized_Exclude_1[(Data_Normalized_Exclude_1[,3]>UP),]
  CorrectionValues_Constants<-data.frame(matrix(1,ncol=NumTemps))
  #Calculates correction constants for each protein and at each temperature
  B_Proteome<-Data_CurveFit1Parameters[1,3]
  T_Proteome<-Data_CurveFit1Parameters[1,4]
  xmid_Proteome<-Data_CurveFit1Parameters[1,5]
  b_Proteome<-Data_CurveFit1Parameters[1,6]
  s_Proteome<-1
  Temps_Correction<-data.frame(t(Data_Normalized_Exclude[1,(4+NumTemps):(3+NumTemps*2)]))
  colnames(Temps_Correction)<-"Temperature"
  CorrectionValues<-B_Proteome+((T_Proteome-B_Proteome)/(1 + 10^(b_Proteome*(xmid_Proteome-log10(Temps_Correction))))^s_Proteome)
  Data_Quantified_Rep<-merge(Temps_Correction,Data_Quantified, by="Temperature")
  CorrectionValues_Constants<-t(as.data.frame((1-((as.numeric(Data_Quantified_Rep[,2])-CorrectionValues)/(as.numeric(Data_Quantified_Rep[,2]))))))
  #Creates data frames for later analysis steps
  Data_Normalized_Exclude_Corrected<-data.frame(matrix(nrow=as.numeric(nrow(Data_Normalized_Exclude)),ncol=NumTemps))
  Data_Normalized_Exclude_Corrected_Parameters<-data.frame(matrix(nrow=as.numeric(nrow(Data_Normalized_Exclude)),ncol=4))
  count<-1
  #Corrects protein abundance for each protein in the experiment based on the correction constants calculated above.
  repeat{
    Data_Normalized_Exclude_Corrected[count,1:NumTemps]<-Data_Normalized_Exclude[count,(7+2*NumTemps):(6+3*NumTemps)]*CorrectionValues_Constants[1,]
    Data_Normalized_Exclude_Corrected_Parameters[count,1:4]<-Data_CurveFit1Parameters[1,3:6]
    count<-count+1
    if(count>nrow(Data_Normalized_Exclude)){
      break
    }
  }
  colnames(Data_Normalized_Exclude_Corrected)<-paste(colnames(Data_Normalized_Exclude[,(7+2*NumTemps):(6+3*NumTemps)]), "Corrected",sep = "_")
  colnames(Data_Normalized_Exclude_Corrected_Parameters)<-colnames(Data_CurveFit1Parameters[,3:6])
  Data_Corrected<-cbind(Data_Normalized_Exclude,cbind(Data_Normalized_Exclude_Corrected,Data_Normalized_Exclude_Corrected_Parameters))
rm(Data_Normalized_Exclude,Data_Normalized_Exclude_1)
return(Data_Corrected)
}
