#' This function determines the median abundance value across the proteome for all experiments together
#' @param Data_Normalized the normalized abundance data calculated in the Normalize function
#' @param NReps the number of replicates to be analyzed
#' @export
#' @return The median abundance data for all experiments at the proteome level
#' @examples
#' \dontrun{
#'      Data_Quantified<-Quantify(Data_Normalized)}

Quantify<-function(Data_Normalized,NReps){
  NumTemps<-as.data.frame(Data_Normalized$NumberTemps)[1,]
  ProteinCounter<-1
  RepCount<-1
#Creates columns in the table that correspond with the temperature treatments in the experiment.
  TemperatureSummary<-data.frame(Data_Normalized[1,(6+NumTemps*2):(5+NumTemps*3)])
  colnames(TemperatureSummary)<-"Temperature"
  Data_Normalized_Quant<-as.data.frame(t(unique(rbind(TemperatureSummary[order(TemperatureSummary$Temperature),],TemperatureSummary[order(TemperatureSummary$Temperature),]))))
  colnames(Data_Normalized_Quant)<-"Temperature"
  rownames(Data_Normalized_Quant)<-1:NROW(Data_Normalized_Quant)
#calculation of median abundance for each temperature in the heat treatment.
  Temperatures<-t(Data_Normalized[1,(4+NumTemps):(3+2*NumTemps)])
  Data_Normalized_Abundance<-Data_Normalized[,(7+NumTemps*2):(6+NumTemps*3)]
  Data_Quantified<-cbind(Temperatures,data.frame(apply(Data_Normalized_Abundance,2,median)))
  colnames(Data_Quantified)<-c("Temperature","Median_Abundance")
#Removes variables from memory that will not be used later in the program.
  rm(Data_Normalized_Abundance)
  return(Data_Quantified)
}
