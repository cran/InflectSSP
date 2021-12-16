#' This function normalizes the abundance values to that measured at the lowest temperature
#' @param Data_Imported the abundance data imported from Import function
#' @export
#' @return Normalized data
#' @examples
#' \dontrun{
#'      Data_Normalized<-Normalize(Data_Imported)}

Normalize<-function(Data_Imported){
  #Determines the number of temperatures in the experiment.
  NumTemps<-as.data.frame(Data_Imported$NumberTemps)[1,]
  #Omits rows that don't have complete data reported and thus can't be analyzed.
  Data_Imported_ZeroOmit<-subset(Data_Imported,Data_Imported$`T1-Abundance`>0)
  #Ensures data is in numeric format so that it can be analyzed
  Data_Imported_ZeroOmit[4:(3+NumTemps)]<-sapply(Data_Imported_ZeroOmit[4:(3+NumTemps)],as.numeric)
  Data_Imported_ZeroOmit[(4+NumTemps):(3+2*NumTemps)]<-sapply(Data_Imported_ZeroOmit[(4+NumTemps):(3+2*NumTemps)],as.numeric)
  #Data normalization calculation
  Data_Imported_Normalized<-Data_Imported_ZeroOmit[,4:(3+NumTemps)]/Data_Imported_ZeroOmit[,4]
  colnames(Data_Imported_Normalized)<-paste(colnames(Data_Imported_Normalized), "Normalized",sep = "_")
  Data_Normalized<-cbind(Data_Imported_ZeroOmit,Data_Imported_Normalized)
  #Removes variables from memory that will not be used later in the program.
  rm(Data_Imported,Data_Imported_Normalized,Data_Imported_ZeroOmit)
  return(Data_Normalized)
}
