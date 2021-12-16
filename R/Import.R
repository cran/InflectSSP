#' This function imports data that will be analyzed in downstream functions.
#' @param Directory the directory where the source data files to be analyzed are saved. This is also the location where the results will be saved.
#' @param NControl the number of Control replicate experiments that are to be analyzed
#' @param NCondition the number of Condition replicate experiments that are to be analyzed
#' @importFrom readxl read_excel
#' @export
#' @return Imported data from all experiments
#' @examples
#' \dontrun{
#' Data_Imported<-Import(NControl,NCondition,Directory)
#' }

Import<-function(NControl,NCondition,Directory){
#Determines the number of replicate experiments based on the maximum number of control or condition experiments
NReps<-max(NControl,NCondition)

Rep<-1

  repeat{
#This section of code imports the data from the Control data files
  if(Rep<=NControl){
  ControlData<-readxl::read_excel(paste(Directory,paste("/Control ",paste(Rep,".xlsx",sep=""),sep=""),sep=""))
  NumTemps<-ncol(ControlData)-3
  ColumnNames_Control<-paste("T",seq(1,NumTemps,by=1),sep="")
  Temps_Import<-data.frame(matrix(nrow=1,ncol=NumTemps))
  Temps_Import<-data.frame(as.numeric(colnames(ControlData)[4:(3+NumTemps)]))
  ControlData[,(4+NumTemps):(3+2*NumTemps)]<-t(Temps_Import)
  ControlData[,(4+2*NumTemps)]<-NumTemps
  ControlData[,(5+2*NumTemps)]<-"Control"
  colnames(ControlData)<-c("Accession","PSM","Unique Peptides",paste(ColumnNames_Control,"Abundance",sep="-"),paste(ColumnNames_Control,"Temperature",sep="-"),"NumberTemps","Condition/Control")
  }

#This section of code imports the data from the Condition data files
  if(Rep<=NCondition){
  ConditionData<-readxl::read_excel(paste(Directory,paste("/Condition ",paste(Rep,".xlsx",sep=""),sep=""),sep=""))
  NumTemps<-ncol(ConditionData)-3
  ColumnNames_Condition<-paste("T",seq(1,NumTemps,by=1),sep="")
  Temps_Import<-data.frame(matrix(nrow=1,ncol=NumTemps))
  Temps_Import<-data.frame(as.numeric(colnames(ConditionData)[4:(3+NumTemps)]))
  ConditionData[,(4+NumTemps):(3+2*NumTemps)]<-t(Temps_Import)
  ConditionData[,(4+2*NumTemps)]<-NumTemps
  ConditionData[,(5+2*NumTemps)]<-"Condition"
  colnames(ConditionData)<-c("Accession","PSM","Unique Peptides",paste(ColumnNames_Condition,"Abundance",sep="-"),paste(ColumnNames_Condition,"Temperature",sep="-"),"NumberTemps","Condition/Control")
  }

#This section of code merges the Control and Condition data into a single spreadsheet with a column that designates which replicate the data is from.
  NumRows_Control<-as.numeric(nrow(ControlData))
  NumRows_Condition<-as.numeric(nrow(ConditionData))
  Replicate_Control<-data.frame(matrix(data=Rep,nrow=NumRows_Control,ncol=1))
  Replicate_Condition<-data.frame(matrix(data=Rep,nrow=NumRows_Condition,ncol=1))
  colnames(Replicate_Control)<-"Replicate"
  colnames(Replicate_Condition)<-"Replicate"
  ControlData_WithRep<-cbind(ControlData,Replicate_Control)
  ConditionData_WithRep<-cbind(ConditionData,Replicate_Condition)
  if(Rep==1){
    ControlData_AllReps<-ControlData_WithRep
  } else {
    ControlData_AllReps<-rbind(ControlData_AllReps,ControlData_WithRep)
  }

  if(Rep==1){
    ConditionData_AllReps<-ConditionData_WithRep
  } else {
    ConditionData_AllReps<-rbind(ConditionData_AllReps,ConditionData_WithRep)
  }

  Rep<-Rep+1
  ControlData[,]<-NA
  ConditionData[,]<-NA

  if(Rep>NReps){
    break
  }

  }

#This step removes rows in data frame that do not have complete data for example where abundance data is missing and then merges
#the control and condition data sets into a single table.
ControlData_AllReps_Complete<-ControlData_AllReps[complete.cases(ControlData_AllReps),]
ConditionData_AllReps_Complete<-ConditionData_AllReps[complete.cases(ConditionData_AllReps),]
Data_Imported<-rbind(ControlData_AllReps_Complete,ConditionData_AllReps_Complete)

#Removes variables from memory that will not be used later in the program.
rm(ControlData_AllReps_Complete,ConditionData_AllReps_Complete,ControlData_AllReps,ConditionData_AllReps,ControlData,ConditionData,ConditionData_WithRep,ControlData_WithRep,Replicate_Condition,Replicate_Control)
return(Data_Imported)
}

