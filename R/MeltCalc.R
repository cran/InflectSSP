#' This function determines melt shifts for all proteins that meet quality criteria and also determines the melt shift p-values
#' @param Data_CurveFit2_Complete_Unique the curve fit data from the CurveFit2 function
#' @param CurveRsq the criteria for melt curve p-values
#' @param PValMelt the criteria for the melt shift p-values
#' @param Directory the directory data is saved to
#' @param MeltLimit the melt shift temperature limit used for determining which proteins are significant
#' @param PValMeltFDR Whether or not the FDR correction for pvalue is used in designation of melts of interest
#' @importFrom utils write.csv
#' @importFrom stats pnorm
#' @importFrom stats p.adjust
#' @export
#' @return Proteins melt shifts
#' @examples
#' \dontrun{
#'      Data_Melts<-MeltCalc(Directory,Data_CurveFit2_Complete_Unique,
#'      CurveRsq,PValMelt,MeltLimit,PValMeltFDR)}

MeltCalc<-function(Directory,Data_CurveFit2_Complete_Unique,CurveRsq,PValMelt,MeltLimit,PValMeltFDR){
#Determines number of temperatures and proteins for further analysis
  NumTemps<-as.data.frame(Data_CurveFit2_Complete_Unique$NumberTemps)[1,]
  Protein<-1
  ControlData<-subset(Data_CurveFit2_Complete_Unique,Data_CurveFit2_Complete_Unique$`Condition/Control`=="Control")
  ConditionData<-subset(Data_CurveFit2_Complete_Unique,Data_CurveFit2_Complete_Unique$`Condition/Control`=="Condition")
  colnames(ControlData)<-paste(colnames(ControlData),"Control",sep="-")
  colnames(ControlData)[1]<-"Accession"
  colnames(ConditionData)<-paste(colnames(ConditionData),"Condition",sep="-")
  colnames(ConditionData)[1]<-"Accession"
  Data_CurveFit2_Merged<-merge(ControlData,ConditionData,by="Accession")
  Melt<-data.frame(matrix(ncol=8,nrow=nrow(Data_CurveFit2_Merged)))
  colnames(Melt)<-c("MeltControl","MeltCondition","MeltShift","Melt_pValue","Melt_pValue_FDRAdj","Good_Curve_Rsq","Good_Melt_pValue","Good_Melt_Magnitude")

  repeat{
    Melt_Control<-(Data_CurveFit2_Merged$`xmid_Const -Protein-Control`[Protein])
    Melt_Condition<-(Data_CurveFit2_Merged$`xmid_Const -Protein-Condition`[Protein])
    Melt_SE_Control<-(Data_CurveFit2_Merged$`Xmid_SE-Control`[Protein])
    Melt_SE_Condition<-(Data_CurveFit2_Merged$`Xmid_SE-Condition`[Protein])
    MeltpValue<-(1-pnorm(abs((Melt_Control-Melt_Condition)/(sqrt(Melt_SE_Control^2 + Melt_SE_Condition^2)))))*2
    Melt[Protein,1]<-10^(Melt_Control)
    Melt[Protein,2]<-10^(Melt_Condition)
    Melt[Protein,3]<-as.data.frame(10^Melt_Condition-10^Melt_Control)
    Melt[Protein,4]<-MeltpValue
    ContRsq<-Data_CurveFit2_Merged$`Rsq-Control`[Protein]
    CondRsq<-Data_CurveFit2_Merged$`Rsq-Condition`[Protein]
  #Determines whether a melt shift meets criteria from user
    if(abs(Melt[Protein,3])>MeltLimit){
      GoodMeltLimit<-"Yes"
    } else {
      GoodMeltLimit<-"No"
    }
  #Determines whether a curve fit meets criteria from user
    if(ContRsq>CurveRsq & CondRsq>CurveRsq){
      GoodCurveRsq<-"Yes"
    } else {
      GoodCurveRsq<-"No"
    }
  #Determines whether a melt shift p-value meets user criteria
    if(MeltpValue<PValMelt){
      GoodMeltPVal<-"Yes"
    } else {
      GoodMeltPVal<-"No"
    }

    Melt[Protein,6]<-GoodCurveRsq

    Melt[Protein,8]<-GoodMeltLimit
    Protein<-Protein+1

    if(Protein>nrow(Data_CurveFit2_Merged)){
      break
    }

  }

colnames(Melt)<-c("Melt Control","Melt Condition","Melt Shift","Melt p-value","Melt_pValue_FDRAdj","Good Rsq","Good Melt p-value","Good Melt Limit")

#calculates FDR corrected p-values
MeltpValue<-data.frame(p.adjust(Melt$`Melt p-value`,method = "fdr"))
Melt[,5]<-MeltpValue

#Determines whether melts satisfy p-value criteria for melts
Protein<-1
repeat{

  if(PValMeltFDR=="No"){
  if(Melt[Protein,4]<PValMelt){
    GoodMeltPVal<-"Yes"
  } else {
    GoodMeltPVal<-"No"
  }
  } else {
    if(Melt[Protein,5]<PValMelt){
      GoodMeltPVal<-"Yes"
    } else {
      GoodMeltPVal<-"No"
    }
  }
  Melt[Protein,7]<-GoodMeltPVal
  Protein<-Protein+1
  if(Protein>nrow(Data_CurveFit2_Merged)){
    break
  }
}



Data_Melts_preordered<-cbind(Data_CurveFit2_Merged,Melt)
MeltCoefficient<-log10(Data_Melts_preordered$`PSM-Control`)*log10(Data_Melts_preordered$`Unique Peptides-Control`)*(Data_Melts_preordered$`Rsq-Control`)*(Data_Melts_preordered$`Rsq-Condition`)*(-log10(Data_Melts_preordered$`Melt p-value`))*abs(Data_Melts_preordered$`Melt Shift`)
Data_Melts_preordered_2<-cbind(Data_Melts_preordered,MeltCoefficient)
Data_Melts<-as.data.frame(Data_Melts_preordered_2[order(Data_Melts_preordered_2$`Melt Shift`),])



  return(Data_Melts)

  }


