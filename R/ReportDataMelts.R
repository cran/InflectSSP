#' This function generates results from the Inflect function after applying criteria input from the user
#' @param Data_Melts abundance and fit data for proteins that meet quality criteria in overall workflow
#' @param Data_CurveFit2_Control the curve fit data from the Curve Fit 2 function
#' @param Data_CurveFit2_Condition the curve fit data from the Curve Fit 2 function
#' @param Directory directory where data is saved
#' @importFrom xlsx write.xlsx
#' @importFrom grDevices dev.off
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics axis
#' @importFrom graphics legend
#' @importFrom plotrix addtable2plot
#' @importFrom data.table data.table
#' @importFrom grDevices pdf
#' @importFrom graphics grid
#' @importFrom graphics par
#' @importFrom utils write.csv
#' @export
#' @return files with summary of data along with melt curve plots for significant proteins
#' @examples
#'      \dontrun{
#'      ReportDataMelts(Data_Melts,Data_CurveFit2_Control,Data_CurveFit2_Condition,Directory)}

ReportDataMelts<-function(Data_Melts,Data_CurveFit2_Control,Data_CurveFit2_Condition,Directory){

#Subsets proteins based on whether they meet Rsq, pValue and melt limit criteria
  row.names(Data_Melts) <- 1:nrow(Data_Melts)
  Data_Melts_GoodRsq<-Data_Melts[Data_Melts$`Good Rsq`=="Yes",]
  Data_Melts_GoodRsq_GoodpVal<-Data_Melts_GoodRsq[Data_Melts_GoodRsq$`Good Melt p-value`=="Yes",]
  Data_Melts_GoodRsq_GoodpVal_GoodMelt<-subset(Data_Melts_GoodRsq_GoodpVal,Data_Melts_GoodRsq_GoodpVal$`Good Melt Limit`=="Yes")
  Data_CurveFit2_Both_Pre<-rbind(Data_CurveFit2_Control,Data_CurveFit2_Condition)
  Data_CurveFit2_Both<-Data_CurveFit2_Both_Pre[order(Data_CurveFit2_Both_Pre$Accession),]
  row.names(Data_CurveFit2_Both) <- 1:nrow(Data_CurveFit2_Both)
  dir.create(paste(Directory,"Result Files",sep="/"))
  dir.create(paste(Directory,"Result Files/Melt Curves",sep="/"))
  if(NROW(Data_Melts)>0){
#Determines number of temperatures and proteins for further analysis
  NumProteins<-as.numeric(nrow(Data_CurveFit2_Both))
  NumTemps<-Data_CurveFit2_Both$NumberTemps[1]
#This section creates melt curves for all proteins that meet curve fit criteria
  Protein<-1
  repeat{
    Abundance_Control<-data.frame(matrix(nrow=1,ncol=NumTemps))
    Abundance_Condition<-data.frame(matrix(nrow=1,ncol=NumTemps))
    Temperature_Control<-data.frame(matrix(nrow=1,ncol=NumTemps))
    Temperature_Condition<-data.frame(matrix(nrow=1,ncol=NumTemps))
    ProteinID<-Data_CurveFit2_Both[Protein,1]
    ControlCount<-0
    ConditionCount<-0
    repeat{

      if(ProteinID==Data_CurveFit2_Both[Protein,1]){

        if(Data_CurveFit2_Both$`Condition/Control`[Protein]=="Control"){
        Abundance_Control<-cbind(Abundance_Control,as.data.frame(Data_CurveFit2_Both[Protein,(7+NumTemps*3):(6+NumTemps*4)]))
        Temperature_Control<-cbind(Temperature_Control,as.data.frame(Data_CurveFit2_Both[Protein,(4+NumTemps):(3+NumTemps*2)]))
        T_Control<-Data_CurveFit2_Both[Protein,(11+NumTemps*4)]
        B_Control<-Data_CurveFit2_Both[Protein,(12+NumTemps*4)]
        xmid_Control<-Data_CurveFit2_Both[Protein,(13+NumTemps*4)]
        scal_Control<-Data_CurveFit2_Both[Protein,(14+NumTemps*4)]
        Temperature_Control<-Temperature_Control[,is.na(Temperature_Control)==FALSE]
        FitTemps<-as.numeric(c((min(Temperature_Control)-5):(max(Temperature_Control)+5)))
        CurveFit_Control<-B_Control+((T_Control-B_Control)/(1 + 10^(scal_Control*(xmid_Control-log10(FitTemps))))^1)
        CurveFit_ControlData<-data.frame(as.numeric(FitTemps),as.numeric(CurveFit_Control))
        colnames(CurveFit_ControlData)<-c("x","y")
        ProteinStart<-Protein
        Protein<-Protein+1
        ControlCount<-ControlCount+1
        } else {
          Abundance_Condition<-cbind(Abundance_Condition,as.data.frame(Data_CurveFit2_Both[Protein,(7+NumTemps*3):(6+NumTemps*4)]))
          Temperature_Condition<-cbind(Temperature_Condition,as.data.frame(Data_CurveFit2_Both[Protein,(4+NumTemps):(3+NumTemps*2)]))
          T_Condition<-Data_CurveFit2_Both[Protein,(11+NumTemps*4)]
          B_Condition<-Data_CurveFit2_Both[Protein,(12+NumTemps*4)]
          xmid_Condition<-Data_CurveFit2_Both[Protein,(13+NumTemps*4)]
          scal_Condition<-Data_CurveFit2_Both[Protein,(14+NumTemps*4)]
          Temperature_Condition<-Temperature_Condition[,is.na(Temperature_Condition)==FALSE]
          FitTemps<-as.numeric(c((min(Temperature_Condition)-5):(max(Temperature_Condition)+5)))
          CurveFit_Condition<-B_Condition+((T_Condition-B_Condition)/(1 + 10^(scal_Condition*(xmid_Condition-log10(FitTemps))))^1)
          CurveFit_ConditionData<-data.frame(as.numeric(FitTemps),as.numeric(CurveFit_Condition))
          colnames(CurveFit_ConditionData)<-c("x","y")
          ProteinStart<-Protein
          Protein<-Protein+1
          ConditionCount<-ConditionCount+1
        }

        if (Protein>NumProteins){
          break
        }

      } else {
        break
      }
    }

    if(ControlCount>0 & ConditionCount>0){
      Abundance_Control<-Abundance_Control[,is.na(Abundance_Control)==FALSE]
      Abundance_Condition<-Abundance_Condition[,is.na(Abundance_Condition)==FALSE]
      Temperature_Control<-Temperature_Control[,is.na(Temperature_Control)==FALSE]
      Temperature_Condition<-Temperature_Condition[,is.na(Temperature_Condition)==FALSE]
      Abundance_Control_Plot<-cbind(as.numeric(t(Temperature_Control)),t(Abundance_Control))
      row.names(Abundance_Control_Plot)<-1:nrow(Abundance_Control_Plot)
      colnames(Abundance_Control_Plot)<-c("Temperature","Abundance")
      CurveFit_ControlData_Plot<-cbind(CurveFit_ControlData,data.frame(matrix(data='Control-Fit',nrow=nrow(CurveFit_ControlData),ncol=1)))
      colnames(CurveFit_ControlData_Plot)<-c("Temperature","Abundance")
      Abundance_Condition_Plot<-cbind(as.numeric(t(Temperature_Condition)),t(Abundance_Condition))
      row.names(Abundance_Condition_Plot)<-1:nrow(Abundance_Condition_Plot)
      colnames(Abundance_Condition_Plot)<-c("Temperature","Abundance")
      CurveFit_ConditionData_Plot<-cbind(CurveFit_ConditionData,data.frame(matrix(data='Control-Fit',nrow=nrow(CurveFit_ConditionData),ncol=1)))
      colnames(CurveFit_ConditionData_Plot)<-c("Temperature","Abundance")

      pdf(paste(Directory,paste("Result Files/Melt Curves",paste(ProteinID,"pdf",sep="."),sep="/"),sep="/"))
      plot(x = Abundance_Control_Plot[,1] , y = Abundance_Control_Plot[,2] ,pch = 2, frame = TRUE,xlab = "Temperature (C)", ylab = "Normalized Abundance",col = "blue",axes=FALSE,ylim=c(0,min(1.5,1.2*max(Data_Melts$`T1-Abundance_Normalized-Control`))),xlim=c(20,95),main=ProteinID)
      lines(CurveFit_ControlData_Plot[,1],CurveFit_ControlData_Plot[,2],col="blue")
      points(Abundance_Condition_Plot[,1], Abundance_Condition_Plot[,2], col="red", pch=8,lty=1)
      lines(CurveFit_ConditionData_Plot[,1],CurveFit_ConditionData_Plot[,2],col="red")
      axis(side=1, at=seq(20,100,by=5))
      axis(side=2, at=seq(0,1.5, by=.1))
      legend(20,0.2,legend=c("Control","Condition"), col=c("blue","red"),pch=c(2,8),lty=c(1,2))
      plottable=matrix(data=NA,nrow=7,ncol=2)
      colnames(plottable)<-c("Parameter","Value")
      plottable[1,1]<-"Control Tm"
      plottable[2,1]<-"Condition Tm"
      plottable[3,1]<-"Delta Tm"
      plottable[4,1]<-"pr(>|z|)"
      plottable[5,1]<-"Adj pr(>|z|)"
      plottable[6,1]<-"Control Rsq"
      plottable[7,1]<-"Condition Rsq"
      Data_Melts_RowNumber<-as.numeric(which(Data_Melts$Accession==ProteinID))
      plottable[1,2]<-round(as.numeric(Data_Melts[Data_Melts_RowNumber,c(NumTemps*8+44)]),2)
      plottable[2,2]<-round(as.numeric(Data_Melts[Data_Melts_RowNumber,c(NumTemps*8+45)]),2)
      plottable[3,2]<-round(as.numeric(Data_Melts[Data_Melts_RowNumber,c(NumTemps*8+46)]),2)
      plottable[4,2]<-formatC(Data_Melts[Data_Melts_RowNumber,c(NumTemps*8+47)], format = "e", digits = 2)
      plottable[5,2]<-formatC(Data_Melts[Data_Melts_RowNumber,c(NumTemps*8+48)], format = "e", digits = 2)
      plottable[6,2]<-round(as.numeric(Data_Melts[Data_Melts_RowNumber,c(NumTemps*4+22)]),2)
      plottable[7,2]<-round(as.numeric(Data_Melts[Data_Melts_RowNumber,c(NumTemps*8+43)]),2)
      addtable2plot(60,0.8,plottable,hlines=TRUE,vlines=TRUE,bty="o",bg="gray")
      dev.off()
    }


  if (Protein>NumProteins){
      break
    }
  }





#Rank order plot that contains all proteins that have met R squared criteria. Significant proteins are in red and are those that meet p-value and melt shift criteria
  pdf(paste(Directory,paste("Result Files","WaterfallPlot.pdf",sep="/"),sep="/"))
  plot(x = 1:nrow(Data_Melts), y = Data_Melts$`Melt Shift`,pch = 2, frame = TRUE,xlab = "Protein #", ylab = "Melt Shift, Condition-Control (C)", cex=0.5)
  points(row.names(Data_Melts_GoodRsq_GoodpVal_GoodMelt),as.numeric(Data_Melts_GoodRsq_GoodpVal_GoodMelt$`Melt Shift`),col="red",pch=8,lty=1)
  legend(0,0.9*max(Data_Melts$`Melt Shift`),legend=c("Melt","Melt of Interest"), col=c("black","red"),pch=c(2,8),lty=c(1,2))
  grid(col = "lightgray", lty = "dotted",lwd = par("lwd"))
  dev.off()

#Writes results to files in the directory chosen by the user.

  Data_Melts_worawdata<-cbind(Data_Melts$Accession,cbind(Data_Melts[,(4*NumTemps+7):(4*NumTemps+22)],Data_Melts[,(8*NumTemps+28):(ncol(Data_Melts))]))
  Data_Melts_GoodRsq_GoodpVal_GoodMelt_worawdata<-cbind(Data_Melts_GoodRsq_GoodpVal_GoodMelt$Accession,cbind(Data_Melts_GoodRsq_GoodpVal_GoodMelt[,(4*NumTemps+7):(4*NumTemps+22)],Data_Melts_GoodRsq_GoodpVal_GoodMelt[,(8*NumTemps+28):(ncol(Data_Melts_GoodRsq_GoodpVal_GoodMelt))]))

  colnames(Data_Melts_worawdata)[1]<-"Accession"
  colnames(Data_Melts_GoodRsq_GoodpVal_GoodMelt_worawdata)[1]<-"Accession"

  write.csv(Data_CurveFit2_Both, file = paste(Directory,paste("Result Files","AnalysisResults.csv",sep="/"),sep="/"))
  write.csv(Data_Melts_worawdata, file = paste(Directory,paste("Result Files","AllMeltShifts.csv",sep="/"),sep="/"))
  write.csv(Data_Melts_GoodRsq_GoodpVal_GoodMelt_worawdata, file = paste(Directory,paste("Result Files","FilteredMeltShifts.csv",sep="/"),sep="/"))



  } else {
    print("No melt shifts with good fit")
  }

}
