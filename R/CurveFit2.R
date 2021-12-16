#' This function determines the best curve fit for each protein using the data post correction and also determines the R squared for each curve fit
#' @param Data_Corrected data that meets exclusion criteria from Exclude function
#' @importFrom xlsx write.xlsx
#' @importFrom stats complete.cases
#' @importFrom stats nls
#' @importFrom stats nls.control
#' @export
#' @return Curve fits and R squared for each protein
#' @examples
#' \dontrun{
#' Data_CurveFit2_Control<-CurveFit2(Data_Corrected_Control)}

CurveFit2<-function(Data_Corrected){
  #Determines number of temperatures and proteins for further analysis
  NumTemps<-as.data.frame(Data_Corrected$NumberTemps)[1,]
  NumProteins<-as.numeric(nrow(Data_Corrected))
  Protein<-1
  CurveFitStats<-data.frame(matrix(ncol = 11, nrow = NumProteins))
  Data_Corrected_Ordered<-Data_Corrected[order(Data_Corrected$Accession),]
  NumberRows<-as.numeric(nrow(Data_Corrected_Ordered))
  row.names(Data_Corrected_Ordered) <- 1:NumberRows
  #Creates a table of abundance values for both condition and control experiments in one repeat and then fits the data to a four parameter fit in the next set of repeat code.
  repeat{
    ProteinCount<-1
    repeat{
      if(ProteinCount==1){
        ProteinReplicateNumber<-1
        Temperature<-as.data.frame(Data_Corrected_Ordered[Protein,(4+NumTemps):(3+NumTemps*2)])
        Abundance<-as.data.frame(Data_Corrected_Ordered[Protein,(7+NumTemps*3):(6+NumTemps*4)])
        ProteinCount<-ProteinCount+1
        ProteinStart<-Protein
        ProteinID<-Data_Corrected_Ordered[Protein,1]
        Protein<-Protein+1
      } else {
        if(ProteinCount>1 & ProteinID==Data_Corrected_Ordered[Protein,1]){
          Temperature<-cbind(Temperature,as.data.frame(Data_Corrected_Ordered[Protein,(4+NumTemps):(3+NumTemps*2)]))
          Abundance<-cbind(Abundance,as.data.frame(Data_Corrected_Ordered[Protein,(7+NumTemps*3):(6+NumTemps*4)]))
          ProteinCount<-ProteinCount+1
          Protein<-Protein+1
        } else {
          ProteinEnd<-Protein-1
          break
        }
      }
      if (Protein>NumProteins){
        ProteinEnd<-(Protein-1)
        break
      }
    }

    Abundance<-as.numeric(t(Abundance))
    Temperature<-as.numeric(t(Temperature))
    AbundTemp<-as.data.frame(cbind(Temperature,Abundance))
    CurveFit2<-nls(Abundance ~ B_Const+((T_Const-B_Const)/(1 + 10^(scal_Const*(xmid_Const-log10(Temperature)))))^1, data = AbundTemp, start = c(T_Const=1,B_Const=.1,xmid_Const=1.7,scal_Const=-10),nls.control(maxiter = 500,warnOnly = TRUE))
    T_Const<-CurveFit2$m$getAllPars()[1]
    B_Const<-CurveFit2$m$getAllPars()[2]
    xmid_Const<-CurveFit2$m$getAllPars()[3]
    scal_Const<-CurveFit2$m$getAllPars()[4]


  #If the protein is able to be fit to a 4PL, the number of parameters in the fit is listed as 4. In the event that the nls function does not converge, a 3PL fit is attempted.
  #If no fit is possible, the protein is skipped and NA is listed for fit parameters. If the melt temperature is less than or greater than the minimum or maximum temperatures in
  #the temperature gradient, a 3PL is attempted.
    if(CurveFit2$convInfo$isConv==TRUE){
      CurveInflection<-10^(CurveFit2$m$getAllPars()[3])
      FitParameters<-4
      CurveFit2_Summary<-summary(CurveFit2)
      CurveSE<-CurveFit2_Summary$sigma
      if(CurveInflection < min(Temperature) || CurveInflection > max(Temperature)) {
        CurveFit2<-nls(Abundance ~ 0+((T_Const-0)/(1 + 10^(scal_Const*(xmid_Const-log10(Temperature)))))^1, data = AbundTemp, start = c(T_Const=1,xmid_Const=1.7,scal_Const=-10),nls.control(maxiter = 500,warnOnly = TRUE))
        if(CurveFit2$convInfo$isConv==TRUE ){
          T_Const<-CurveFit2$m$getAllPars()[1]
          B_Const<-0
          xmid_Const<-CurveFit2$m$getAllPars()[2]
          scal_Const<-CurveFit2$m$getAllPars()[3]
          CurveFit2_Summary<-summary(CurveFit2)
          CurveSE<-CurveFit2_Summary$sigma
          PVal<-data.frame(t(CurveFit2_Summary$parameters[,4]))
          Xmid_SE<-data.frame(CurveFit2_Summary$parameters[3,2])
          Data_CurveFit2Pars<-cbind(T_Const,cbind(B_Const,cbind(xmid_Const,scal_Const)))
          colnames(Data_CurveFit2Pars)<-paste(colnames(Data_CurveFit2Pars),"-Protein")
          colnames(PVal)<-paste("P-value",colnames(PVal))
          colnames(Xmid_SE)<-"Xmid_SE"
          FitParameters<-3
        } else {
          FitParameters<-NA
          T_Const<-NA
          B_Const<-NA
          xmid_Const<-NA
          scal_Const<-NA
          CurveSE<-NA
          PVal<-data.frame(t(c(NA,NA,NA,NA)))
          colnames(PVal)<-c("T_Const","B_Const","xmid_Const","scal_Const")
        }
      } else {

        CurveFit2_Summary<-summary(CurveFit2)
        CurveSE<-CurveFit2_Summary$sigma
        PVal<-data.frame(t(CurveFit2_Summary$parameters[,4]))
        Xmid_SE<-data.frame(CurveFit2_Summary$parameters[3,2])
        Data_CurveFit2Pars<-cbind(T_Const,cbind(B_Const,cbind(xmid_Const,scal_Const)))
        colnames(Data_CurveFit2Pars)<-paste(colnames(Data_CurveFit2Pars),"-Protein")
        colnames(PVal)<-paste("P-value",colnames(PVal))
        colnames(Xmid_SE)<-"Xmid_SE"
      }

    } else {


      CurveFit2<-nls(Abundance ~ 0+((T_Const-0)/(1 + 10^(scal_Const*(xmid_Const-log10(Temperature)))))^1, data = AbundTemp, start = c(T_Const=1,xmid_Const=1.7,scal_Const=-10),nls.control(maxiter = 500,warnOnly = TRUE))
      if(CurveFit2$convInfo$isConv==TRUE ){
        T_Const<-CurveFit2$m$getAllPars()[1]
        B_Const<-0
        xmid_Const<-CurveFit2$m$getAllPars()[2]
        scal_Const<-CurveFit2$m$getAllPars()[3]
        CurveFit2_Summary<-summary(CurveFit2)
        CurveSE<-CurveFit2_Summary$sigma
        PVal<-data.frame(t(CurveFit2_Summary$parameters[,4]))
        Xmid_SE<-data.frame(CurveFit2_Summary$parameters[3,2])
        Data_CurveFit2Pars<-cbind(T_Const,cbind(B_Const,cbind(xmid_Const,scal_Const)))
        colnames(Data_CurveFit2Pars)<-paste(colnames(Data_CurveFit2Pars),"-Protein")
        colnames(PVal)<-paste("P-value",colnames(PVal))
        colnames(Xmid_SE)<-"Xmid_SE"
        FitParameters<-3
      } else {
        FitParameters<-NA
        T_Const<-NA
        B_Const<-NA
        xmid_Const<-NA
        scal_Const<-NA
        CurveSE<-NA
        PVal<-data.frame(t(c(NA,NA,NA,NA)))
        colnames(PVal)<-c("T_Const","B_Const","xmid_Const","scal_Const")
        Xmid_SE<-data.frame(c(NA))
        Data_CurveFit2Pars<-cbind(T_Const,cbind(B_Const,cbind(xmid_Const,scal_Const)))
        colnames(Data_CurveFit2Pars)<-paste(colnames(Data_CurveFit2Pars),"-Protein")
        colnames(PVal)<-paste("P-value",colnames(PVal))
        colnames(Xmid_SE)<-"Xmid_SE"
      }


    }

    Data_CurveFit2_Summary<-cbind(Data_CurveFit2Pars,cbind(FitParameters,cbind(CurveSE,cbind(PVal,Xmid_SE))))
    CurveFitStats[ProteinStart:ProteinEnd,]<-Data_CurveFit2_Summary

    if (Protein>NumProteins){
      break
    }

  }

  colnames(CurveFitStats)<-colnames(Data_CurveFit2_Summary)
  Data_CurveFit2<-cbind(Data_Corrected_Ordered,CurveFitStats)
  Data_CurveFit2_woNA<-Data_CurveFit2[complete.cases(Data_CurveFit2),]


  #Determines the R squared for each protein curve fit.

  Rsq<-data.frame(matrix(ncol = 1, nrow = NROW(Data_CurveFit2_woNA)))


  if(NROW(Data_CurveFit2_woNA)>0){

    row.names(Data_CurveFit2_woNA)<-1:as.numeric(nrow(Data_CurveFit2_woNA))
    NumProteins<-as.numeric(nrow(Data_CurveFit2_woNA))
    NumTemps<-as.data.frame(Data_CurveFit2_woNA$NumberTemps)[1,]
    Protein<-1
    repeat{

      ProteinCount<-1

      repeat{

        if(ProteinCount==1){
          Temperature<-as.data.frame(Data_CurveFit2_woNA[Protein,(4+NumTemps):(3+NumTemps*2)])
          Abundance<-as.data.frame(Data_CurveFit2_woNA[Protein,(7+NumTemps*3):(6+NumTemps*4)])
          ProteinID<-Data_CurveFit2_woNA[Protein,1]
          T_Const<-Data_CurveFit2_woNA[Protein,(11+NumTemps*4)]
          B_Const<-Data_CurveFit2_woNA[Protein,(12+NumTemps*4)]
          xmid_Const<-Data_CurveFit2_woNA[Protein,(13+NumTemps*4)]
          scal_Const<-Data_CurveFit2_woNA[Protein,(14+NumTemps*4)]
          CurveFit_Predict<-B_Const+((T_Const-B_Const)/(1 + 10^(scal_Const*(xmid_Const-log10(Temperature))))^1)
          CurveFit_Actual<-t(Abundance)
          CurveFit_Mean<-mean(CurveFit_Actual)
          CurveFitData_Rsq<-data.frame(as.numeric(Temperature),as.numeric(CurveFit_Predict),as.numeric(CurveFit_Mean),as.numeric(CurveFit_Actual))
          ProteinStart<-Protein
          ProteinCount<-ProteinCount+1
          Protein<-Protein+1
        } else {
          if(ProteinCount>1 & ProteinID==Data_CurveFit2_woNA[Protein,1]){
            Temperature<-cbind(Temperature,as.data.frame(Data_CurveFit2_woNA[Protein,(4+NumTemps):(3+NumTemps*2)]))
            Abundance<-cbind(Abundance,as.data.frame(Data_CurveFit2_woNA[Protein,(7+NumTemps*3):(6+NumTemps*4)]))
            CurveFit_Predict<-B_Const+((T_Const-B_Const)/(1 + 10^(scal_Const*(xmid_Const-log10(Temperature))))^1)
            CurveFit_Actual<-t(Abundance)
            CurveFit_Mean<-mean(CurveFit_Actual)
            CurveFitData_Rsq<-data.frame(as.numeric(Temperature),as.numeric(CurveFit_Predict),as.numeric(CurveFit_Mean),as.numeric(CurveFit_Actual))
            ProteinCount<-ProteinCount+1
            Protein<-Protein+1

          } else {

            break
          }
        }
        if (Protein>NumProteins){
          break
        }

      }

      ProteinEnd<-Protein-1

      TempCount<-1
      RSE<-0
      TSE<-0
      repeat{
        RSE<-((CurveFitData_Rsq[TempCount,2]-CurveFitData_Rsq[TempCount,4])^2)+RSE
        TSE<-((CurveFitData_Rsq[TempCount,3]-CurveFitData_Rsq[TempCount,4])^2)+TSE
        TempCount<-TempCount+1
        if(TempCount>NROW(CurveFitData_Rsq)){
          break
        }
      }

      Rsq[ProteinStart:ProteinEnd,]<-1-(RSE/TSE)
      if (Protein>NumProteins){
        break
      }
    }
  }

  colnames(Rsq)<-"Rsq"
  Data_CurveFit2_PreComplete<-cbind(Data_CurveFit2_woNA,Rsq)
  Data_CurveFit2<-Data_CurveFit2_PreComplete[complete.cases(Data_CurveFit2_PreComplete),]

rm(CurveFitStats,Data_Corrected,Data_Corrected_Ordered,Data_CurveFit2_woNA,Data_CurveFit2_PreComplete)

  return(Data_CurveFit2)
}
