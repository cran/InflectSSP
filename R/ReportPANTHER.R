#' This function generates a PANTHER GO term based plot using the significant melt shifts from analysis
#' @param Data_Melts abundance and fit data for proteins
#' @param Directory directory where results are saved
#' @param Species species number for bioinformatics search
#' @param PANTHERpvalue pvalue for PANTHER analysis
#' @param PValMeltFDR Whether or not the FDR correction for pvalue is used in designation of melts of interest
#' @importFrom grDevices pdf
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 ggsave
#' @export
#' @return Excel files with summary of data along with melt curve plots for significant proteins
#' @examples
#'      \dontrun{
#'      ReportPANTHER(Data_Melts,Directory,Species,PANTHERpvalue,PValMeltFDR)}

ReportPANTHER<-function(Data_Melts,Directory,Species,PANTHERpvalue,PValMeltFDR){

  #Subsets proteins based on whether they meet Rsq, pValue and melt limit criteria
  Data_Melts_GoodRsq<-Data_Melts[Data_Melts$`Good Rsq`=="Yes",]
  Data_Melts_Ordered_GoodFit<-as.data.frame(Data_Melts_GoodRsq[order(Data_Melts_GoodRsq$`Melt Shift`),])
  Data_Melts_Ordered_GoodFit_GoodpVal<-Data_Melts_Ordered_GoodFit[Data_Melts_Ordered_GoodFit$`Good Melt p-value`=="Yes",]
  Data_Melts_Ordered_GoodFit_GoodpVal_GoodMelt<-Data_Melts_Ordered_GoodFit_GoodpVal[Data_Melts_Ordered_GoodFit_GoodpVal$`Good Melt Limit`=="Yes",]

  `Melt Shift`<-NULL
  term<-NULL
  `Melt p-value`<-NULL
  `PANTHER_pvalue`<-NULL
  `Melt_pValue_FDRAdj`<-NULL

  #Collects protein IDs that are queried in the PANTHER database
  if(NROW(Data_Melts_Ordered_GoodFit_GoodpVal_GoodMelt)>0){

    Reference<-as.data.frame(Data_Melts_Ordered_GoodFit$Accession)
    Data_Melts_Ordered_GoodFit_GoodpVal_GoodMelt_AbsOrder<-Data_Melts_Ordered_GoodFit_GoodpVal_GoodMelt[order(Data_Melts_Ordered_GoodFit_GoodpVal_GoodMelt$`Melt Shift`,decreasing=TRUE),]
    Accession<-as.data.frame(unique(Data_Melts_Ordered_GoodFit_GoodpVal_GoodMelt_AbsOrder$Accession))

    Count=1
    repeat{
      if(Count==1){
        Reference_Text<-Reference[1,]
        Count<-Count+1
      } else {
        Reference_Text<-paste(Reference_Text,Reference[Count,],sep=",")
        Count<-Count+1
      }
      if(Count>NROW(Reference)){
        break
      }
    }


    Count=1
    repeat{
      if(Count==1){
        Input_Text<-Accession[1,]
        Count<-Count+1
      } else {
        Input_Text<-paste(Input_Text,Accession[Count,],sep=",")
        Count<-Count+1
      }
      if(Count>NROW(Accession)){
        break
      }
    }

    #Queries the PANTHER database. ‘GO:0008150 biological_process’, ‘GO:0003674 molecular_function’, or ‘GO:0005575 cellular_component’.
    Output_BP<-httr::GET(paste(paste(paste(paste(paste(paste(paste(paste("http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=",Input_Text,sep=""),"&organism=",sep=""),Species,sep=""),"&refInputList=",sep=""),Reference_Text,sep=""),"&refOrganism=",sep=""),Species,sep=""),"&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR",sep=""))
    OutputAll_BP<-jsonlite::fromJSON(rawToChar(Output_BP$content))
    Data_BP<-data.frame(OutputAll_BP$results$result)

    Output_MF<-httr::GET(paste(paste(paste(paste(paste(paste(paste(paste("http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=",Input_Text,sep=""),"&organism=",sep=""),Species,sep=""),"&refInputList=",sep=""),Reference_Text,sep=""),"&refOrganism=",sep=""),Species,sep=""),"&annotDataSet=GO%3A0003674&enrichmentTestType=FISHER&correction=FDR",sep=""))
    OutputAll_MF<-jsonlite::fromJSON(rawToChar(Output_MF$content))
    Data_MF<-as.data.frame(OutputAll_MF$results$result)

    Output_CC<-httr::GET(paste(paste(paste(paste(paste(paste(paste(paste("http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=",Input_Text,sep=""),"&organism=",sep=""),Species,sep=""),"&refInputList=",sep=""),Reference_Text,sep=""),"&refOrganism=",sep=""),Species,sep=""),"&annotDataSet=GO%3A0005575&enrichmentTestType=FISHER&correction=FDR",sep=""))
    OutputAll_CC<-jsonlite::fromJSON(rawToChar(Output_CC$content))
    Data_CC<-as.data.frame(OutputAll_CC$results$result)

    Output_GOtoAcc<-httr::GET(paste(paste(paste("http://pantherdb.org//services/oai/pantherdb/geneinfo?geneInputList=",Input_Text,sep=""),"&organism=",sep=""),Species,sep=""))
    OutputAll_GOtoAcc<-jsonlite::fromJSON(rawToChar(Output_GOtoAcc$content))
    Data_GOtoAcc<-data.frame(OutputAll_GOtoAcc$search$mapped_genes$gene)

    if(length(OutputAll_GOtoAcc$search$error)>0){

    if(OutputAll_GOtoAcc$search$error == "Required parameter geneInputList has more items than supported"){
      print("Required parameter geneInputList has more items than supported. Increase stringency")
    }
    } else {

    GOtoAccRowCount<-as.numeric(nrow(Data_GOtoAcc))

    Count<-1
    repeat{

      GOtoAcc_Accession<-sub(".*UniProtKB=", "", Data_GOtoAcc$accession[Count])

      if(Count==1){
      if(length(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0003674",arr.ind = FALSE))>0){
      MFRowID<-as.numeric(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0003674",arr.ind = FALSE))
      if(is.null(ncol(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]]))==TRUE){
      GOtoAcc_GOTerm_MF_All<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[MFRowID]])
      } else {
        GOtoAcc_GOTerm_MF_All<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[MFRowID,1]])
      }
      colnames(GOtoAcc_GOTerm_MF_All)<-"term"
      GOtoAcc_GOTerm_MF_All[,2]<-GOtoAcc_Accession
      colnames(GOtoAcc_GOTerm_MF_All)[2]<-"Accession"
      }

      if(length(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0008150",arr.ind = FALSE))>0){
      BPRowID<-as.numeric(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0008150",arr.ind = FALSE))
      if(is.null(ncol(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]]))==TRUE){
        GOtoAcc_GOTerm_BP_All<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[BPRowID]])
      } else {
        GOtoAcc_GOTerm_BP_All<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[BPRowID,1]])
      }
      colnames(GOtoAcc_GOTerm_BP_All)<-"term"
      GOtoAcc_GOTerm_BP_All[,2]<-GOtoAcc_Accession
      colnames(GOtoAcc_GOTerm_BP_All)[2]<-"Accession"
      }

      if(length(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0005575",arr.ind = FALSE))>0){
      CCRowID<-as.numeric(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0005575",arr.ind = FALSE))
      if(is.null(ncol(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]]))==TRUE){
        GOtoAcc_GOTerm_CC_All<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[CCRowID]])
      } else {
        GOtoAcc_GOTerm_CC_All<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[CCRowID,1]])
      }
      colnames(GOtoAcc_GOTerm_CC_All)<-"term"
      GOtoAcc_GOTerm_CC_All[,2]<-GOtoAcc_Accession
      colnames(GOtoAcc_GOTerm_CC_All)[2]<-"Accession"
      }

      Count<-Count+1




      } else {

        GOtoAcc_Accession<-sub(".*UniProtKB=", "", Data_GOtoAcc$accession[Count])

        if(length(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0003674",arr.ind = FALSE))>0){
        MFRowID<-as.numeric(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0003674",arr.ind = FALSE))
        if(is.null(ncol(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]]))==TRUE){
          GOtoAcc_GOTerm_MF<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[MFRowID]])
        } else {
          GOtoAcc_GOTerm_MF<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[MFRowID,1]])
        }
        colnames(GOtoAcc_GOTerm_MF)<-"term"
        GOtoAcc_GOTerm_MF[,2]<-GOtoAcc_Accession
        colnames(GOtoAcc_GOTerm_MF)[2]<-"Accession"
        GOtoAcc_GOTerm_MF_All<-rbind(GOtoAcc_GOTerm_MF_All,GOtoAcc_GOTerm_MF)
        }

        if(length(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0008150",arr.ind = FALSE))>0){
        BPRowID<-as.numeric(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0008150",arr.ind = FALSE))
        if(is.null(ncol(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]]))==TRUE){
          GOtoAcc_GOTerm_BP<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[BPRowID]])
        } else {
          GOtoAcc_GOTerm_BP<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[BPRowID,1]])
        }
        colnames(GOtoAcc_GOTerm_BP)<-"term"
        GOtoAcc_GOTerm_BP[,2]<-GOtoAcc_Accession
        colnames(GOtoAcc_GOTerm_BP)[2]<-"Accession"
        GOtoAcc_GOTerm_BP_All<-rbind(GOtoAcc_GOTerm_BP_All,GOtoAcc_GOTerm_BP)
        }

        if(length(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0005575",arr.ind = FALSE))>0){
        CCRowID<-as.numeric(which(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["content"]]=="GO:0005575",arr.ind = FALSE))
        if(is.null(ncol(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]]))==TRUE){
          GOtoAcc_GOTerm_CC<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[CCRowID]])
        } else {
          GOtoAcc_GOTerm_CC<-data.frame(Data_GOtoAcc$annotation_type_list$annotation_data_type[[Count]][["annotation_list"]][["annotation"]][[CCRowID,1]])
        }
        colnames(GOtoAcc_GOTerm_CC)<-"term"
        GOtoAcc_GOTerm_CC[,2]<-GOtoAcc_Accession
        colnames(GOtoAcc_GOTerm_CC)[2]<-"Accession"
        GOtoAcc_GOTerm_CC_All<-rbind(GOtoAcc_GOTerm_CC_All,GOtoAcc_GOTerm_CC)
        }

          Count<-Count+1

        }

      if(Count>GOtoAccRowCount){
        break
      }

      }


  Data_BP_PValnumeric<-data.frame(as.numeric(Data_BP$pValue))
  Data_BP_Subset<-cbind(Data_BP$fdr,cbind(Data_BP$term$label,Data_BP_PValnumeric))
  colnames(Data_BP_Subset)<-c("fdr","term","PANTHER_pvalue")
  Data_BP_Merged<-merge(Data_BP_Subset,GOtoAcc_GOTerm_BP_All,by="term")
  Data_BP_Merged_WMelts<-merge(Data_Melts_Ordered_GoodFit_GoodpVal_GoodMelt,Data_BP_Merged,by="Accession")
  Data_BP_Melts_Subset<-subset(Data_BP_Merged_WMelts[Data_BP_Merged_WMelts$PANTHER_pvalue<=PANTHERpvalue,])

  write.csv(Data_BP_Melts_Subset, file = paste(Directory,paste("Result Files","Data_BP_Melts_Subset.csv",sep="/"),sep="/"))

  Data_MF_PValnumeric<-data.frame(as.numeric(Data_MF$pValue))
  Data_MF_Subset<-cbind(Data_MF$fdr,cbind(Data_MF$term$label,Data_MF_PValnumeric))
  colnames(Data_MF_Subset)<-c("fdr","term","PANTHER_pvalue")
  Data_MF_Merged<-merge(Data_MF_Subset,GOtoAcc_GOTerm_MF_All,by="term")
  Data_MF_Merged_WMelts<-merge(Data_Melts_Ordered_GoodFit_GoodpVal_GoodMelt,Data_MF_Merged,by="Accession")
  Data_MF_Melts_Subset<-subset(Data_MF_Merged_WMelts[Data_MF_Merged_WMelts$PANTHER_pvalue<=PANTHERpvalue,])

  write.csv(Data_MF_Melts_Subset, file = paste(Directory,paste("Result Files","Data_MF_Melts_Subset.csv",sep="/"),sep="/"))

  Data_CC_PValnumeric<-data.frame(as.numeric(Data_CC$pValue))
  Data_CC_Subset<-cbind(Data_CC$fdr,cbind(Data_CC$term$label,Data_CC_PValnumeric))
  colnames(Data_CC_Subset)<-c("fdr","term","PANTHER_pvalue")
  Data_CC_Merged<-merge(Data_CC_Subset,GOtoAcc_GOTerm_CC_All,by="term")
  Data_CC_Merged_WMelts<-merge(Data_Melts_Ordered_GoodFit_GoodpVal_GoodMelt,Data_CC_Merged,by="Accession")
  Data_CC_Melts_Subset<-subset(Data_CC_Merged_WMelts[Data_CC_Merged_WMelts$PANTHER_pvalue<=PANTHERpvalue,])

  write.csv(Data_CC_Melts_Subset, file = paste(Directory,paste("Result Files","Data_CC_Melts_Subset.csv",sep="/"),sep="/"))

    #Creates plots that use both melt shift data and PANTHER GO term data

  if(PValMeltFDR=="Yes"){

    if(NROW(Data_BP_Melts_Subset)>0){

      MaxBP<-1.1*max(abs(Data_BP_Melts_Subset$`Melt Shift`))

      ggplot2::ggplot(Data_BP_Melts_Subset, ggplot2::aes(x=`Melt Shift`, y=term, size=Melt_pValue_FDRAdj, color=PANTHER_pvalue),) +
      ggplot2::geom_point(alpha=0.7)+
      ggplot2::xlab("Melt Shift, Condition-Control, (C)")+
      ggplot2::ylab("PANTHER Biological Process Terms")+
      ggplot2::scale_x_continuous(limits=c(-(MaxBP),MaxBP))+
      ggplot2::scale_size_continuous(range = c(6,1))
      ggplot2::ggsave(paste(Directory,paste("Result Files","PANTHER_BP.svg",sep="/"),sep="/"),width=7, height=10)

    } else {

      print("Insufficient number of terms for complete Biological Process analysis")
    }

    if(NROW(Data_MF_Melts_Subset)>0){

      MaxMF<-1.1*max(abs(Data_MF_Melts_Subset$`Melt Shift`))

      ggplot2::ggplot(Data_MF_Melts_Subset, ggplot2::aes(x=`Melt Shift`, y=term, size=Melt_pValue_FDRAdj, color=PANTHER_pvalue),) +
        ggplot2::geom_point(alpha=0.7)+
        ggplot2::xlab("Melt Shift, Condition-Control, (C)")+
        ggplot2::ylab("PANTHER Molecular Function Terms")+
        ggplot2::scale_x_continuous(limits=c(-(MaxMF),MaxMF))+
        ggplot2::scale_size_continuous(range = c(6,1))
      ggplot2::ggsave(paste(Directory,paste("Result Files","PANTHER_MF.svg",sep="/"),sep="/"),width=7, height=10)

    } else {

      print("Insufficient number of terms for complete Molecular Function analysis")
    }


    if(NROW(Data_CC_Melts_Subset)>0){

      MaxCC<-1.1*max(abs(Data_CC_Melts_Subset$`Melt Shift`))
      ggplot2::ggplot(Data_CC_Melts_Subset, ggplot2::aes(x=`Melt Shift`, y=term, size=Melt_pValue_FDRAdj, color=PANTHER_pvalue),) +
        ggplot2::geom_point(alpha=0.7)+
        ggplot2::xlab("Melt Shift, Condition-Control, (C)")+
        ggplot2::ylab("PANTHER Cellular Component Terms")+
        ggplot2::scale_x_continuous(limits=c(-(MaxCC),MaxCC))+
        ggplot2::scale_size_continuous(range = c(6,1))
      ggplot2::ggsave(paste(Directory,paste("Result Files","PANTHER_CC.svg",sep="/"),sep="/"),width=7, height=10)

    } else {

      print("Insufficient number of terms for complete Cellular Component analysis")
    }


  } else {

    if(NROW(Data_BP_Melts_Subset)>0){

      MaxBP<-1.1*max(abs(Data_BP_Melts_Subset$`Melt Shift`))

      ggplot2::ggplot(Data_BP_Melts_Subset, ggplot2::aes(x=`Melt Shift`, y=term, size=`Melt p-value`, color=PANTHER_pvalue),) +
        ggplot2::geom_point(alpha=0.7)+
        ggplot2::xlab("Melt Shift, Condition-Control, (C)")+
        ggplot2::ylab("PANTHER Biological Process Terms")+
        ggplot2::scale_x_continuous(limits=c(-(MaxBP),MaxBP))+
        ggplot2::scale_size_continuous(range = c(6,1))
      ggplot2::ggsave(paste(Directory,paste("Result Files","PANTHER_BP.svg",sep="/"),sep="/"),width=7, height=10)

    } else {

      print("Insufficient number of terms for complete Biological Process analysis")
    }

    if(NROW(Data_MF_Melts_Subset)>0){

      MaxMF<-1.1*max(abs(Data_MF_Melts_Subset$`Melt Shift`))

      ggplot2::ggplot(Data_MF_Melts_Subset, ggplot2::aes(x=`Melt Shift`, y=term, size=`Melt p-value`, color=PANTHER_pvalue),) +
        ggplot2::geom_point(alpha=0.7)+
        ggplot2::xlab("Melt Shift, Condition-Control, (C)")+
        ggplot2::ylab("PANTHER Molecular Function Terms")+
        ggplot2::scale_x_continuous(limits=c(-(MaxMF),MaxMF))+
        ggplot2::scale_size_continuous(range = c(6,1))
      ggplot2::ggsave(paste(Directory,paste("Result Files","PANTHER_MF.svg",sep="/"),sep="/"),width=7, height=10)

    } else {

      print("Insufficient number of terms for complete Molecular Function analysis")
    }


    if(NROW(Data_CC_Melts_Subset)>0){

      MaxCC<-1.1*max(abs(Data_CC_Melts_Subset$`Melt Shift`))
      ggplot2::ggplot(Data_CC_Melts_Subset, ggplot2::aes(x=`Melt Shift`, y=term, size=`Melt p-value`, color=PANTHER_pvalue),) +
        ggplot2::geom_point(alpha=0.7)+
        ggplot2::xlab("Melt Shift, Condition-Control, (C)")+
        ggplot2::ylab("PANTHER Cellular Component Terms")+
        ggplot2::scale_x_continuous(limits=c(-(MaxCC),MaxCC))+
        ggplot2::scale_size_continuous(range = c(6,1))
      ggplot2::ggsave(paste(Directory,paste("Result Files","PANTHER_CC.svg",sep="/"),sep="/"),width=7, height=10)

    } else {

      print("Insufficient number of terms for complete Cellular Component analysis")
    }

  }


}



  } else {

    print("Insufficient number of proteins with good fit to complete PANTHER analysis")

  }



}


