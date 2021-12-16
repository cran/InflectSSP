#' This function generates a STRING based network using the significant melt shifts from analysis
#' @param Data_Melts abundance and fit data for proteins that meet quality criteria in overall workflow
#' @param STRINGScore the STRING score that is used to determine whether an interaction is significant
#' @param Directory directory where results are saved
#' @param Species species taxon number for bioinformatics search
#' @importFrom grDevices pdf
#' @importFrom httr GET
#' @importFrom jsonlite fromJSON
#' @importFrom GGally ggnet2
#' @importFrom network network
#' @importFrom network %v%
#' @importFrom ggplot2 ggsave
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 cut_number
#' @export
#' @return Excel files with summary of data along with melt curve plots for significant proteins
#' @examples
#' \dontrun{
#' ReportSTRING(Data_Melts,STRINGScore,Directory,Species)
#' }


ReportSTRING<-function(Data_Melts,STRINGScore,Directory,Species){


  if(NROW(Data_Melts)>0){
    #Subsets proteins based on whether they meet Rsq, pValue and melt limit criteria
    Data_Melts_GoodRsq<-Data_Melts[Data_Melts$`Good Rsq`=="Yes",]
    Data_Melts_GoodRsq_GoodpVal<-Data_Melts_GoodRsq[Data_Melts_GoodRsq$`Good Melt p-value`=="Yes",]
    Data_Melts_GoodRsq_GoodpVal_GoodMelt<-Data_Melts_GoodRsq_GoodpVal[Data_Melts_GoodRsq_GoodpVal$`Good Melt Limit`=="Yes",]

    Accession<-as.data.frame(Data_Melts_GoodRsq_GoodpVal_GoodMelt$Accession)
    Melt<-as.data.frame(Data_Melts_GoodRsq_GoodpVal_GoodMelt$`Melt Shift`)
    Meltp<-Data_Melts_GoodRsq_GoodpVal_GoodMelt$Melt_pValue_FDRAdj
    Accession_Melt<-as.data.frame(cbind(cbind(Accession,Melt),Meltp))
    colnames(Accession_Melt)<-c("Accession","Melt","Melt_pValue_FDRAdj")

    Count<-1
    repeat{
      if(Count==1){
        Output1<-httr::GET(paste(paste(paste(paste("https://string-db.org/api/json/interaction_partners?identifiers=",Accession[Count,],sep=""),"&species=",sep=""),Species,sep=""),"&limit=1000",sep=""))
        if(Output1$status_code==200){
          OutputAll<-jsonlite::fromJSON(rawToChar(Output1$content))
          MeltValue<-as.data.frame(Accession_Melt[Accession==Accession[Count,],2])
          MeltpValue<-as.data.frame(Accession_Melt[Accession==Accession[Count,],3])
          OutputAll[,14]<-MeltValue[1,1]
          OutputAll[,15]<-MeltpValue[1,1]
        } else {

        }

      } else {
        Output<-httr::GET(paste(paste(paste(paste("https://string-db.org/api/json/interaction_partners?identifiers=",Accession[Count,],sep=""),"&species=",sep=""),Species,sep=""),"&limit=1000",sep=""))
        if(Output$status_code==200){
          Output_values <- jsonlite::fromJSON(rawToChar(Output$content))
          MeltValue<-as.data.frame(Accession_Melt[Accession==Accession[Count,],2])
          MeltpValue<-as.data.frame(Accession_Melt[Accession==Accession[Count,],3])
          Output_values[,14]<-MeltValue[1,1]
          Output_values[,15]<-MeltpValue[1,1]
          OutputAll<-rbind(OutputAll,Output_values)

        } else {

        }

      }
      Count<-Count+1
      if(Count>NROW(Accession)){
        break
      }
    }

    OutputAll_scoresubset<-OutputAll[OutputAll$score>=STRINGScore,]

    if(NROW(OutputAll_scoresubset)>0){
      Contains1<-OutputAll_scoresubset[,4] %in% OutputAll_scoresubset[,3]
      OutputAll_scoresubset[Contains1,16]<-"Yes"
      OutputAll_scoresubset_2<-subset(OutputAll_scoresubset,OutputAll_scoresubset[,16]=="Yes")
      if(NROW(OutputAll_scoresubset_2)>0){

        Contains2<-OutputAll_scoresubset_2[,3] %in% OutputAll_scoresubset_2[,4]
        OutputAll_scoresubset_2[Contains2,17]<-"Yes"

        OutputAll_Interactors<-subset(OutputAll_scoresubset_2,OutputAll_scoresubset_2[,17]=="Yes")

        OutputAllMatrix<-tapply(OutputAll_Interactors$score, OutputAll_Interactors[c("preferredName_A", "preferredName_B")], mean)
        OutputAllMatrix[is.na(OutputAllMatrix)] <- 0
        OutputAllMatrix[OutputAllMatrix>=STRINGScore] <- 1


        net<-network::network(OutputAllMatrix, directed=FALSE)

        Nodes = as.data.frame(unique(OutputAll_Interactors$preferredName_A))
        colnames(Nodes)<-"Accession"
        OutputAll_Interactors_2 <- unique(cbind(cbind(OutputAll_Interactors[,3],OutputAll_Interactors[,14]),OutputAll_Interactors[,15]))
        colnames(OutputAll_Interactors_2)<-c("Accession","Melt","MeltpValue")

        NodeMelts<-as.data.frame(merge(Nodes,OutputAll_Interactors_2,by="Accession"))
        network::`%v%`(net,"MeltShift") <- round(as.numeric(NodeMelts$Melt),2)


        network::`%v%`(net,"range")  <- ifelse(net %v% "MeltShift" <= (-5),"Melt<=-5",ifelse(net %v% "MeltShift" <= (-2.5),"-5<Melt<=-2.5",ifelse(net %v% "MeltShift" <= (0),"2.5<Melt<=0",ifelse(net %v% "MeltShift" <= (2.5),"0<Melt<=2.5",ifelse(net %v% "MeltShift" < (5),"2.5<Melt<5","Melt>=5")))))

        GGally::ggnet2(net,label=TRUE,label.size=5,node.size=5,color="range",color.legend="Melt Shift (C)",palette = c("Melt<=-5" = "firebrick4", "-5<Melt<=-2.5" = "firebrick1","2.5<Melt<=0" = "pink","0<Melt<=2.5" = "lightblue1", "2.5<Melt<5"="dodgerblue","Melt>=5"="darkblue"),edge.size=.5)
        ggplot2::ggsave(paste(Directory,paste("Result Files","STRINGNetwork.svg",sep="/"),sep="/"))


        write.csv(NodeMelts, file = paste(Directory,paste("Result Files","NodeMelts.csv",sep="/"),sep="/"))



      } else {

        print("STRING analysis not possible due to lack of interactions")

      }

    } else {
      print("Insufficient number of proteins with good fit to complete STRING analysis")
    }

  } else {

    print("Insufficient number of proteins with good fit to complete STRING analysis")
  }

}
