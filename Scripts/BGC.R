###################################### Supplementary Figure5 a ######################################################
setwd("D:/metagenome/New bigscape")
library(ggplot2)
library(stringr)
library(Hmisc)
sj_type<-c("Hot","Pressure","Saline_alkaline","Acid","Cold")
bgclist<-c("NRPS","PKSI","PKS-NRP_Hybrids","PKSother","RiPPs","Saccharides","Terpene","Others")

data_all<-matrix(nrow=1,ncol=4)
colnames(data_all)<-c("Habitat","NMDC_ID","Region_filename","BGC_group")
for(mm in 1:length(sj_type)){
  data<-read.table(paste0("new_",sj_type[mm],"_all_antismash_bigscape_region_length_morethan5kb_addNMDC.txt"),sep="\t",header=T)
  data$Type[which(data$Type=="New")]<-"Unannotated"
  data3<-cbind(matrix(rep(sj_type[mm],nrow(data)),ncol=1),data[,c("NMDC_ID","Region_filename","BiG.SCAPE_class")])
  colnames(data3)<-colnames(data_all)
  data_all<-rbind(data_all,data3)
}
data_all<-data_all[-1,]
data_all$info=paste(data_all$NMDC_ID,data_all$Region_filename)

data_all3<-aggregate(data_all$info,by=list(data_all$Habitat,data_all$BGC_group),length)
colnames(data_all3)<-c("Habitat","BGC_group","Number_of_BGCs")
data_all3$BGC_group<-factor(data_all3$BGC_group,levels = c("NRPS","PKSI","PKS-NRP_Hybrids","PKSother","RiPPs","Saccharides","Terpene","Others"))

palette<-c("NRPS"="#377EB8","PKSI"="#4DAF4A","PKS-NRP_Hybrids"="#984EA3","PKSother"="#FF7F00","RiPPs"="#E41A1C","Saccharides"="#A65628","Terpene"="#F781BF","Others"="#999999")
ggplot(data_all3,aes(x=Habitat,y=Number_of_BGCs,fill=BGC_group))+geom_bar(stat = "identity",position="stack")+
  theme_classic()+ylab("Number of BGCs")+
  scale_fill_manual(values =palette)+
  theme(axis.line=element_line(color="black"),axis.text=element_text(size=10,color="black"),axis.title.x = element_blank(),axis.title=element_text(size=12,color="black"),legend.text=element_text(size=10,color="black"),legend.title=element_text(size=10,color="black"))
ggsave(filename = "BGC_regionNum_in_all_habitats_barplot_20240409.pdf",width=9,height=6)




###################################### Supplementary Figure5 bcdef ######################################################
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(patchwork)

setwd("D:/工作/metagenome/New bigscape")
sj_type<-c("Hot","Pressure","Saline_alkaline","Acid","Cold")
sj_type2<-c("Hot","Pressure","Saline_alkaline","Acid","Cold")
all_all<-matrix(ncol=4,nrow=1)
colnames(all_all)<-c("Habitat","Phyla","Type","Number_of_proteins")
all_all2<-matrix(ncol=4,nrow=1)
colnames(all_all2)<-c("Habitat","Phyla","BiG.SCAPE_class","protein_number")

####### hot #######
for(mm in 1){
  data<-read.table(paste0("new_",sj_type[mm],"_all_antismash_bigscape_region_length_morethan5kb_addNMDC.txt"),header=T,sep="\t",quote='')
  data$Type[which(data$Type=="New")]<-"Unannotated"
  gtdb<-gsub(";","",gsub("p__","",str_extract(data$GTDB,"p__(.*?);")))
  data2<-cbind(data,gtdb)
  if(length(unique(data2$gtdb))>=10){
    needlist<-names(sort(table(data2$gtdb),decreasing = T))[1:10]
  }else{
    needlist<-names(sort(table(data2$gtdb),decreasing = T))
  }
  data3<-data2[which(!is.na(match(data2$gtdb,needlist))),]
  data3$BiG.SCAPE_class<-as.factor(data3$BiG.SCAPE_class)

  data_plot<-aggregate(data3$Region_filename,by=list(data3$gtdb,data3$BiG.SCAPE_class),length)
  colnames(data_plot)<-c("gtdb","BiG.SCAPE_class","region_number")
  bgclist<-rev(c("NRPS","PKSI","PKS-NRP_Hybrids","PKSother","RiPPs","Saccharides","Terpene","Others"))


  data4<-unique(data3[,c("Region_filename","gtdb")])
  gtdbnumber2<-aggregate(data4$Region_filename,by=list(data4$gtdb),length)
  gtdbnumber3<-gtdbnumber2[order(as.numeric(gtdbnumber2[,2]),decreasing = T),]
  colnames(gtdbnumber3)<-c("GTDB","Number_of_MAGs")
  gtdbnumber4<-matrix(nrow=nrow(gtdbnumber3)*2,ncol=3)
  gtdbnumber4[,1]<-rep(gtdbnumber3$GTDB,each=2)
  gtdbnumber4[,2]<-rep(c("Annotated","Unannotated"),times=nrow(gtdbnumber3))
  colnames(gtdbnumber4)<-c("Phyla","Type","Number_of_regions")
  for(i in 1:nrow(gtdbnumber4)){
    index<-which(data3$Type==gtdbnumber4[i,2] & data3$gtdb==gtdbnumber4[i,1])
    if(length(index)>0){
      gtdbnumber4[i,3]<-length(unique(data3[index,"Region_filename"]))
    }else{
      gtdbnumber4[i,3]<-0
    }
  }
  gtdbnumber4<-as.data.frame(gtdbnumber4)
  gtdbnumber4$Phyla<-factor(gtdbnumber4$Phyla,levels=rev(gtdbnumber3$GTDB))
  gtdbnumber4$Type<-factor(gtdbnumber4$Type,levels = rev(c("Annotated","Unannotated")))
  gtdbnumber4$Number_of_regions<-as.numeric(gtdbnumber4$Number_of_regions)

  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(gtdbnumber4)),ncol=1),gtdbnumber4)
  colnames(all_all0)<-colnames(all_all)
  all_all<-rbind(all_all,all_all0)

  for (i in seq(1, 2)) {
    for(j in 1:length(unique(gtdbnumber4[,1]))){
      if (i == 1) {
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = 1
      }else{
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = sum(gtdbnumber4$Number_of_regions[(2*(j-1)+1):(2*(j-1)+i-1)]) + gtdbnumber4$Number_of_regions[2*(j-1)+i]/4
      }
    }
  }
  p1<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    theme(axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))


  data_plot2<-matrix(nrow=length(needlist)*length(bgclist),ncol=ncol(data_plot))
  colnames(data_plot2)<-colnames(data_plot)
  data_plot2[,1]<-rep(as.character(gtdbnumber3$GTDB),times=length(bgclist))
  data_plot2[,2]<-rep(bgclist,each=length(needlist))
  for(i in 1:nrow(data_plot2)){
    index<-which(data_plot$gtdb==data_plot2[i,1] & data_plot$BiG.SCAPE_class==data_plot2[i,2])
    if(length(index)>0){
      data_plot2[i,3]<-data_plot[index,3]
    }else{
      data_plot2[i,3]<-0
    }
  }
  data_plot2<-as.data.frame(data_plot2)
  data_plot2$region_number<-as.numeric(data_plot2$region_number)

  gtdbnumber<-aggregate(data_plot2$region_number,by=list(data_plot2$gtdb),sum)
  for(i in 1:nrow(data_plot2)){
    data_plot2[i,1]<-paste0(data_plot2[i,1]," (",gtdbnumber[which(gtdbnumber[,1]==data_plot2[i,1]),2],")")
  }

  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(data_plot2)),ncol=1),data_plot2)
  colnames(all_all0)<-colnames(all_all2)
  all_all2<-rbind(all_all2,all_all0)

  data_plot2$gtdb<-factor(data_plot2$gtdb,levels = data_plot2$gtdb[length(unique(data_plot2$gtdb)):1])
  p2<-ggplot(data_plot2,aes(x=gtdb,y=region_number,fill=factor(BiG.SCAPE_class,levels=bgclist)))+geom_col(position="fill",width=1,color="black")+
    labs(x="Phyla",y="BGC content(%)",fill="Group")+
    theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=25,color="black"),legend.text=element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_y_continuous(breaks=seq(0,1,0.25),labels=c('0','25','50','75','100'),expand=c(0,0))+
    coord_flip()+
    scale_fill_manual(values=c("NRPS"="#377EB8","PKSI"="#4DAF4A","PKS-NRP_Hybrids"="#984EA3","PKSother"="#FF7F00","RiPPs"="#E41A1C","Saccharides"="#A65628","Terpene"="#F781BF","Others"="#999999"))

  p1_2<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    xlab(NULL)+theme(axis.line=element_line(color="black"),axis.line.y=element_blank(),axis.ticks.y = element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))
  p3_1<-p2+p1_2+ plot_layout(guides = 'collect')+ plot_annotation(title = sj_type2[mm]) & theme(plot.title = element_text(size = 27,hjust=0.5))
  p4_1<-p2
}

####### pressure #######
for(mm in 2){
  data<-read.table(paste0("new_",sj_type[mm],"_all_antismash_bigscape_region_length_morethan5kb_addNMDC.txt"),header=T,sep="\t",quote='')
  data$Type[which(data$Type=="New")]<-"Unannotated"
  gtdb<-gsub(";","",gsub("p__","",str_extract(data$GTDB,"p__(.*?);")))
  data2<-cbind(data,gtdb)
  if(length(unique(data2$gtdb))>=10){
    needlist<-names(sort(table(data2$gtdb),decreasing = T))[1:10]
  }else{
    needlist<-names(sort(table(data2$gtdb),decreasing = T))
  }
  data3<-data2[which(!is.na(match(data2$gtdb,needlist))),]
  data3$BiG.SCAPE_class<-as.factor(data3$BiG.SCAPE_class)

  data_plot<-aggregate(data3$Region_filename,by=list(data3$gtdb,data3$BiG.SCAPE_class),length)
  colnames(data_plot)<-c("gtdb","BiG.SCAPE_class","region_number")
  bgclist<-rev(c("NRPS","PKSI","PKS-NRP_Hybrids","PKSother","RiPPs","Saccharides","Terpene","Others"))

  data4<-unique(data3[,c("Region_filename","gtdb")])
  gtdbnumber2<-aggregate(data4$Region_filename,by=list(data4$gtdb),length)
  gtdbnumber3<-gtdbnumber2[order(as.numeric(gtdbnumber2[,2]),decreasing = T),]
  colnames(gtdbnumber3)<-c("GTDB","Number_of_MAGs")
  gtdbnumber4<-matrix(nrow=nrow(gtdbnumber3)*2,ncol=3)
  gtdbnumber4[,1]<-rep(gtdbnumber3$GTDB,each=2)
  gtdbnumber4[,2]<-rep(c("Annotated","Unannotated"),times=nrow(gtdbnumber3))
  colnames(gtdbnumber4)<-c("Phyla","Type","Number_of_regions")
  for(i in 1:nrow(gtdbnumber4)){
    index<-which(data3$Type==gtdbnumber4[i,2] & data3$gtdb==gtdbnumber4[i,1])
    if(length(index)>0){
      gtdbnumber4[i,3]<-length(unique(data3[index,"Region_filename"]))
    }else{
      gtdbnumber4[i,3]<-0
    }
  }
  gtdbnumber4<-as.data.frame(gtdbnumber4)
  gtdbnumber4$Phyla<-factor(gtdbnumber4$Phyla,levels=rev(gtdbnumber3$GTDB))
  gtdbnumber4$Type<-factor(gtdbnumber4$Type,levels = rev(c("Annotated","Unannotated")))
  gtdbnumber4$Number_of_regions<-as.numeric(gtdbnumber4$Number_of_regions)

  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(gtdbnumber4)),ncol=1),gtdbnumber4)
  colnames(all_all0)<-colnames(all_all)
  all_all<-rbind(all_all,all_all0)

  for (i in seq(1, 2)) {
    for(j in 1:length(unique(gtdbnumber4[,1]))){
      if (i == 1) {
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = 1
      }else{
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = sum(gtdbnumber4$Number_of_regions[(2*(j-1)+1):(2*(j-1)+i-1)]) + gtdbnumber4$Number_of_regions[2*(j-1)+i]/4
      }
    }
  }
  p1<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    theme(axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))

  data_plot2<-matrix(nrow=length(needlist)*length(bgclist),ncol=ncol(data_plot))
  colnames(data_plot2)<-colnames(data_plot)
  data_plot2[,1]<-rep(as.character(gtdbnumber3$GTDB),times=length(bgclist))
  data_plot2[,2]<-rep(bgclist,each=length(needlist))
  for(i in 1:nrow(data_plot2)){
    index<-which(data_plot$gtdb==data_plot2[i,1] & data_plot$BiG.SCAPE_class==data_plot2[i,2])
    if(length(index)>0){
      data_plot2[i,3]<-data_plot[index,3]
    }else{
      data_plot2[i,3]<-0
    }
  }
  data_plot2<-as.data.frame(data_plot2)
  data_plot2$region_number<-as.numeric(data_plot2$region_number)

  gtdbnumber<-aggregate(data_plot2$region_number,by=list(data_plot2$gtdb),sum)
  for(i in 1:nrow(data_plot2)){
    data_plot2[i,1]<-paste0(data_plot2[i,1]," (",gtdbnumber[which(gtdbnumber[,1]==data_plot2[i,1]),2],")")
  }

  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(data_plot2)),ncol=1),data_plot2)
  colnames(all_all0)<-colnames(all_all2)
  all_all2<-rbind(all_all2,all_all0)

  data_plot2$gtdb<-factor(data_plot2$gtdb,levels = data_plot2$gtdb[length(unique(data_plot2$gtdb)):1])
  p2<-ggplot(data_plot2,aes(x=gtdb,y=region_number,fill=factor(BiG.SCAPE_class,levels=bgclist)))+geom_col(position="fill",width=1,color="black")+
    labs(x="Phyla",y="BGC content(%)",fill="Group")+
    theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=25,color="black"),legend.text=element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_y_continuous(breaks=seq(0,1,0.25),labels=c('0','25','50','75','100'),expand=c(0,0))+
    coord_flip()+
    scale_fill_manual(values=c("NRPS"="#377EB8","PKSI"="#4DAF4A","PKS-NRP_Hybrids"="#984EA3","PKSother"="#FF7F00","RiPPs"="#E41A1C","Saccharides"="#A65628","Terpene"="#F781BF","Others"="#999999"))

  p1_2<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    xlab(NULL)+theme(axis.line=element_line(color="black"),axis.line.y=element_blank(),axis.ticks.y = element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))
  p3_2<-p2+p1_2+ plot_layout(guides = 'collect')+ plot_annotation(title = sj_type2[mm]) & theme(plot.title = element_text(size = 27,hjust=0.5))
  p4_2<-p2
}


####### saline_alkaline #######
for(mm in 3){
  data<-read.table(paste0("new_",sj_type[mm],"_all_antismash_bigscape_region_length_morethan5kb_addNMDC.txt"),header=T,sep="\t",quote='')
  data$Type[which(data$Type=="New")]<-"Unannotated"
  gtdb<-gsub(";","",gsub("p__","",str_extract(data$GTDB,"p__(.*?);")))
  data2<-cbind(data,gtdb)
  if(length(unique(data2$gtdb))>=10){
    needlist<-names(sort(table(data2$gtdb),decreasing = T))[1:10]
  }else{
    needlist<-names(sort(table(data2$gtdb),decreasing = T))
  }
  data3<-data2[which(!is.na(match(data2$gtdb,needlist))),]
  data3$BiG.SCAPE_class<-as.factor(data3$BiG.SCAPE_class)

  data_plot<-aggregate(data3$Region_filename,by=list(data3$gtdb,data3$BiG.SCAPE_class),length)
  colnames(data_plot)<-c("gtdb","BiG.SCAPE_class","region_number")
  bgclist<-rev(c("NRPS","PKSI","PKS-NRP_Hybrids","PKSother","RiPPs","Saccharides","Terpene","Others"))

  data4<-unique(data3[,c("Region_filename","gtdb")])
  gtdbnumber2<-aggregate(data4$Region_filename,by=list(data4$gtdb),length)
  gtdbnumber3<-gtdbnumber2[order(as.numeric(gtdbnumber2[,2]),decreasing = T),]
  colnames(gtdbnumber3)<-c("GTDB","Number_of_MAGs")
  gtdbnumber4<-matrix(nrow=nrow(gtdbnumber3)*2,ncol=3)
  gtdbnumber4[,1]<-rep(gtdbnumber3$GTDB,each=2)
  gtdbnumber4[,2]<-rep(c("Annotated","Unannotated"),times=nrow(gtdbnumber3))
  colnames(gtdbnumber4)<-c("Phyla","Type","Number_of_regions")
  for(i in 1:nrow(gtdbnumber4)){
    index<-which(data3$Type==gtdbnumber4[i,2] & data3$gtdb==gtdbnumber4[i,1])
    if(length(index)>0){
      gtdbnumber4[i,3]<-length(unique(data3[index,"Region_filename"]))
    }else{
      gtdbnumber4[i,3]<-0
    }
  }
  gtdbnumber4<-as.data.frame(gtdbnumber4)
  gtdbnumber4$Phyla<-factor(gtdbnumber4$Phyla,levels=rev(gtdbnumber3$GTDB))
  gtdbnumber4$Type<-factor(gtdbnumber4$Type,levels = rev(c("Annotated","Unannotated")))
  gtdbnumber4$Number_of_regions<-as.numeric(gtdbnumber4$Number_of_regions)

  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(gtdbnumber4)),ncol=1),gtdbnumber4)
  colnames(all_all0)<-colnames(all_all)
  all_all<-rbind(all_all,all_all0)

  for (i in seq(1, 2)) {
    for(j in 1:length(unique(gtdbnumber4[,1]))){
      if (i == 1) {
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = 1
      }else{
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = sum(gtdbnumber4$Number_of_regions[(2*(j-1)+1):(2*(j-1)+i-1)]) + gtdbnumber4$Number_of_regions[2*(j-1)+i]/4
      }
    }
  }
  p1<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    theme(axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))

  data_plot2<-matrix(nrow=length(needlist)*length(bgclist),ncol=ncol(data_plot))
  colnames(data_plot2)<-colnames(data_plot)
  data_plot2[,1]<-rep(as.character(gtdbnumber3$GTDB),times=length(bgclist))
  data_plot2[,2]<-rep(bgclist,each=length(needlist))
  for(i in 1:nrow(data_plot2)){
    index<-which(data_plot$gtdb==data_plot2[i,1] & data_plot$BiG.SCAPE_class==data_plot2[i,2])
    if(length(index)>0){
      data_plot2[i,3]<-data_plot[index,3]
    }else{
      data_plot2[i,3]<-0
    }
  }
  data_plot2<-as.data.frame(data_plot2)
  data_plot2$region_number<-as.numeric(data_plot2$region_number)

  gtdbnumber<-aggregate(data_plot2$region_number,by=list(data_plot2$gtdb),sum)
  for(i in 1:nrow(data_plot2)){
    data_plot2[i,1]<-paste0(data_plot2[i,1]," (",gtdbnumber[which(gtdbnumber[,1]==data_plot2[i,1]),2],")")
  }

  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(data_plot2)),ncol=1),data_plot2)
  colnames(all_all0)<-colnames(all_all2)
  all_all2<-rbind(all_all2,all_all0)

  data_plot2$gtdb<-factor(data_plot2$gtdb,levels = data_plot2$gtdb[length(unique(data_plot2$gtdb)):1])
  p2<-ggplot(data_plot2,aes(x=gtdb,y=region_number,fill=factor(BiG.SCAPE_class,levels=bgclist)))+geom_col(position="fill",width=1,color="black")+
    labs(x="Phyla",y="BGC content(%)",fill="Group")+
    theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=25,color="black"),legend.text=element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_y_continuous(breaks=seq(0,1,0.25),labels=c('0','25','50','75','100'),expand=c(0,0))+
    coord_flip()+
    scale_fill_manual(values=c("NRPS"="#377EB8","PKSI"="#4DAF4A","PKS-NRP_Hybrids"="#984EA3","PKSother"="#FF7F00","RiPPs"="#E41A1C","Saccharides"="#A65628","Terpene"="#F781BF","Others"="#999999"))

  p1_2<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    xlab(NULL)+theme(axis.line=element_line(color="black"),axis.line.y=element_blank(),axis.ticks.y = element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))
  p3_3<-p2+p1_2+ plot_layout(guides = 'collect')+ plot_annotation(title = sj_type2[mm]) & theme(plot.title = element_text(size = 27,hjust=0.5))
  p4_3<-p2
}


####### acid#######
for(mm in 4){
  data<-read.table(paste0("new_",sj_type[mm],"_all_antismash_bigscape_region_length_morethan5kb_addNMDC.txt"),header=T,sep="\t",quote='')
  data$Type[which(data$Type=="New")]<-"Unannotated"
  gtdb<-gsub(";","",gsub("p__","",str_extract(data$GTDB,"p__(.*?);")))
  data2<-cbind(data,gtdb)
  if(length(unique(data2$gtdb))>=10){
    needlist<-names(sort(table(data2$gtdb),decreasing = T))[1:10]
  }else{
    needlist<-names(sort(table(data2$gtdb),decreasing = T))
  }
  data3<-data2[which(!is.na(match(data2$gtdb,needlist))),]
  data3$BiG.SCAPE_class<-as.factor(data3$BiG.SCAPE_class)

  data_plot<-aggregate(data3$Region_filename,by=list(data3$gtdb,data3$BiG.SCAPE_class),length)
  colnames(data_plot)<-c("gtdb","BiG.SCAPE_class","region_number")
  bgclist<-rev(c("NRPS","PKSI","PKS-NRP_Hybrids","PKSother","RiPPs","Saccharides","Terpene","Others"))

  data4<-unique(data3[,c("Region_filename","gtdb")])
  gtdbnumber2<-aggregate(data4$Region_filename,by=list(data4$gtdb),length)
  gtdbnumber3<-gtdbnumber2[order(as.numeric(gtdbnumber2[,2]),decreasing = T),]
  colnames(gtdbnumber3)<-c("GTDB","Number_of_MAGs")
  gtdbnumber4<-matrix(nrow=nrow(gtdbnumber3)*2,ncol=3)
  gtdbnumber4[,1]<-rep(gtdbnumber3$GTDB,each=2)
  gtdbnumber4[,2]<-rep(c("Annotated","Unannotated"),times=nrow(gtdbnumber3))
  colnames(gtdbnumber4)<-c("Phyla","Type","Number_of_regions")
  for(i in 1:nrow(gtdbnumber4)){
    index<-which(data3$Type==gtdbnumber4[i,2] & data3$gtdb==gtdbnumber4[i,1])
    if(length(index)>0){
      gtdbnumber4[i,3]<-length(unique(data3[index,"Region_filename"]))
    }else{
      gtdbnumber4[i,3]<-0
    }
  }
  gtdbnumber4<-as.data.frame(gtdbnumber4)
  gtdbnumber4$Phyla<-factor(gtdbnumber4$Phyla,levels=rev(gtdbnumber3$GTDB))
  gtdbnumber4$Type<-factor(gtdbnumber4$Type,levels = rev(c("Annotated","Unannotated")))
  gtdbnumber4$Number_of_regions<-as.numeric(gtdbnumber4$Number_of_regions)
  
  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(gtdbnumber4)),ncol=1),gtdbnumber4)
  colnames(all_all0)<-colnames(all_all)
  all_all<-rbind(all_all,all_all0)

  for (i in seq(1, 2)) {
    for(j in 1:length(unique(gtdbnumber4[,1]))){
      if (i == 1) {
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = 1
      }else{
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = sum(gtdbnumber4$Number_of_regions[(2*(j-1)+1):(2*(j-1)+i-1)]) + gtdbnumber4$Number_of_regions[2*(j-1)+i]/4
      }
    }
  }
  p1<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    theme(axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))

  data_plot2<-matrix(nrow=length(needlist)*length(bgclist),ncol=ncol(data_plot))
  colnames(data_plot2)<-colnames(data_plot)
  data_plot2[,1]<-rep(as.character(gtdbnumber3$GTDB),times=length(bgclist))
  data_plot2[,2]<-rep(bgclist,each=length(needlist))
  for(i in 1:nrow(data_plot2)){
    index<-which(data_plot$gtdb==data_plot2[i,1] & data_plot$BiG.SCAPE_class==data_plot2[i,2])
    if(length(index)>0){
      data_plot2[i,3]<-data_plot[index,3]
    }else{
      data_plot2[i,3]<-0
    }
  }
  data_plot2<-as.data.frame(data_plot2)
  data_plot2$region_number<-as.numeric(data_plot2$region_number)

  gtdbnumber<-aggregate(data_plot2$region_number,by=list(data_plot2$gtdb),sum)
  for(i in 1:nrow(data_plot2)){
    data_plot2[i,1]<-paste0(data_plot2[i,1]," (",gtdbnumber[which(gtdbnumber[,1]==data_plot2[i,1]),2],")")
  }

  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(data_plot2)),ncol=1),data_plot2)
  colnames(all_all0)<-colnames(all_all2)
  all_all2<-rbind(all_all2,all_all0)

  data_plot2$gtdb<-factor(data_plot2$gtdb,levels = data_plot2$gtdb[length(unique(data_plot2$gtdb)):1])
  p2<-ggplot(data_plot2,aes(x=gtdb,y=region_number,fill=factor(BiG.SCAPE_class,levels=bgclist)))+geom_col(position="fill",width=1,color="black")+
    labs(x="Phyla",y="BGC content(%)",fill="Group")+
    theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=25,color="black"),legend.text=element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_y_continuous(breaks=seq(0,1,0.25),labels=c('0','25','50','75','100'),expand=c(0,0))+
    coord_flip()+
    scale_fill_manual(values=c("NRPS"="#377EB8","PKSI"="#4DAF4A","PKS-NRP_Hybrids"="#984EA3","PKSother"="#FF7F00","RiPPs"="#E41A1C","Saccharides"="#A65628","Terpene"="#F781BF","Others"="#999999"))

  p1_2<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    xlab(NULL)+theme(axis.line=element_line(color="black"),axis.line.y=element_blank(),axis.ticks.y = element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))
  p3_4<-p2+p1_2+ plot_layout(guides = 'collect')+ plot_annotation(title = sj_type2[mm]) & theme(plot.title = element_text(size = 27,hjust=0.5))
  p4_4<-p2
}


####### cold #######
for(mm in 5){
  data<-read.table(paste0("new_",sj_type[mm],"_all_antismash_bigscape_region_length_morethan5kb_addNMDC.txt"),header=T,sep="\t",quote='')
  data$Type[which(data$Type=="New")]<-"Unannotated"
  gtdb<-gsub(";","",gsub("p__","",str_extract(data$GTDB,"p__(.*?);")))
  data2<-cbind(data,gtdb)
  if(length(unique(data2$gtdb))>=10){
    needlist<-names(sort(table(data2$gtdb),decreasing = T))[1:10]
  }else{
    needlist<-names(sort(table(data2$gtdb),decreasing = T))
  }
  data3<-data2[which(!is.na(match(data2$gtdb,needlist))),]
  data3$BiG.SCAPE_class<-as.factor(data3$BiG.SCAPE_class)

  data_plot<-aggregate(data3$Region_filename,by=list(data3$gtdb,data3$BiG.SCAPE_class),length)
  colnames(data_plot)<-c("gtdb","BiG.SCAPE_class","region_number")
  bgclist<-rev(c("NRPS","PKSI","PKS-NRP_Hybrids","PKSother","RiPPs","Saccharides","Terpene","Others"))

  data4<-unique(data3[,c("Region_filename","gtdb")])
  gtdbnumber2<-aggregate(data4$Region_filename,by=list(data4$gtdb),length)
  gtdbnumber3<-gtdbnumber2[order(as.numeric(gtdbnumber2[,2]),decreasing = T),]
  colnames(gtdbnumber3)<-c("GTDB","Number_of_MAGs")
  gtdbnumber4<-matrix(nrow=nrow(gtdbnumber3)*2,ncol=3)
  gtdbnumber4[,1]<-rep(gtdbnumber3$GTDB,each=2)
  gtdbnumber4[,2]<-rep(c("Annotated","Unannotated"),times=nrow(gtdbnumber3))
  colnames(gtdbnumber4)<-c("Phyla","Type","Number_of_regions")
  for(i in 1:nrow(gtdbnumber4)){
    index<-which(data3$Type==gtdbnumber4[i,2] & data3$gtdb==gtdbnumber4[i,1])
    if(length(index)>0){
      gtdbnumber4[i,3]<-length(unique(data3[index,"Region_filename"]))
    }else{
      gtdbnumber4[i,3]<-0
    }
  }
  gtdbnumber4<-as.data.frame(gtdbnumber4)
  gtdbnumber4$Phyla<-factor(gtdbnumber4$Phyla,levels=rev(gtdbnumber3$GTDB))
  gtdbnumber4$Type<-factor(gtdbnumber4$Type,levels = rev(c("Annotated","Unannotated")))
  gtdbnumber4$Number_of_regions<-as.numeric(gtdbnumber4$Number_of_regions)

  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(gtdbnumber4)),ncol=1),gtdbnumber4)
  colnames(all_all0)<-colnames(all_all)
  all_all<-rbind(all_all,all_all0)

  for (i in seq(1, 2)) {
    for(j in 1:length(unique(gtdbnumber4[,1]))){
      if (i == 1) {
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = 1
      }else{
        gtdbnumber4$Number_of_regions2[2*(j-1)+i] = sum(gtdbnumber4$Number_of_regions[(2*(j-1)+1):(2*(j-1)+i-1)]) + gtdbnumber4$Number_of_regions[2*(j-1)+i]/4
      }
    }
  }
  p1<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    theme(axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))

  data_plot2<-matrix(nrow=length(needlist)*length(bgclist),ncol=ncol(data_plot))
  colnames(data_plot2)<-colnames(data_plot)
  data_plot2[,1]<-rep(as.character(gtdbnumber3$GTDB),times=length(bgclist))
  data_plot2[,2]<-rep(bgclist,each=length(needlist))
  for(i in 1:nrow(data_plot2)){
    index<-which(data_plot$gtdb==data_plot2[i,1] & data_plot$BiG.SCAPE_class==data_plot2[i,2])
    if(length(index)>0){
      data_plot2[i,3]<-data_plot[index,3]
    }else{
      data_plot2[i,3]<-0
    }
  }
  data_plot2<-as.data.frame(data_plot2)
  data_plot2$region_number<-as.numeric(data_plot2$region_number)

  gtdbnumber<-aggregate(data_plot2$region_number,by=list(data_plot2$gtdb),sum)
  for(i in 1:nrow(data_plot2)){
    data_plot2[i,1]<-paste0(data_plot2[i,1]," (",gtdbnumber[which(gtdbnumber[,1]==data_plot2[i,1]),2],")")
  }

  all_all0<-cbind(matrix(rep(sj_type2[mm],nrow(data_plot2)),ncol=1),data_plot2)
  colnames(all_all0)<-colnames(all_all2)
  all_all2<-rbind(all_all2,all_all0)

  data_plot2$gtdb<-factor(data_plot2$gtdb,levels = data_plot2$gtdb[length(unique(data_plot2$gtdb)):1])
  p2<-ggplot(data_plot2,aes(x=gtdb,y=region_number,fill=factor(BiG.SCAPE_class,levels=bgclist)))+geom_col(position="fill",width=1,color="black")+
    labs(x="Phyla",y="BGC content(%)",fill="Group")+
    theme(panel.background=element_blank(),panel.grid=element_blank(),axis.line=element_line(color="black"),axis.text.y=element_text(size=20,color="black"),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=25,color="black"),legend.text=element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_y_continuous(breaks=seq(0,1,0.25),labels=c('0','25','50','75','100'),expand=c(0,0))+
    coord_flip()+
    scale_fill_manual(values=c("NRPS"="#377EB8","PKSI"="#4DAF4A","PKS-NRP_Hybrids"="#984EA3","PKSother"="#FF7F00","RiPPs"="#E41A1C","Saccharides"="#A65628","Terpene"="#F781BF","Others"="#999999"))

  p1_2<-ggplot(gtdbnumber4,aes(x=Phyla,y=Number_of_regions,fill=Type))+geom_bar(stat="identity",width = 1,color="black")+coord_flip()+
    geom_text(aes(x=Phyla,y=as.numeric(Number_of_regions2),label = Number_of_regions),vjust=0.5,hjust=0, size=5, color="black")+theme_classic()+
    xlab(NULL)+theme(axis.line=element_line(color="black"),axis.line.y=element_blank(),axis.ticks.y = element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(size=20,color="black",vjust=0.5,hjust=0.5,angle=0),axis.title=element_text(size=20,color="black"),legend.text = element_text(size=20,color="black"),legend.title = element_text(size=20,color="black"))+
    scale_fill_manual(values = c("Annotated"="#F8766D","Unannotated"="#00BFC4"),breaks=c("Annotated","Unannotated"))
  p3_5<-p2+p1_2+ plot_layout(guides = 'collect')+ plot_annotation(title = sj_type2[mm]) & theme(plot.title = element_text(size = 27,hjust=0.5))
  p4_5<-p2
}


######## combine all habitats #############
setwd("D:/工作/metagenome/New bigscape")
all_all<-all_all[-1,]
write.table(all_all,"BGC_top10_phyla_protein_number_annotedORunannoted_all_habitats_20240409.txt",sep="\t",col.names = T,row.names = F,quote = F)
all_all2<-all_all2[-1,]
write.table(all_all2,"BGC_top10_phyla_protein_number_inBGCgroup_all_habitats_20240409.txt",sep="\t",col.names = T,row.names = F,quote = F)

pdf("BGC_top10_phyla_protein_number_inBGCgroup_all_habitats_annotedORunannoted_20240409.pdf",width=60,height=25)
((p3_4|p3_5)/(p3_1|p3_2)/(p3_3|plot_spacer()))+ plot_layout(guides = 'collect')+plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size=30))
dev.off()

pdf("BGC_top10_phyla_protein_number_inBGCgroup_each_habitat_barplot_fill_20231007.pdf",width=26,height=17)
((p4_4|p4_5)/(p4_1|p4_2)/(p4_3|plot_spacer()))+ plot_layout(guides = 'collect')+plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size=30))
dev.off()