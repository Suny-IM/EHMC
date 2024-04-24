setwd("D:/metagenome/New bigscape")
sj_type<-c("Hot","Pressure","Saline_alkaline","Acid","Cold")
sj_type2<-c("Hot","Pressure","Saline_alkaline","Acid","Cold")
all_all<-matrix(ncol=4,nrow=1)
colnames(all_all)<-c("Habitat","Phyla","Type","Number_of_proteins")
all_all2<-matrix(ncol=4,nrow=1)
colnames(all_all2)<-c("Habitat","Phyla","BiG.SCAPE_class","protein_number")

for(mm in 1:length(sj_type)){
  data<-read.table(paste0("new_",sj_type[mm],"_all_antismash_bigscape_region_length_morethan5kb_addNMDC.txt"),header=T,sep="\t",quote='')
  data$Type[which(data$Type=="New")]<-"Unannotated"
  gtdb<-gsub(";","",gsub("p__","",str_extract(data$GTDB,"p__(.*?);")))
  data2<-cbind(data,gtdb)

  needlist<-unique(data2$gtdb)
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
}
all_all2<-all_all2[-1,]
all_all2$Phyla<-sapply(strsplit(all_all2$Phyla,split="\\("),"[",1)

setwd("D:/metagenome/New bigscape")
bgclist<-c("NRPS","PKSI","PKS-NRP_Hybrids","PKSother","RiPPs","Saccharides","Terpene","Others")
all_all3<-read.table("All_habitats_BGC_top10Phyla_BGCgroup_proteinNum_20240409.txt",sep="\t",header = T,quote = '')
all_all3$Phyla<-sapply(strsplit(all_all3$Phyla,split="\\("),"[",1)
phyla<-unique(all_all3$Phyla)

plot_data<-matrix(nrow=length(phyla)*length(unique(all_all3$Habitat))*length(bgclist),ncol=4)
plot_data[,1]<-rep(unique(all_all3$Habitat),each=length(phyla)*length(bgclist))
plot_data[,2]<-rep(rep(unique(all_all3$Phyla),each=length(bgclist)),times=length(unique(all_all3$Habitat)))
plot_data[,3]<-rep(rep(bgclist,times=length(unique(all_all3$Phyla))),times=length(unique(all_all3$Habitat)))
for(i in 1:nrow(plot_data)){
  index<-which(all_all2$Habitat==plot_data[i,1] & all_all2$Phyla==plot_data[i,2] & all_all2$BiG.SCAPE_class==plot_data[i,3])
  if(length(index)==1){
    plot_data[i,4]<-all_all2[index,4]
  }else if(length(index)==0){
    plot_data[i,4]<-0
  }else{
    print ("ERROR")
    print (i)
  }
}
plot_data<-as.data.frame(plot_data)
colnames(plot_data)<-colnames(all_all2)
plot_data$BiG.SCAPE_class<-factor(plot_data$BiG.SCAPE_class,levels=bgclist)
plot_data$Phyla<-factor(plot_data$Phyla,levels = sort(unique(plot_data$Phyla),decreasing = T))
plot_data$protein_number<-as.numeric(plot_data$protein_number)
plot_data2<-plot_data[-which(plot_data$protein_number==0),]
write.table(plot_data2,"Number_of_BGC_proteins_in_top10_BGC_phyla_bubble_chart_20240409.txt",col.names = T,row.names = F,sep = "\t")
ggplot(data=plot_data2,aes(x=BiG.SCAPE_class,y=Phyla,fill=BiG.SCAPE_class))+geom_point(aes(size=protein_number),shape=21, colour="black", stroke = 0.5)+scale_size_continuous(limits = c(1,max(plot_data2$protein_number)))+
  facet_grid(. ~ Habitat)+theme_classic()+
  theme(axis.line=element_line(color="black"),axis.ticks = element_line(color="black"),axis.text.x=element_text(size=10,color="black",vjust=1,hjust=1,angle=45),axis.text.y=element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"),legend.text = element_text(size=10,color="black"),legend.title = element_text(size=10,color="black"))+
  scale_color_manual(values=c("NRPS"="#377EB8","PKSI"="#4DAF4A","PKS-NRP_Hybrids"="#984EA3","PKSother"="#FF7F00","RiPPs"="#E41A1C","Saccharides"="#A65628","Terpene"="#F781BF","Others"="#999999"))
ggsave("Number_of_BGC_proteins_in_top10_BGC_phyla_bubble_chart_20240409.pdf",width=12,height = 6)