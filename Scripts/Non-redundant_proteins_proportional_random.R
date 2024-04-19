## hot
setwd("D:/metagenome/protein_analysis/hot")
data<-read.table("MMseqs2_clusterRes_allfaa_statis.txt",sep="\t",header=T)
data$Method[which(data$Method=="single_merge")]<-"single_co"
data<-cbind(data,matrix(paste(data$Method,data$Identity),ncol=1))
colnames(data)[ncol(data)]<-"group"
data[,ncol(data)]<-paste0(data[,ncol(data)],"%identity")
data$group<-as.factor(data$group)
class(data$Selectfaa_cluster_num)

data0<-data
data<-data[which(data$Method=="single_co"),]
data$group<-gsub("single_co ","",data$group)

datause<-data[which(data$Proportion==100),]
xnum<-rep("98",nrow(datause))
datause1<-cbind(datause,xnum)
datause1$xnum<-as.numeric(datause1$xnum)
datause2<-cbind(datause1,datause1$Selectfaa_cluster_num)
colnames(datause2)[ncol(datause2)]<-"ynum"
datause2[which(datause2$Method=="single"),"ynum"]<-datause2[which(datause2$Method=="single"),"ynum"]-200000
datause2[which(datause2$Method=="single_co"),"ynum"]<-datause2[which(datause2$Method=="single_co"),"ynum"]+200000

library(ggplot2)
p1<-ggplot(data,aes(x=Proportion,y=Selectfaa_cluster_num,color=group))+geom_point()+geom_line(stat = "identity")+theme_bw()+
  labs(title="Hot",x="Proportion of proteins subsampled (%)",y="Number of proteins clusters")+
  guides(color=guide_legend(title = "Group"))+
  theme(axis.text=element_text(color="black",size=12),axis.title = element_text(size=15,color="black"),plot.title = element_text(size=18,color="black",hjust=0.5),legend.text = element_text(size=12),legend.title = element_text(size=15))+
  annotate("text",x=datause2$xnum,y=datause2$ynum,label=datause2$Selectfaa_cluster_num,size=4)

## pressure
setwd("D:/metagenome/protein_analysis/pressure")
data<-read.table("MMseqs2_clusterRes_allfaa_statis.txt",sep="\t",header=T)
data$Method[which(data$Method=="single_merge")]<-"single_co"
data<-cbind(data,matrix(paste(data$Method,data$Identity),ncol=1))
colnames(data)[ncol(data)]<-"group"
data[,ncol(data)]<-paste0(data[,ncol(data)],"%identity")
data$group<-as.factor(data$group)
class(data$Selectfaa_cluster_num)

data0<-data
data<-data[which(data$Method=="single_co"),]
data$group<-gsub("single_co ","",data$group)

datause<-data[which(data$Proportion==100),]
xnum<-rep("98",nrow(datause))
datause1<-cbind(datause,xnum)
datause1$xnum<-as.numeric(datause1$xnum)
datause2<-cbind(datause1,datause1$Selectfaa_cluster_num)
colnames(datause2)[ncol(datause2)]<-"ynum"
datause2[which(datause2$Method=="single"),"ynum"]<-datause2[which(datause2$Method=="single"),"ynum"]-200000
datause2[which(datause2$Method=="single_co"),"ynum"]<-datause2[which(datause2$Method=="single_co"),"ynum"]+200000

library(ggplot2)
p2<-ggplot(data,aes(x=Proportion,y=Selectfaa_cluster_num,color=group))+geom_point()+geom_line(stat = "identity")+theme_bw()+
  labs(title="Pressure",x="Proportion of proteins subsampled (%)",y="Number of proteins clusters")+
  guides(color=guide_legend(title = "Group"))+
  theme(axis.text=element_text(color="black",size=12),axis.title = element_text(size=15,color="black"),plot.title = element_text(size=18,color="black",hjust=0.5),legend.text = element_text(size=12),legend.title = element_text(size=15))+
  annotate("text",x=datause2$xnum,y=datause2$ynum,label=datause2$Selectfaa_cluster_num,size=4)

## saline_alkaline
setwd("D:/metagenome/protein_analysis/saline_alkaline")
data<-read.table("MMseqs2_clusterRes_allfaa_statis.txt",sep="\t",header=T)
data$Method[which(data$Method=="single_merge")]<-"single_co"
data<-cbind(data,matrix(paste(data$Method,data$Identity),ncol=1))
colnames(data)[ncol(data)]<-"group"
data[,ncol(data)]<-paste0(data[,ncol(data)],"%identity")
data$group<-as.factor(data$group)
class(data$Selectfaa_cluster_num)

data0<-data
data<-data[which(data$Method=="single_co"),]
data$group<-gsub("single_co ","",data$group)

datause<-data[which(data$Proportion==100),]
xnum<-rep("98",nrow(datause))
datause1<-cbind(datause,xnum)
datause1$xnum<-as.numeric(datause1$xnum)
datause2<-cbind(datause1,datause1$Selectfaa_cluster_num)
colnames(datause2)[ncol(datause2)]<-"ynum"
datause2[which(datause2$Method=="single"),"ynum"]<-datause2[which(datause2$Method=="single"),"ynum"]-200000
datause2[which(datause2$Method=="single_co"),"ynum"]<-datause2[which(datause2$Method=="single_co"),"ynum"]+200000

library(ggplot2)
p3<-ggplot(data,aes(x=Proportion,y=Selectfaa_cluster_num,color=group))+geom_point()+geom_line(stat = "identity")+theme_bw()+
  labs(title="Saline_alkaline",x="Proportion of proteins subsampled (%)",y="Number of proteins clusters")+
  guides(color=guide_legend(title = "Group"))+
  theme(axis.text=element_text(color="black",size=12),axis.title = element_text(size=15,color="black"),plot.title = element_text(size=18,color="black",hjust=0.5),legend.text = element_text(size=12),legend.title = element_text(size=15))+
  annotate("text",x=datause2$xnum,y=datause2$ynum,label=datause2$Selectfaa_cluster_num,size=4)

## cold
setwd("D:/metagenome/protein_analysis/cold")
data<-read.table("MMseqs2_clusterRes_allfaa_statis.txt",sep="\t",header=T)
data$Method[which(data$Method=="single_merge")]<-"single_co"
data<-cbind(data,matrix(paste(data$Method,data$Identity),ncol=1))
colnames(data)[ncol(data)]<-"group"
data[,ncol(data)]<-paste0(data[,ncol(data)],"%identity")
data$group<-as.factor(data$group)
class(data$Selectfaa_cluster_num)

data0<-data
data<-data[which(data$Method=="single_co"),]
data$group<-gsub("single_co ","",data$group)

datause<-data[which(data$Proportion==100),]
xnum<-rep("98",nrow(datause))
datause1<-cbind(datause,xnum)
datause1$xnum<-as.numeric(datause1$xnum)
datause2<-cbind(datause1,datause1$Selectfaa_cluster_num)
colnames(datause2)[ncol(datause2)]<-"ynum"
datause2[which(datause2$Method=="single"),"ynum"]<-datause2[which(datause2$Method=="single"),"ynum"]-200000
datause2[which(datause2$Method=="single_co"),"ynum"]<-datause2[which(datause2$Method=="single_co"),"ynum"]+200000

library(ggplot2)
p4<-ggplot(data,aes(x=Proportion,y=Selectfaa_cluster_num,color=group))+geom_point()+geom_line(stat = "identity")+theme_bw()+
  labs(title="Cold",x="Proportion of proteins subsampled (%)",y="Number of proteins clusters")+
  guides(color=guide_legend(title = "Group"))+
  theme(axis.text=element_text(color="black",size=12),axis.title = element_text(size=15,color="black"),plot.title = element_text(size=18,color="black",hjust=0.5),legend.text = element_text(size=12),legend.title = element_text(size=15))+
  annotate("text",x=datause2$xnum,y=datause2$ynum,label=datause2$Selectfaa_cluster_num,size=4)

## acid
setwd("D:/metagenome/protein_analysis/acid")
data<-read.table("MMseqs2_clusterRes_allfaa_statis.txt",sep="\t",header=T)
data$Method[which(data$Method=="single_merge")]<-"single_co"
data<-cbind(data,matrix(paste(data$Method,data$Identity),ncol=1))
colnames(data)[ncol(data)]<-"group"
data[,ncol(data)]<-paste0(data[,ncol(data)],"%identity")
data$group<-as.factor(data$group)
class(data$Selectfaa_cluster_num)

data0<-data
data<-data[which(data$Method=="single_co"),]
data$group<-gsub("single_co ","",data$group)

datause<-data[which(data$Proportion==100),]
xnum<-rep("98",nrow(datause))
datause1<-cbind(datause,xnum)
datause1$xnum<-as.numeric(datause1$xnum)
datause2<-cbind(datause1,datause1$Selectfaa_cluster_num)
colnames(datause2)[ncol(datause2)]<-"ynum"
datause2[which(datause2$Method=="single"),"ynum"]<-datause2[which(datause2$Method=="single"),"ynum"]-100000
datause2[which(datause2$Method=="single_co"),"ynum"]<-datause2[which(datause2$Method=="single_co"),"ynum"]+100000

library(ggplot2)
p5<-ggplot(data,aes(x=Proportion,y=Selectfaa_cluster_num,color=group))+geom_point()+geom_line(stat = "identity")+theme_bw()+
  labs(title="Acid",x="Proportion of proteins subsampled (%)",y="Number of proteins clusters")+
  guides(color=guide_legend(title = "Group"))+
  theme(axis.text=element_text(color="black",size=12),axis.title = element_text(size=15,color="black"),plot.title = element_text(size=18,color="black",hjust=0.5),legend.text = element_text(size=12),legend.title = element_text(size=15))+
  annotate("text",x=datause2$xnum,y=datause2$ynum,label=datause2$Selectfaa_cluster_num,size=4)


##### combine all habitat
setwd("D:/metagenome/protein_analysis/")
library(patchwork)
pdf("Non-redundant_proteins_in_proportional_random_all_habitats_line_chart.pdf",width = 16,height = 16)
p5+p4+p1+p2+p3+plot_layout(ncol = 2)  + guide_area() + plot_layout(guides = 'collect')
dev.off()
