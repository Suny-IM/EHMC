setwd("D:/metagenome")
library(UpSetR)
library(RColorBrewer)
data0<-read.table("Biome5_OTU_16963_rongyu_table.txt",sep="\t")
tempdata<-aggregate(data0$V1,by=list(data0$V2,data0$V3),length)
colnames(tempdata)<-c("Habitat","OTU_Cluster","MAGs_number")

sj_type<-unique(data0[,2])
otu<-unique(data0[,3])
data<-matrix(ncol=(length(sj_type)+1),nrow=length(otu))
colnames(data)<-c("OTU_Cluster",sj_type)

for(i in 1:nrow(data)){
  for(j in 2:ncol(data)){
    data[i,1]<-otu[i]
    index<-which(tempdata$OTU_Cluster==otu[i] & tempdata$Habitat==colnames(data)[j])
    if(length(index)>0){
      data[i,j]<-tempdata[index,3]
    }else{
      data[i,j]<-0
    }
  }
}
data2<-as.data.frame(data)
for(j in 2:ncol(data2)){
  data2[,j]<-as.numeric(data2[,j])
}
data3<-data2
for(j in 2:ncol(data3)){
  index<-which(data3[,j]>1)
  if(length(index)>0){
    data3[index,j]<-1
  }
}


pdf("OTU_UpSet.pdf",width = 11,height = 6)
upset(data3,
      nsets=ncol(data3),
      nintersects=100,
      sets=c("Hot","Acid","Cold","Saline-alkaline","Pressure"),
      number.angles = 0,
      point.size=4,
      line.size=1,
      mainbar.y.label="Intersection size",
      main.bar.color = 'black',
      matrix.color="black",
      sets.x.label="Nubmer of species",
      sets.bar.color=brewer.pal(5,"Set1"),
      mb.ratio = c(0.7, 0.3),
      order.by = "freq",
      text.scale=c(1.5,1.5,1.5,1.5,1.5,1.5),
      shade.color="red"
)
dev.off()
