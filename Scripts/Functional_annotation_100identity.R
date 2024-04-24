###################################### Supplementary Figure4 a ######################################################
library(ggplot2)
setwd("D:/metagenome/protein_analysis/")
sj_type<-c("hot","pressure","saline_alkaline","acid","cold")
data<-read.table(paste0(sj_type[1],"/clusterRes_allfaa_allMAGs_100identity_represent_function_id40_cov80_result_anno_statis.txt"),sep="\t",header=T)
for(i in 2:length(sj_type)){
  data2<-read.table(paste0(sj_type[i],"/clusterRes_allfaa_allMAGs_100identity_represent_function_id40_cov80_result_anno_statis.txt"),sep="\t",header=T)
  data<-rbind(data,data2)
}
data2<-matrix(nrow=nrow(data)*2,ncol=4)
colnames(data2)<-c("Habitat","Function","Group","Protein_number")
data2[,1]<-rep(data[,1],each=2)
data2[,2]<-rep(data[,2],each=2)
data2[,3]<-rep(c("Annotated","Unannotated"),times=nrow(data))
for(i in 1:nrow(data2)){
  index<-which(data$Habitat==data2[i,1] & data$Function==data2[i,2])
  if(data2[i,3]=="Annotated"){
    data2[i,4]<-data[index,4]
  }else{
    data2[i,4]<-data[index,5]
  }
}
data2<-as.data.frame(data2)
data2$Protein_number<-as.numeric(data2$Protein_number)

temp0<-read.table("all5habitat/clusterRes_allfaa_allMAGs_100identity_represent_function_id40_cov80_result_anno_statis.txt",header = T,sep="\t")
temp0_1<-rbind(cbind(rep("Annotated",5),temp0[,"Annoted_num"]),cbind(rep("Unannotated",5),temp0[,"Unannoted_num"]))
temp0_2<-rbind(temp0[,1:2],temp0[,1:2])
temp<-cbind(temp0_2,temp0_1)
colnames(temp)<-colnames(data2)
temp$Habitat<-"All"

library(Hmisc)
data3<-rbind(data2,temp)
data3<-as.data.frame(data3)
data3$Habitat<-capitalize(data3$Habitat)
data3$Habitat<-factor(data3$Habitat,levels = c("All","Acid","Cold","Hot","Pressure","Saline_alkaline"))
data3$Protein_number<-as.numeric(data3$Protein_number)
data3$Function<-as.factor(data3$Function)
data3$Group<-as.factor(data3$Group)

data3_0<-data3
data3$Group<-gsub("Annotated","No",data3$Group)
data3$Group<-gsub("Unannotated","Yes",data3$Group)

ggplot(data3,aes(x=Function,y=Protein_number,fill=Group))+geom_bar(stat="identity")+facet_grid(.~Habitat)+
  theme_classic()+labs(fill = "Unannotated Protein")+
  theme(axis.line=element_line(color="black"),axis.text.y=element_text(size=10,color="black"),axis.text.x=element_text(size=10,color="black",vjust=1,hjust=1,angle=45),axis.title=element_text(size=12,color="black"),legend.text=element_text(size=10,color="black"),legend.title=element_text(size=10,color="black"))
ggsave(filename = "Number_of_protein_annotations_in_each_database_barplot_all_habitats_id40cov80_20240412.pdf",width=12,height=6.75)
write.table(data3_0,"Number_of_protein_annotations_in_each_database_barplot_all_habitats_id40cov80_20240412.txt",col.names = T,row.names = F,sep="\t",quote = F)



###################################### Supplementary Figure4 b ######################################################
setwd("D:/metagenome/protein_analysis/all5habitat")
data<-read.table("clusterRes_allfaa_allMAGs_100identity_represent_function_id40_cov80_result_anno_statis.txt",sep="\t",header=T,quote='')
nohitnum<-as.numeric(data[which(data$Function=="COG"),"Unannoted_num"])
coginfo00<-read.table("D:/metagenome/cog_classify_matched.list",sep="\t")
coginfo01<-coginfo00[order(coginfo00[,1]),]
coginfo<-coginfo01[order(coginfo01[,3]),]
allnum_annoprotein4<-read.table("function_id40_cov80/COG_100identity_protein_number_all_habitats_20240412.txt",sep="\t",header = T,quote='')
data_statis<-aggregate(allnum_annoprotein4$Annoprotein_number,by=list(allnum_annoprotein4$COG_functional_category),sum)
data_2_statis<-merge(data_statis,coginfo,by.x="Group.1",by.y="V1")
colnames(data_2_statis)<-c("COG_functional_category","Number_of_protein","COG_functional_category_description","COG_function_group")
data_2_statis2<-data_2_statis[order(data_2_statis$COG_function_group),]
a<-matrix(c("No hit",nohitnum,"No hit","No hit"),nrow=1)
colnames(a)<-colnames(data_2_statis2)
data_2_statis3<-rbind(data_2_statis2,a)
data_2_statis3$Number_of_protein<-as.numeric(data_2_statis3$Number_of_protein)
data_2_statis3$COG_functional_category_description<-factor(data_2_statis3$COG_functional_category_description,levels = rev(data_2_statis3$COG_functional_category_description))
p<-ggplot(data_2_statis3,aes(x=1,y=Number_of_protein,fill=COG_functional_category_description))+
  geom_bar(position = position_stack(), stat ="identity",colour="black") +
  labs(x='Representative proteins',y='Number of proteins')+
  scale_fill_manual(values = c("Cell cycle control, cell division, chromosome partitioning"="#FF0000","Cell wall/membrane/envelope biogenesis"= "#FF2400","Cell motility"= "#FF4900","Posttranslational modification, protein turnover, chaperones"= "#FF6D00",
                               "Signal transduction mechanisms"="#FF9200","Intracellular trafficking, secretion, and vesicular transport"= "#FFB600","Defense mechanisms"= "#FFDB00","Extracellular structures"= "#FFFF00","Nuclear structure"= "#FFFF40","Cytoskeleton"= "#FFFFBF",
                               "RNA processing and modification"= "#00441B","Chromatin structure and dynamics"= "#006D2C","Translation, ribosomal structure and biogenesis"= "#238B45","Transcription"= "#41AE76","Replication, recombination and repair"= "#66C2A4",
                               "Energy production and conversion"= "#08306B","Amino acid transport and metabolism"= "#08519C","Nucleotide transport and metabolism"= "#2171B5","Carbohydrate transport and metabolism"= "#4292C6","Coenzyme transport and metabolism"= "#6BAED6",
                               "Lipid transport and metabolism"="#9ECAE1","Inorganic ion transport and metabolism"= "#C6DBEF","Secondary metabolites biosynthesis, transport and catabolism"= "#DEEBF7","Mobilome: prophages, transposons"= "#F7FBFF","General function prediction only"= "grey70","Function unknown"= "grey85","No hit"="grey95"))+
  guides(fill = guide_legend(title = 'COG functional category',reverse=TRUE))+
  coord_flip()+theme_test()+
  theme(legend.position="bottom",legend.text = element_text(size=14),legend.title = element_text(size=25),axis.text=element_text(color="black",size=20),axis.title = element_text(size=25,color="black"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        strip.background = element_rect(colour = "white", fill = "grey"))+
  theme(legend.title=element_blank())
ggsave(p,filename = "COG_100identity_protein_number_all_habitats_id40_cov80_barplot_20240412.pdf",width=26,height=10)
