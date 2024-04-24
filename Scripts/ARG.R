setwd("D:/metagenome/protein_analysis")
data<-read.table("Habitat_8_rgi_arotypenum_id60_drugclass.txt",sep="\t",header = T,quote = '',check.names = F)
data$info<-paste(data$Habitat_id,data$MAG,sep="--")
data2<-data[1:19887,]
mgedata<-read.table("Habitat_8_mge_rgi_arotypenum_id60_drugclass.txt",sep="\t",header = T,quote = '',check.names = F)
mgedata$info<-paste(mgedata$Habitat_id,mgedata$MAG,sep="--")
mgedata2<-mgedata[1:19887,]

qualityinfo<-read.table("D:/metagenome/extreme_5biomedrep_16963_add_path_rongyu_forMGE_20240403_addInfo.txt",sep="\t",header = T,quote='')
qualityinfo$info<-paste(qualityinfo$Habitat,qualityinfo$MAG,sep="--")
for(i in 1:nrow(qualityinfo)){
  if(length(grep("^C",qualityinfo[i,"MAG"]))>0){
    qualityinfo[i,"info"]<-paste(qualityinfo$Habitat[i],qualityinfo$NMDC_ID[i],sep="--")
  }
}
qualityinfo2<-qualityinfo[which(qualityinfo$Quality=="High-quality"),]
data3<-data2[which(!is.na(match(data2$info,qualityinfo2$info))),c(1:2,9:ncol(data2))]
mgedata3<-mgedata2[which(!is.na(match(mgedata2$info,qualityinfo2$info))),c(1:2,9:ncol(mgedata2))]
colnames(data3)<-gsub("RGI-","",colnames(data3))
colnames(mgedata3)<-gsub("MGE_RGI-","",colnames(mgedata3))


mgedata3_1<-cbind(mgedata3,matrix(0,nrow=nrow(mgedata3),ncol=length(setdiff(colnames(data3),colnames(mgedata3))),dimnames=list(rownames(mgedata3),setdiff(colnames(data3),colnames(mgedata3)))))
mgedata4<-mgedata3_1[,colnames(data3)]
rownames(data3)<-data3$info
rownames(mgedata4)<-mgedata4$info
data4<-data3[rownames(mgedata4),]

## Figure 5 e 
mgedata5<-mgedata4[,3:(ncol(mgedata4)-1)]
dff<-mgedata5
dff_temp<-apply(dff,1,sum)
dff2<-dff[-which(dff_temp==0),]


library(ComplexHeatmap)

annot_df2 = data.frame(Habitat_group=sapply(strsplit(rownames(dff2),split = "--"),"[",1))
annot_col2 = list(Habitat_group=c("Hot"="#E41A1C","Acid"="#377EB8","Cold"="#4DAF4A","Saline-alkaline"="#984EA3","Pressure"="#FF7F00"))
ha2 <- rowAnnotation(df = annot_df2, col = annot_col2, show_legend = T,annotation_name_rot = 45,
                     annotation_legend_param = list(
                       Habitat_group = list(ncol=1,at=c("Acid","Cold","Hot","Pressure","Saline-alkaline"),labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20))))

library(circlize)
mycol <- colorRamp2(c(0,1,2,3), c("white","#F8766D","#BD0026","#00BFC4"))

pdf(file="D:/metagenome/High_Quality_MAGs_mgeAsso_arotypenumber_in_antibiotic_heatmap.pdf",width=20,height=10)
ht_list<-Heatmap(dff2, name="Types of ARGs",col=mycol,
                 left_annotation = ha2, row_title = "High Quality MAGs",column_title = "Antibiotic",
                 column_names_rot = 45,
                 column_title_side ="bottom",
                 row_title_side ="left",
                 row_names_side = "right",
                 rect_gp = gpar(col = "black", lwd = 1),
                 heatmap_legend_param = list(direction = "vertical",labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20)),
                 column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 row_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 show_row_names = F, show_column_names = T,
                 cluster_row_slices = FALSE,
                 cluster_column_slices = FALSE,
                 cluster_rows = T, cluster_columns = F,
                 row_split = annot_df2,column_split = annot_df,
                 show_row_dend = F,
)
draw(ht_list, heatmap_legend_side = "right",annotation_legend_side = "right")
dev.off()

## Supplementary Figure 6
data5<-data4[,3:(ncol(data4)-1)]
dff<-data5
dff_temp<-apply(dff,1,sum)
dff2<-dff[-which(dff_temp==0),]

library(ComplexHeatmap)
annot_df2 = data.frame(Habitat_group=sapply(strsplit(rownames(dff2),split = "--"),"[",1))
annot_col2 = list(Habitat_group=c("Hot"="#E41A1C","Acid"="#377EB8","Cold"="#4DAF4A","Saline-alkaline"="#984EA3","Pressure"="#FF7F00"))
ha2 <- rowAnnotation(df = annot_df2, col = annot_col2, show_legend = T,annotation_name_rot = 45,
                     annotation_legend_param = list(
                       Habitat_group = list(ncol=1,at=c("Acid","Cold","Hot","Pressure","Saline-alkaline"),labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20))))

library(circlize)
mycol <- colorRamp2(c(0,1,2,12), c("white","#F8766D","#BD0026","#00BFC4"))

pdf(file="D:/metagenome/High_Quality_MAGs_arotypenumber_in_antibiotic_heatmap.pdf",width=20,height=25)
ht_list<-Heatmap(dff2, name="Number",col=mycol,
                 left_annotation = ha2, row_title = "MAGs", column_title = "Antibiotic",
                 column_names_rot = 45,
                 column_title_side ="bottom",
                 row_title_side ="left",
                 row_names_side = "right",
                 heatmap_legend_param = list(direction = "vertical",labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20)),
                 column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 row_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 show_row_names = F, show_column_names = T,
                 cluster_row_slices = FALSE,
                 cluster_column_slices = FALSE,
                 cluster_rows = T, cluster_columns = F,
                 row_split = annot_df2,column_split = annot_df,
                 show_row_dend = F
)
draw(ht_list, heatmap_legend_side = "right",annotation_legend_side = "right")
dev.off()
