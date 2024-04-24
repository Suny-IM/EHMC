setwd("D:/metagenome")
result4<-read.table("High_Quality_MAGs_Phylum_core_dispensable_unique_moduleNum.txt",sep="\t",header=T,check.names = F)
result4$info<-paste(result4$Habitat,result4$MAGs,sep="--")
rownames(result4)<-result4$info
result4_2<-result4[,3:4]

library(ComplexHeatmap)
annot_df2 = data.frame(Habitat_group=result4$Habitat)
annot_col2 = list(Habitat_group=c("Hot"="#E41A1C","Acid"="#377EB8","Cold"="#4DAF4A","Saline-alkaline"="#984EA3","Pressure"="#FF7F00"))
ha2 <- rowAnnotation(df = annot_df2, col = annot_col2, show_legend = T,annotation_name_rot = 45,
                     annotation_name_gp = gpar(fontsize=15),
                     annotation_legend_param = list(
                       Habitat_group = list(ncol=1,at=c("Acid","Cold","Hot","Pressure","Saline-alkaline"),labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20))))

library(circlize)
mycol <- colorRamp2(c(0, 10, 20,30,40,50,60,70,80,90,100,150,200), c("#221557","#313695","#323b81","#036db9","#2ca8e2","#1ECFE9","#a2d9f2","#fee5a7","#E5B42C","#fdac55","#e85d20","#b52e13","#8e190b"))

pdf(file="High_Quality_MAGs_Phylum_core_dispensable_moduleNum_heatmap.pdf",width=9,height=12)
ht_list<-Heatmap(result4_2, name="Number",col=mycol,
                 left_annotation = ha2, row_title = "High Quality MAGs", column_title = "Type",
                 column_names_rot = 0,
                 column_title_side ="bottom",
                 row_title_side ="left",
                 row_names_side = "right",
                 heatmap_legend_param = list(direction = "vertical",labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20)),
                 column_names_gp = gpar(fontsize = 15),
                 column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 row_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 show_row_names = F, show_column_names = T,
                 cluster_row_slices = FALSE,
                 cluster_column_slices = FALSE,
                 cluster_rows = T, cluster_columns = F,
                 row_split = annot_df2,
                 show_row_dend = F
)
draw(ht_list, heatmap_legend_side = "left",annotation_legend_side = "left")
dev.off()
