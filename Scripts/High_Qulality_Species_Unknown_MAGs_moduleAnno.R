setwd("D:/metagenome")

#################################  Figure 6a #################################################
data<-read.table("Supplement_allfaa_KEGG_result_module_absence_use_highquality_speciesunknown.txt",sep="\t",quote = '',header=T,check.names = F)
level2list<-unique(sapply(strsplit(colnames(data)[5:ncol(data)],split="--"),"[",2))
data_plot<-matrix(nrow=nrow(data),ncol=(length(level2list)+4))
data_plot[,1]<-data[,1]
data_plot[,2]<-data[,2]
data_plot[,3]<-data[,3]
data_plot[,4]<-data[,4]
colnames(data_plot)<-c(colnames(data)[1:4],level2list)
for(i in 1:nrow(data_plot)){
  for(j in 5:ncol(data_plot)){
    index<-grep(paste0(colnames(data_plot)[j],"--"),colnames(data))
    index2<-which(data[i,index]>0)
    data_plot[i,j]<-length(index2)/length(index)
  }
}

needlevel0<-matrix(c(sapply(strsplit(colnames(data)[5:ncol(data)],split="--"),"[",1),sapply(strsplit(colnames(data)[5:ncol(data)],split="--"),"[",2)),ncol=2)
needlevel<-unique(needlevel0)
colnames(needlevel)<-c("First_level","Second_level")
needlevel<-as.data.frame(needlevel)
needlevel<-needlevel[order(needlevel[,1]),]

data_plot2<-data_plot[,5:ncol(data_plot)]
rownames(data_plot2)<-data_plot[,2]
data_plot2<-as.data.frame(data_plot2)
for(j in 1:ncol(data_plot2)){
  data_plot2[,j]<-as.numeric(data_plot2[,j])
}
data_plot3<-data_plot2[,needlevel[,2]]


library(ComplexHeatmap)
annot_df2 = data.frame(Habitat_group=data_plot[,1])
annot_col2 = list(Habitat_group=c("Hot"="#E41A1C","Acid"="#377EB8","Cold"="#4DAF4A","Saline-alkaline"="#984EA3","Pressure"="#FF7F00"))
ha2 <- rowAnnotation(df = annot_df2, col = annot_col2, show_legend = T,annotation_name_rot = 45,
                     annotation_name_gp = gpar(fontsize=15),
                     annotation_legend_param = list(
                       Habitat_group = list(ncol=1,at=c("Acid","Cold","Hot","Pressure","Saline-alkaline"),labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20))))
annot_df = data.frame(First_level=needlevel[,1],Second_level=needlevel[,2])
annot_col = list(First_level=c("Amino acid metabolism"="#7FC97F","Biosynthesis of other secondary metabolites"="#BEAED4","Biosynthesis of terpenoids and polyketides"="#FDC086","Carbohydrate metabolism"="#FFFF99",
                               "Energy metabolism"="#386CB0","Glycan metabolism"="#BF5B17",
                               "Lipid metabolism"="#666666","Metabolism of cofactors and vitamins"="#1B9E77",
                               "Nucleotide metabolism"="#7570B3","Xenobiotics biodegradation"="#E7298A","Heavy metal metabolism"="#66A61E","Oxidative stress"="#E6AB02"),
                 Second_level=c("Arginine and proline metabolism"="#7FC97F","Aromatic amino acid metabolism"="#BEAED4","Aromatics degradation"="#FDC086","Arsenic"="#CCCCCC","ATP synthesis"="#FFFF99",
                                "Biosynthesis of beta-lactams"="#386CB0","Biosynthesis of other antibiotics"="#F0027F","Biosynthesis of other bacterial compounds"="#BF5B17",
                                "Biosynthesis of other fungal compounds"="#666666","Branched-chain amino acid metabolism"="#1B9E77","Carbon fixation"="#D95F02","Catalase"="#E41A1C",
                                "Central carbohydrate metabolism"="#7570B3","Cofactor and vitamin metabolism"="#E7298A","Copper"="#377EB8","Cysteine and methionine metabolism"="#66A61E","Drug resistance"="#E6AB02",
                                "Enediyne biosynthesis"="#A6761D","Fatty acid metabolism"="#A6CEE3","Glycan biosynthesis"="#1F78B4","Glycosaminoglycan metabolism"="#B2DF8A",
                                "Histidine metabolism"="#33A02C","Lipid metabolism"="#FB9A99","Lipopolysaccharide metabolism"="#E31A1C","Lysine metabolism"="#FDBF6F",
                                "Macrolide biosynthesis"="#FF7F00","Mercury"="#4DAF4A","Metabolic capacity"="#CAB2D6","Methane metabolism"="#6A3D9A","Nitrogen metabolism"="#B15928",
                                "Other amino acid metabolism"="#FBB4AE","Other carbohydrate metabolism"="#B3CDE3","Pathogenicity"="#CCEBC5","Photosynthesis"="#DECBE4","Plant pathogenicity"="#FED9A6",
                                "Plant terpenoid biosynthesis"="#FFFFCC","Polyamine biosynthesis"="#E5D8BD","Polyketide sugar unit biosynthesis"="#FDDAEC","Purine metabolism"="#F2F2F2",
                                "Pyrimidine metabolism"="#B3E2CD","Selenium"="#984EA3","Serine and threonine metabolism"="#FDCDAC","Sterol biosynthesis"="#CBD5E8","Sulfur metabolism"="#F4CAE4","Superoxide dismutase"="#FFFF33","Symbiosis"="#E6F5C9",
                                "Terpenoid backbone biosynthesis"="#FFF2AE","Type II polyketide biosynthesis"="#F1E2CC"))

ha <- HeatmapAnnotation(df = annot_df, col = annot_col, #gp = gpar(col = "grey80"),
                        show_legend = T,
                        annotation_legend_param = list(First_level = list(ncol=1,at=unique(needlevel[,1]),labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20)),
                                                       Second_level=list(ncol=1,at=unique(needlevel[,2]),labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20))))

pdf(file="Each_habitat_high-quality_species-unknown_MAGs_module_percentage_heatmap.pdf",width=16,height=16)
ht_list<-Heatmap(data_plot3, name="Proportion",top_annotation = ha,
                 left_annotation = ha2, row_title = "MAGs", column_title = "Module",
                 column_names_rot = 45,
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
                 cluster_rows = T, cluster_columns = T,
                 row_split = annot_df2,
                 column_split = annot_df$First_level,
                 show_row_dend = F,show_column_dend = F,
                 column_dend_height = unit(2, "cm"), row_dend_width = unit(1, "cm")
)
draw(ht_list, heatmap_legend_side = "right",annotation_legend_side = "right")
dev.off()


#################################  Figure 6b #################################################
acidmags<-data_plot[which(data_plot[,1]=="Acid"),2]
data_plot3_acid<-data_plot3[acidmags,]
annot_df = data.frame(First_level=needlevel[,1],Second_level=needlevel[,2])
annot_col = list(First_level=c("Amino acid metabolism"="#7FC97F","Biosynthesis of other secondary metabolites"="#BEAED4","Biosynthesis of terpenoids and polyketides"="#FDC086","Carbohydrate metabolism"="#FFFF99",
                               "Energy metabolism"="#386CB0","Glycan metabolism"="#BF5B17",
                               "Lipid metabolism"="#666666","Metabolism of cofactors and vitamins"="#1B9E77",
                               "Nucleotide metabolism"="#7570B3","Xenobiotics biodegradation"="#E7298A","Heavy metal metabolism"="#66A61E","Oxidative stress"="#E6AB02"),
                 Second_level=c("Arginine and proline metabolism"="#7FC97F","Aromatic amino acid metabolism"="#BEAED4","Aromatics degradation"="#FDC086","Arsenic"="#CCCCCC","ATP synthesis"="#FFFF99",
                                "Biosynthesis of beta-lactams"="#386CB0","Biosynthesis of other antibiotics"="#F0027F","Biosynthesis of other bacterial compounds"="#BF5B17",
                                "Biosynthesis of other fungal compounds"="#666666","Branched-chain amino acid metabolism"="#1B9E77","Carbon fixation"="#D95F02","Catalase"="#E41A1C",
                                "Central carbohydrate metabolism"="#7570B3","Cofactor and vitamin metabolism"="#E7298A","Copper"="#377EB8","Cysteine and methionine metabolism"="#66A61E","Drug resistance"="#E6AB02",
                                "Enediyne biosynthesis"="#A6761D","Fatty acid metabolism"="#A6CEE3","Glycan biosynthesis"="#1F78B4","Glycosaminoglycan metabolism"="#B2DF8A",
                                "Histidine metabolism"="#33A02C","Lipid metabolism"="#FB9A99","Lipopolysaccharide metabolism"="#E31A1C","Lysine metabolism"="#FDBF6F",
                                "Macrolide biosynthesis"="#FF7F00","Mercury"="#4DAF4A","Metabolic capacity"="#CAB2D6","Methane metabolism"="#6A3D9A","Nitrogen metabolism"="#B15928",
                                "Other amino acid metabolism"="#FBB4AE","Other carbohydrate metabolism"="#B3CDE3","Pathogenicity"="#CCEBC5","Photosynthesis"="#DECBE4","Plant pathogenicity"="#FED9A6",
                                "Plant terpenoid biosynthesis"="#FFFFCC","Polyamine biosynthesis"="#E5D8BD","Polyketide sugar unit biosynthesis"="#FDDAEC","Purine metabolism"="#F2F2F2",
                                "Pyrimidine metabolism"="#B3E2CD","Selenium"="#984EA3","Serine and threonine metabolism"="#FDCDAC","Sterol biosynthesis"="#CBD5E8","Sulfur metabolism"="#F4CAE4","Superoxide dismutase"="#FFFF33","Symbiosis"="#E6F5C9",
                                "Terpenoid backbone biosynthesis"="#FFF2AE","Type II polyketide biosynthesis"="#F1E2CC"))

ha <- HeatmapAnnotation(df = annot_df, col = annot_col,
                        show_legend = T,
                        annotation_legend_param = list(First_level = list(nrow=3,at=unique(needlevel[,1]),labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20)),
                                                       Second_level=list(nrow=8,at=unique(needlevel[,2]),labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20))))

pdf(file="Acid_high-quality_species-unknown_MAGs_module_percentage_heatmap.pdf",width=12,height=26)
ht_list<-Heatmap(data_plot3_acid, name="Proportion",top_annotation = ha,
                 row_title = "MAGs", column_title = "Module",
                 column_names_rot = 45,
                 column_title_side ="bottom",
                 row_title_side ="left",
                 row_names_side = "left",
                 heatmap_legend_param = list(direction = "horizontal",labels_gp=gpar(fontsize=15),title_gp=gpar(fontsize=20)),
                 column_names_gp = gpar(fontsize = 15),
                 column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 row_title_gp = gpar(fontsize = 20, fontface = "bold"),
                 show_row_names = T, show_column_names = T,
                 cluster_row_slices = FALSE,
                 cluster_column_slices = FALSE,
                 cluster_rows = T, cluster_columns = T,
                 column_split = annot_df$First_level,
                 show_row_dend = F,show_column_dend = F,
                 column_dend_height = unit(2, "cm"), row_dend_width = unit(1, "cm")
)
draw(ht_list, heatmap_legend_side = "bottom",annotation_legend_side = "bottom")
dev.off()
p1<-draw(ht_list, heatmap_legend_side = "bottom",annotation_legend_side = "bottom")

aaa<-row_order(draw(ht_list))
bbb<-rownames(data_plot3_acid)[aaa]




#################################  Figure 6cd #################################################
### BGC
sj_type<-c("Hot","Pressure","Saline_alkaline","Acid","Cold")
bgclist<-rev(c("NRPS","PKSI","PKS-NRP_Hybrids","PKSother","RiPPs","Saccharides","Terpene","Others"))

all_all<-matrix(ncol=4,nrow=1)
colnames(all_all)<-c("Habitat","NMDC","BGC_Group","Region_Number")
for(mm in 1:5){
  data<-read.table(paste0("New bigscape/new_",sj_type[mm],"_all_antismash_bigscape_region_length_morethan5kb_addNMDC.txt"),header=T,sep="\t",quote='')
  data2<-aggregate(data$Region_filename,by=list(data$NMDC_ID,data$BiG.SCAPE_class),length)
  colnames(data2)<-c("NMDC","BGC_Group","Region_Number")
  data3<-cbind(matrix(rep(sj_type[mm],nrow(data2)),ncol=1),data2)
  colnames(data3)[1]<-"Habitat"
  all_all<-rbind(all_all,data3)
}
all_all$Habitat<-gsub("Saline_alkaline","Saline-alkaline",all_all$Habitat)  
all_all<-all_all[-1,]

drepdata<-read.table("Supplement_allfaa_KEGG_result_module_absence_use.txt",header = T,sep="\t",quote = '',check.names = F)
dff<-matrix(nrow=nrow(drepdata),ncol=10)
dff[,1]<-drepdata$Habitat
dff[,2]<-drepdata$MAG
colnames(dff)<-c("Habitat","NMDC",bgclist)
for(i in 1:nrow(dff)){
  for(j in 3:ncol(dff)){
    tempindex<-which(all_all$Habitat==dff[i,1] & all_all$NMDC==dff[i,2] & all_all$BGC_Group==colnames(dff)[j])
    if(length(tempindex)>0){
      dff[i,j]<-all_all[tempindex,4]
    }else{
      dff[i,j]<-0
    }
  }
}
dff2<-cbind(dff,matrix(paste(dff2[,1],dff2[,2],sep="_"),ncol=1))
colnames(dff2)[ncol(dff2)]<-"Info"
data_info<-drepdata[,1:4]
data_info2<-cbind(data_info,matrix(paste(data_info[,1],data_info[,2],sep="_"),ncol=1))
colnames(data_info2)[ncol(data_info2)]<-"Info"
dff3<-merge(dff2,data_info2,by="Info")
dff4<-dff3[,c(12:15,4:11)]
colnames(dff4)[1]<-"Habitat"
write.table(dff4,"All_representative_MAGs_BGC_region_numbers.txt",sep="\t",quote = F,col.names = T,row.names = F)

dff4_1<-dff4[which(dff4$Quality=="High-quality"),]
dff4_2<-dff4_1[grep(";s__$",dff4_1$GTDB),]
acid_data<-dff4_2[which(dff4_2$Habitat=="Acid"),]
acid_data_plot<-matrix(nrow=(nrow(acid_data)*8),ncol=3)
acid_data_plot[,1]<-rep(acid_data$MAG,each=8)
acid_data_plot[,2]<-rep(bgclist,times=nrow(acid_data))
colnames(acid_data_plot)<-c("NMDC","BGC_Group","Region_Number")
for(i in 1:nrow(acid_data_plot)){
  indexa<-which(acid_data$MAG==acid_data_plot[i,1])
  indexb<-which(colnames(acid_data)==acid_data_plot[i,2])
  acid_data_plot[i,3]<-acid_data[indexa,indexb]
}

acid_data_plot2<-acid_data_plot[which(!is.na(match(acid_data_plot[,1],bbb))),]
acid_data_plot2<-as.data.frame(acid_data_plot2)
acid_data_plot2$NMDC<-factor(acid_data_plot2$NMDC,levels = rev(bbb))
acid_data_plot2$BGC_Group<-factor(acid_data_plot2$BGC_Group,levels = rev(bgclist))
acid_data_plot2$Region_Number<-as.numeric(acid_data_plot2$Region_Number)

library(ggplot2)
palette<-c("NRPS"="#377EB8","PKSI"="#4DAF4A","PKS-NRP_Hybrids"="#984EA3","PKSother"="#FF7F00","RiPPs"="#E41A1C","Saccharides"="#A65628","Terpene"="#F781BF","Others"="#999999")
p2<-ggplot(acid_data_plot2,aes(y=NMDC,x=Region_Number,fill=BGC_Group))+geom_bar(stat = "identity",position="stack")+
  theme_classic()+xlab("Number of BGCs")+
  scale_fill_manual(values =palette)+
  theme(axis.line=element_line(color="black"),axis.text=element_text(size=12,color="black"),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title=element_text(size=15,color="black"),legend.text=element_text(size=12,color="black"),legend.title=element_text(size=12,color="black"),legend.position = "bottom",legend.box = "horizontal")


###ARO
data<-read.table("Habitat_5_drepMAGs_noassomge_assomge_aronum_use.txt",sep="\t",quote = '',header=T,check.names = F)
data_use<-data[which(!is.na(match(data$MAGs,bbb))),]
data_use2<-data_use[match(bbb,data_use$MAGs),]
rownames(data_use2)<-data_use2$MAGs
data_use2_1<-matrix(c(as.matrix(data_use2[,c(2,5)]),rep("NotAssoWithMGE",nrow(data_use2))),ncol=3)
data_use2_2<-matrix(c(as.matrix(data_use2[,c(2,6)]),rep("AssoWithMGE",nrow(data_use2))),ncol=3)
data_use3<-rbind(data_use2_1,data_use2_2)
colnames(data_use3)<-c("MAGs","Aro_Number","Group")
data_use3<-as.data.frame(data_use3)
data_use3$MAGs<-factor(data_use3$MAGs,levels = rev(bbb))
data_use3$Aro_Number<-as.numeric(data_use3$Aro_Number)

p3<-ggplot(data_use3,aes(y=MAGs,x=Aro_Number,fill=Group))+geom_bar(stat = "identity",position="stack")+
  theme_classic()+xlab("Number of Aro")+
  theme(axis.line=element_line(color="black"),axis.text=element_text(size=12,color="black"),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title=element_text(size=15,color="black"),legend.text=element_text(size=12,color="black"),legend.title=element_text(size=12,color="black"),legend.position = "bottom",legend.box = "horizontal")


library(patchwork)
p2+p3
ggsave(filename = "Acid_high-quality_species-unknown_MAGs_BGC_ARO_heatmap.pdf",width=12,height=26)
