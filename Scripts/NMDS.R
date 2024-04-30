setwd("YouPath")
library(vegan)
library(ggplot2)
library("ggsci")
library(ggpubr)
library(ggalt)
library(dplyr)

##########unique lineage##########
data<-read.table("combine5_unique_lineage.txt",header=T,row.names = 1,sep="\t",na.strings = "NA",encoding = "UTF-8",fill=T,comment.char = "",quote = '')
zero_counts_per_row <- apply(data, 1, function(row) sum(row == 0))
rows_to_remove <- which(zero_counts_per_row >1100)
data <- data[-rows_to_remove, ]

##########core lineage############
data<-read.table("combine5_core_lineage.txt",header=T,row.names = 1,sep="\t",na.strings = "NA",encoding = "UTF-8",fill=T,comment.char = "",quote = '')

#################################
data1<- data.frame(t(data))
group<-read.delim("top10_isolate_source_group.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
group<-group[,c(1,2,5)]

data2 <- data1[(rownames(data1) %in% group[,1]), ]
data2<-data2[which(rowSums(data2) > 0),]

distance <- vegdist(data2, method = 'bray')

set.seed(123)

df_nmds <- metaMDS(distance, k =2,try=10,trymax=20)

df_nmds_stress <- df_nmds$stress
#stressplot(df_nmds)

df_points <- as.data.frame(df_nmds$points)

df_points$samples <- row.names(df_points)

names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
colnames(group) <- c("samples","Group","Habitat")
df <- merge(df_points,group,by="samples")

df$Group <- factor(df$Group, levels = c('Antarctic saline lake','Permafrost','Hot spring','Hydrothermal related',
                                        'Acid mine drainage','Marine','Saline lake','Halite','Antarctic soil and sand','Cryoconite','Others'))

library(RColorBrewer)
library(scales)

color=c("#1F77B4FF","#2CA02CFF","#C0C0C0","#D62728FF","#FF7F0EFF","#9467BDFF")

pdf(file = "·ÇÈßÓà»ùÒò¾ØÕó-Í¼/new-core_gene_lineage_stress0.188_atlest10samples20240406.pdf",width =6,height = 5)

ggplot(data=df,aes(x=NMDS1,y=NMDS2))+
  theme_bw()+
  geom_point(aes(color = Habitat), size=1)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", linewidth = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", linewidth = 1, color = 'grey50')+
  scale_color_manual(values = color) +
  scale_fill_manual(values = color)+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid=element_blank())+
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))
  
dev.off()

######PERMANOVA###########
adonis <-  adonis2(data2 ~ Habitat, data =group, permutations = 999, method = "bray")
print(adonis)

#####################################################
####################subgroup#########################
#######unique lineage######
data<-read.table("combine5_unique_lineage.txt",header=T,row.names = 1,sep="\t",na.strings = "NA",encoding = "UTF-8",fill=T,comment.char = "",quote = '')
zero_counts_per_row <- apply(data, 1, function(row) sum(row == 0))
rows_to_remove <- which(zero_counts_per_row >1100)
data <- data[-rows_to_remove, ]

#######core lineage######
data<-read.table("combine5_core_lineage.txt",header=T,row.names = 1,sep="\t",na.strings = "NA",encoding = "UTF-8",fill=T,comment.char = "",quote = '')

data1<- data.frame(t(data))
group<-read.delim("Five_raw_group.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)

group1<-group[which(group$Isolate_Source1=="Saline lake"|group$Isolate_Source1=="Antarctic saline lake"),]
#group1<-group[which(group$Isolate_Source1=="Permafrost"),]
#group1<-group[which(group$Isolate_Source1=="Hot spring"),]
#group1<-group[which(group$Isolate_Source1=="Hydrothermal related"),]
#group1<-group[which(group$Isolate_Source1=="Acid mine drainage"),]
a<-matrix(c(names(table(group1$Isolate_Source3)),table(group1$Isolate_Source3)),ncol=2)
a<-as.data.frame(a)
b<-a[order(as.numeric(a$V2),decreasing = T),]
b<-b[1:10,]
group2 <- merge(x = group1, y = b, by.x = 4, by.y = 1)
group<-group2[,c(2,1,5)]

data2 <- data1[(rownames(data1) %in% group[,1]), ]
data2<-data2[which(rowSums(data2) > 0),]

distance <- vegdist(data2, method = 'bray')

set.seed(123)
df_nmds <- metaMDS(distance, k =2,try=10,trymax=20)

df_nmds_stress <- df_nmds$stress
#stressplot(df_nmds)

df_points <- as.data.frame(df_nmds$points)

df_points$samples <- row.names(df_points)

names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
colnames(group) <- c("samples","Group","Habitat")
df <- merge(df_points,group,by="samples")


 subgroup<-unique(group$Group)
 df$Group <- factor(df$Group, levels = c(subgroup))


library(RColorBrewer)
library(scales)

color=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00",
        "#CAB2D6","#6A3D9A")

pdf(file = "core_gene_lineages_Saline_lake.pdf",width =10,height = 6)
ggplot(data=df,aes(x=NMDS1,y=NMDS2))+
  theme_bw()+
  geom_point(aes(color = Group,shape=Habitat), size=1)+
  scale_shape_manual(values = c(16,1,18))+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", linewidth = 0.5, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", linewidth = 0.5, color = 'grey50')+
  scale_color_manual(values = color) +
  scale_fill_manual(values = color)+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid=element_blank())+
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))

dev.off()
write.table(df,"core_gene_lineage_Saline_lake.txt",sep="\t",col.names = T,row.names = F,quote = F)

