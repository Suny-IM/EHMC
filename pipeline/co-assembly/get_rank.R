library(NMF)
library(cluster)
library(ConsensusClusterPlus)
library(pheatmap)
library(randomcoloR)


data <- read.delim("combine.g.normalize.nohomo.map.txt",header = T, row.names = 1)

args <- commandArgs(trailingOnly = TRUE)

num <- as.numeric(args[1])

res2 <- nmf(data, rank=num, 'nsNMF',seed=123456)

group1 <- predict(res2,'columns')

palette1 <- distinctColorPalette(num)
palette2 <- distinctColorPalette(num)
palette3 <- distinctColorPalette(num)


pdf("heatmap_cluster.pdf", width = 12, height=6)
coefmap(res2,annRow=list(Signature = ":basis"),annCol=list(':basis',Group=group1),annColors=list(Signature=palette1,Group=palette2))
basismap(res2,annRow=list(Signature = ":basis"),annColors=list(Signature=palette1))
dev.off()


group <- as.data.frame(group1)
group$group1 <- paste0('Cluster',group$group)
write.csv(group,'expr_group.csv')


colorname <-c(Cluster=palette3)
ann_color = list(group1=colorname)


h<-coef(res2)



test<-t(h)

pheatmap(test,scale='row',col=colorRampPalette(rev(c("red","white","blue")))(100),cluster_cols=F,main='nsNMF_rank_normalize',fontsize_row=6,width=5,height=max(5,nrow(test)*0.1+3),annotation_row=group,annotation_colors = ann_color,filename="nsNMF_rank9_normalize.png")

pheatmap(test,scale='row',col=colorRampPalette(rev(c("red","white","blue")))(100),cluster_cols=F,main='nsNMF_rank_normalize',fontsize_row=6,width=5,show_rownames = F,annotation_row=group,annotation_colors = ann_color,filename="nsNMF_rank9_normalize_small.png")

