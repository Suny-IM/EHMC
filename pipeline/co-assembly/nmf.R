library(NMF)
library(cluster)
library(ConsensusClusterPlus)

data <- read.delim("combine.g.map.nohomo.txt",header = T, row.names = 1)
row<-rownames(data)
col<-colnames(data)
data2<-as.matrix(sapply(data,as.numeric))
data3<-normalize.quantiles(data2)
rownames(data3)<-row
colnames(data3)<-col
write.table(data3, "combine.g.normalize.map.nohomo.txt", sep = "\t")

data <- read.delim("combine.g.normalize.map.nohomo.txt",header = T, row.names = 1)

##different clusters
res <- nmf(data, 2:25, 'nsNMF', maxIter=1000)
pdf("consensusmap.pdf",width = 22, height=15)
consensusmap(res)
dev.off()

pdf("nmf_rank_survey.pdf")
plot(res)
dev.off()
write.table(res$measures, "measures.txt", sep = "\t")
summary(res)
