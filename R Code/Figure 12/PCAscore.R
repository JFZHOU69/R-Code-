
library(limma)
library(ggpubr)
expFile="PCDGexp.txt"     
setwd("")    

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

pca=prcomp(data, scale=TRUE)
value=predict(pca)
PCDscore=value[,1]+value[,2]
PCDscore=as.data.frame(PCDscore)
scoreOut=rbind(id=colnames(PCDscore), PCDscore)
write.table(scoreOut, file="PCDscore.txt", sep="\t", quote=F, col.names=F)


PCDCluFile="cluster.txt"        
scoreFile="PCDscore.txt"          
    
PCDClu=read.table(PCDCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(PCDClu), row.names(score))
data=cbind(score[sameSample,,drop=F], PCDClu[sameSample,c("PCDcluster"),drop=F])

data$PCDcluster=factor(data$PCDcluster, levels=levels(factor(data$PCDcluster)))
group=levels(factor(data$PCDcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$PCDcluster)))]

boxplot=ggboxplot(data, x="PCDcluster", y="PCDscore", color="PCDcluster",
                  xlab="PCDcluster",
                  ylab="PCDscore",
                  legend.title="PCDcluster",
                  palette=bioCol,
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)


pdf(file="PCDcluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
