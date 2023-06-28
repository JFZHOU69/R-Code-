library(ConsensusClusterPlus)      
expFile="PCDGexp.txt"          
workDir=""     
setwd(workDir)     

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="treat"]

maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123456,
                             plot="png")

calcICL(results, title="consensusScore", plot="png")

clusterNum=2       
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("PCDcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$PCDcluster))
cluster$PCDcluster=letter[match(cluster$PCDcluster, uniqClu)]
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="cluster.txt", sep="\t", quote=F, col.names=F)
