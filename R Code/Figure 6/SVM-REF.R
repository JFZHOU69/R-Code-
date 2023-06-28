
#set.seed(12345)
library(e1071)

inputFile="PCDGexp.txt"     
setwd("")     
source("msvmRFE.R")

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.data.frame(t(data))
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=cbind(group, data)
data$group=factor(data$group, levels=c("Control","Treat"))

svmRFE(data, k=10, halve.above=50)
nfold=10
geneNum=nrow(data)
folds=rep(1:nfold, len=geneNum)[sample(geneNum)]
folds=lapply(1:nfold, function(x) which(folds == x))
results=lapply(folds, svmRFE.wrap, data, k=10, halve.above=50)

top.features=WriteFeatures(results, data, save=F)
write.table(top.features, file="feature_svm.txt", sep="\t", quote=F,row.names=F)

featsweep=lapply(1:26, FeatSweep.wrap, results, data)

no.info=min(prop.table(table(data[,1])))
errors=sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

pdf(file="errors.pdf", width=5, height=5)
PlotErrors(errors, no.info=no.info)
dev.off()

pdf(file="accuracy.pdf", width=5, height=5)
Plotaccuracy(1-errors, no.info=no.info)
dev.off()

featureGenes=top.features[1:which.min(errors),1,drop=F]
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)
