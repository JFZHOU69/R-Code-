library(glmnet)
library(pROC)

expFile="nomoscore.txt"      
geneFile="interGenes.txt"     
setwd("")   
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1,)
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="con", 0, 1)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
bioCol=rainbow(nrow(geneRT), s=0.9, v=0.9)    
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
	k=k+1
	roc1=roc(y, as.numeric(rt[x,]))    
	if(k==1){
		pdf(file="ROC.genes.pdf", width=5, height=5)
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}else{
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}
}
legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()

rt=rt[as.vector(geneRT[,1]),]
rt=as.data.frame(t(rt))
logit=glm(y ~ ., family=binomial(link='logit'), data=rt)
pred=predict(logit, newx=rt)     

roc1=roc(y, as.numeric(pred))      
ci1=ci.auc(roc1, method="bootstrap")    
ciVec=as.numeric(ci1)
pdf(file="ROC.model.pdf", width=5, height=4.75)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Model")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()
