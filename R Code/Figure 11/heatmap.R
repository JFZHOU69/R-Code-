
library(pheatmap)
inputFile="ssGSEA.result.txt"      
setwd("")    

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

con=grepl("_con", colnames(rt), ignore.case=T)
treat=grepl("_treat", colnames(rt), ignore.case=T)
conData=rt[,con]
treatData=rt[,treat]
conNum=ncol(conData)
treatNum=ncol(treatData)
data=cbind(conData,treatData)

Type=c(rep("Normal",conNum), rep("OA",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=8, height=5)
pheatmap(data, 
         annotation=Type, 
         color=colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
         cluster_cols=F,
         show_colnames=F,
         scale="row",
         fontsize = 6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()

