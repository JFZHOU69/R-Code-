
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

diffFile="all.txt"           
gmtFile="Hallmark.all.v2022.1.Hs.symbols.gmt"    
setwd("")     

rt=read.table(diffFile, header=T, sep="\t", check.names=F)
rt=rt[order(rt$logFC, decreasing=T),]
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])

gmt=read.gmt(gmtFile)

kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff=1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

termNum=5   
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]       
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in treat group")
	pdf(file="GSEA.treat.pdf", width=8, height=6)
	print(gseaplot)
	dev.off()
}

termNum=5   
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]      
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in control group")
	pdf(file="GSEA.con.pdf", width=8, height=6)
	print(gseaplot)
	dev.off()
}


library(stringr)
#kk2@result$Description=gsub("HALLMARK_","",kk2@result$Description)
pdf(paste0("4.","all_ridgeplot_GSEA.pdf"),width = 6,height = 20)
ridgeplot(kk)
dev.off()

