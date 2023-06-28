
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(ggplot2)

setwd("") 
gene_symbol=read.table("diff.txt",sep="\t",check.names=F,header=T)
gene_name=as.vector(gene_symbol[,1])
geneID <- mget(gene_name, org.Hs.egSYMBOL2EG, ifnotfound=NA)
geneID <- as.character(geneID)
data=cbind(id=gene_symbol[,1],entrezID=geneID)
data=as.data.frame(data)
write.table(data,"name_id.txt",sep="\t",quote = F,row.names = F)
id=as.data.frame(data[,-2])
write.table(id,"id.txt",sep="\t",quote = F,row.names = F,col.names = F)

#GO
go <- enrichGO(gene=data$entrezID,
               OrgDb = org.Hs.eg.db, ont='ALL',pvalueCutoff = 0.05)#MF CC BP
write.csv(go,"go.csv",row.names =F)

pdf("GO1.pdf",12)
barplot(go, drop = T, showCategory =5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

pdf("GO2.pdf",9,9)
dotplot(go,showCategory = 5,split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale='free')
dev.off()
##KEGG
kegg <- enrichKEGG(gene = data$entrezID,organism ="human",pvalueCutoff = 0.05)
write.csv(kegg,"KEGG.csv",row.names =F)

pdf("KEGG1.pdf",10,8)
barplot(kegg, showCategory =10)
dev.off()

pdf("KEGG2.pdf",9,9)
dotplot(kegg,showCategory = 20)
dev.off()

#DO
do<-enrichDO(data$entrezID, ont = "DO",
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH")
write.csv(do,"DO.csv",row.names =F)

pdf("DO1.pdf",12,9)
barplot(do, showCategory =20)
dev.off()

pdf("DO2.pdf",9,9)
dotplot(do,showCategory = 20)
dev.off()

