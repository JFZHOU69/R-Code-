
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)

pFilter=0.05
adjPfilter=0.05

inputFile="input.txt"
setwd("")           
rt=read.table(inputFile,sep="\t",header=T,check.names=F)       
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)       
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
rt=rt[is.na(rt[,"entrezID"])==F,]                              
gene=rt$entrezID

GO=enrichGO(gene = gene,
            OrgDb = org.Hs.eg.db, 
            pvalueCutoff =1, 
            qvalueCutoff = 1,
            ont="all",
            readable =T)
GO=as.data.frame(GO)
GO=GO[(GO$pvalue<pFilter & GO$p.adjust<adjPfilter),]
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)       

go=data.frame(Category = GO$ONTOLOGY,ID = GO$ID,Term = GO$Description, Genes = gsub("/", ", ", GO$geneID), adj_pval = GO$p.adjust)

genelist=data.frame(ID = rt$gene, logFC = rt$logFC)
row.names(genelist)=genelist[,1]

circ <- circle_dat(go, genelist)
termNum = 11                                     
termNum=ifelse(nrow(go)<termNum,nrow(go),termNum)
geneNum = nrow(genelist)                       


chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="GOcircos.pdf",width = 15,height =17 )
GOChord(chord, 
        space = 0.0000001,           
        gene.order = 'logFC',   
        gene.space = 0.20,       
        gene.size = 6,          
        border.size = 0.1,     
        process.label = 7)      
dev.off()


pdf(file="GOcluster.pdf",width = 15,height = 14)
GOCluster(circ, as.character(go[1:termNum,3]))
dev.off()

