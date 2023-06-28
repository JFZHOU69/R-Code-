
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")
library("ggnewscale")
library("DOSE")
library(stringr)

pvalueFilter=0.05      
qvalueFilter=0.05      
showNum=20            
highcol="grey"      
lowcol="red3"          

rt=read.table("intersectGenes.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)
R.utils::setOption( "clusterProfiler.download.method",'auto' )

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.xls",sep="\t",quote=F,row.names = F)

rt=KEGG[1:showNum,c(1,2,9,3,5,7)]       
names(rt)=c("ID","Term","Count","Ratio","pvalue","qvalue")

for (i in 1:nrow(rt)) {
  rt[i,2]=paste0(rt[i,1],":",rt[i,2])
}

split_b<-str_split(rt$Ratio,"/")
b<-sapply(split_b,"[",1)
c<-sapply(split_b,"[",2)
rt$Ratio=as.numeric(rt$Count)/as.numeric(c[1])

labels=rt[order(rt$Ratio),"Term"]
rt$Term = factor(rt$Term,levels=labels)

p = ggplot(rt,aes(Ratio, Term)) + 
  geom_point(aes(size=Count, color=qvalue))
p1 = p + 
  scale_colour_gradient(high=highcol, low = lowcol) +
  labs(color="qvalue",size="Count",x="Gene ratio",y="Term")+     
  theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) + 
  scale_size_continuous(range=c(4,9))+     
  theme_bw()
ggsave("KEGG_bubble.pdf", width=7.5, height=6)      

rt=KEGG[1:showNum,c(1,2,9,3,5,7)]          
names(rt)=c("ID","Term","Count","Ration","pvalue","qvalue")

for (i in 1:nrow(rt)) {
  rt[i,2]=paste0(rt[i,1],":",rt[i,2])
}

labels=rt[order(rt$qvalue,decreasing =T),"Term"]
rt$Term = factor(rt$Term,levels=labels)

p=ggplot(data=rt)+geom_bar(aes(x=Term, y=Count, fill=qvalue), stat='identity')+
  coord_flip() + scale_fill_gradient(high=highcol, low = lowcol) +   
  xlab("Term") + ylab("Gene count") +         
  theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) + 
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))+
  theme_bw()
print(p)
ggsave("KEGG_barplot.pdf", width=7.5, height=6)       

