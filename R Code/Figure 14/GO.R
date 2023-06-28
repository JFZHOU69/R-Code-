
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
shownum=7            
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

kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="GO.xls",sep="\t",quote=F,row.names = F)
GO=GO[c(which(GO$ONTOLOGY=="BP")[1:shownum],which(GO$ONTOLOGY=="CC")[1:shownum],which(GO$ONTOLOGY=="MF")[1:shownum]),]

rt=GO[,c(1,2,3,10,4,6,8)]       
names(rt)=c("Ontology","ID","Term","Count","Ratio","pvalue","qvalue")

rt$ID=gsub(":","",rt$ID)
for (i in 1:nrow(rt)) {
  rt[i,3]=paste0(rt[i,2],":",rt[i,3])
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
  theme_bw()+
  facet_grid( Ontology~. ,scales="free")
ggsave("GO_bubble.pdf", width=7.5, height=6)      

rt=GO[,c(1,2,3,10,4,6,8)]        
names(rt)=c("Ontology","ID","Term","Count","Ratio","pvalue","qvalue")

rt$ID=gsub(":","",rt$ID)
for (i in 1:nrow(rt)) {
  rt[i,3]=paste0(rt[i,2],":",rt[i,3])
}

labels=rt[order(rt$qvalue,decreasing =T),"Term"]
rt$Term = factor(rt$Term,levels=labels)

p=ggplot(data=rt)+geom_bar(aes(x=Term, y=Count, fill=qvalue), stat='identity')+
  coord_flip() + scale_fill_gradient(high=highcol, low = lowcol) +    
  xlab("Term") + ylab("Gene count") +         
  theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) + 
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))+
  theme_bw()+
  facet_grid( Ontology~. ,scales="free")
print(p)
ggsave("GO_barplot.pdf", width=7.5, height=6)       
