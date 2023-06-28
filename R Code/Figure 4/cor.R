
options(stringsAsFactors=F)
library(corrplot)
library(circlize)
library(limma)
library(PerformanceAnalytics)
setwd("")

expFile="normalize.txt"     
hub="PCD.txt"         
C="C"                     
higcol="red"                   
midcol="white"               
lowcol="green"                
showtype="upper"     
tlcex=0.45         
numbercex=0.5       
method="pearson"  

rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
rt=rt[,rownames(sample)[-(1:sum(sample[,1]==C))]]

VENN=intersect(rownames(rt),read.table(hub,sep = "\t",header = F,check.names = F)[,1])
rt1=rt[c(VENN),]

rt=apply(rt1[,-1],2,as.numeric)
rownames(rt)=rownames(rt1)
rt=t(rt) 
M=cor(rt,method = method)   
res <- cor.mtest(rt)

pdf(file="corpot1.pdf",width=7,height=7)
corrplot(M,
         method = "circle",
         order = "hclust", 
         type=showtype,
         col=colorRampPalette(c(lowcol, midcol, higcol))(50)
)
dev.off()

pdf(file="corpot2.pdf",width=8,height=8)
corrplot(M,
         order="original",
         method = "color",
         number.cex = 0.7, 
         addCoef.col = "black",
         diag = TRUE,
         type=showtype,
         tl.col="black",
         col=colorRampPalette(c(lowcol, midcol, higcol))(50))
dev.off()

pdf(file="corpot3.pdf",width=7,height=7)
corrplot(M,
         type="upper",
         tl.pos="tp",
         order="AOE",
         col=colorRampPalette(c(lowcol, midcol, higcol))(50),
         tl.cex = tlcex)
corrplot(M,add=TRUE, 
         type="lower",
         method="number",
         order="AOE",
         diag=FALSE,
         tl.pos="n",
         cl.pos="n",
         col=colorRampPalette(c(lowcol, midcol, higcol))(50),
         number.cex=numbercex
         )
dev.off()

pdf(file="corpot4.pdf",width=7,height=7)
corrplot(
  M, 
  order = 'AOE',method = 'pie',
  type = 'lower', 
  tl.pos = 'd',
  col=colorRampPalette(c(lowcol, midcol, higcol))(50),
  tl.cex = tlcex)
corrplot(
  M, 
  add = TRUE, 
  type = 'upper',  
  method = 'number',
  order = 'AOE', 
  diag = FALSE,  
  tl.pos = 'n', 
  cl.pos = 'n',
  col=colorRampPalette(c(lowcol, midcol, higcol))(50),
  number.cex=numbercex)
dev.off()

pdf(file="corpot5.pdf",width=8,height=8)
chart.Correlation(rt,method = method)
dev.off()

pdf(file="corpot6.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type=showtype,
         tl.cex=0.8, pch=T,
         p.mat = res$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c(lowcol, midcol, higcol))(50),
         tl.col="black")
dev.off()

col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
M[M==1]=0  
c1 = ifelse(c(M)>=0,rgb(1,0,0,abs(M)),rgb(0,1,0,abs(M)))
col1 = matrix(c1,nc=ncol(rt))
pdf(file="corpot7.pdf",width=8,height=8)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(M)-1)), start.degree = 180)
chordDiagram(M, grid.col=rainbow(ncol(rt)), col=col1, transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))      
dev.off()
circos.clear()
