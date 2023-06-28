
library(corrplot)     
inputFile="ssGSEA.result.txt"    
setwd("")     

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
rt=data[,group=="treat",drop=F]

rt=t(rt)
M=cor(rt)

res1=cor.mtest(rt, conf.level = .95)

pdf(file="corpot.pdf", width=12, height=12)
corrplot(M,
         type = "upper",      
         method = "circle", 
         col=colorRampPalette(c('blue', 'white', 'red'),alpha = TRUE)(100), tl.pos="lt",   
         p.mat=res1$p, insig="label_sig", sig.level = c(.001, .01, .05), pch.cex = 0.85)  
corrplot(M, type="lower", add=TRUE, method="number",col=colorRampPalette(c('blue', 'white', 'red'), alpha = TRUE)(100), tl.pos = "n", cl.pos="n", diag=FALSE, number.cex = 0.6) 
dev.off()
