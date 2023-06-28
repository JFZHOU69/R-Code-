
library(vioplot)                 
inputFile="ssGSEA.result.txt"     
setwd("")    

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

con=grepl("_con", colnames(rt), ignore.case=T)
treat=grepl("_treat", colnames(rt), ignore.case=T)
conData=rt[,con]
treatData=rt[,treat]
conNum=ncol(conData)
treatNum=ncol(treatData)
rt=t(cbind(conData,treatData))

outTab=data.frame()
pdf(file="vioplot.pdf", width=18, height=8.5)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,(3*ncol(rt)-3)),ylim=c(min(rt),max(rt)+0.05),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

for(i in 1:ncol(rt)){
	  if(sd(rt[1:conNum,i])==0){
	    rt[1,i]=0.00001
	  }
	  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
	    rt[(conNum+1),i]=0.00001
	  }
	  conData=rt[1:conNum,i]
	  treatData=rt[(conNum+1):(conNum+treatNum),i]
	  vioplot(conData,at=3*(i-1),lty=1,add = T,col = 'blue')
	  vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  wilcoxTest=wilcox.test(conData,treatData)
	  p=wilcoxTest$p.value
	  if(p<0.05){
	      cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
		  outTab=rbind(outTab,cellPvalue)
	  }
	  mx=max(c(conData,treatData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("Normal", "OA"),
       lwd=3,bty="n",cex=1,
       col=c("blue","red"))
text(seq(1,(3*ncol(rt)-2),3),-0.075,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()

write.table(outTab,file="immuneDiff.xls",sep="\t",row.names=F,quote=F)
