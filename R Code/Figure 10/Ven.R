
library(VennDiagram)
setwd("")
myfiles=dir()
myfiles=grep("txt",myfiles,value=T)
mylist=list()
for(i in 1:length(myfiles)){
    onefile=myfiles[i]
    data=read.table(onefile,header=F)
    hd=unlist(strsplit(onefile,"\\.|\\-"))
    mylist[[hd[1]]]=as.vector(data[,1])
    onel=length(unique(as.vector(data[,1])))
    print(paste(hd[1],onel,sep=" "))
}

vedata=Reduce(intersect,mylist)
write.table(file="Venn.txt",vedata,sep="\t",quote=F,col.names=F,row.names=F)
venn.plot <-venn.diagram(mylist,fill=c("red","blue","yellow"), filename=NULL)
pdf("venn.pdf");
grid.draw(venn.plot);
dev.off()