library(dplyr)
library(pROC)
library(ggplot2)
library(survival)
library(regplot)
library(rms)
library(ggsci)
library(survminer)
library(timeROC)
library(ggDCA)
library(limma)
library(rmda)
inputF="PCDexp.txt" 
inputFile="normalize.txt"       
hub="Hub.txt"        
setwd("")
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
sample=read.table("sample.txt",sep="\t",header=F,check.names=F)
colnames(sample)=c("ID","Type")
data=data[sample$ID,]
aSAH1=data[,read.table(hub, header=F, sep="\t", check.names=F)[,1]]
aSAH=cbind(sample,aSAH1)

aflist=roc(Type~TNFAIP3+JUN+PPP1R15A+INHBB+DDIT4, data = aSAH)
g3 <- ggroc(aflist, size = 1.2,alpha=.6,)
g5=g3+ggsci::scale_color_lancet()
print(g5)

dd <- datadist(aSAH)
options(datadist="dd")
fit <- lrm(formula = Type ~ TNFAIP3+JUN+PPP1R15A+INHBB+DDIT4, data =aSAH, x=T, y=T,)
print(fit)
coef=as.data.frame(fit$coefficients)[-1,,drop=F]
coefout=cbind(ID=rownames(coef),coef)
write.table(coefout,file="coefficients.txt",sep="\t",quote=F,row.names = F)

pdf(file="nomogram.pdf", width=9, height=7.5)
plot(nomogram(fit,fun.at = seq(0.05,0.95,0.05)),funlabel = "nomogram model")
dev.off()

plot(regplot(fit,plots=c("density","boxes"), observation=T, title="Prediction Nomogram", clickable=T, points=TRUE,droplines=TRUE))

nomoscore=predict(fit, data=t(aSAH))
aSAH$nomoscore=nomoscore
write.table(aSAH,file="nomoscore.txt",sep="\t",quote=F,row.names = F)

cali=calibrate(fit, method="boot", B=1000)
pdf("Calibration.pdf", width=6, height=6)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()

data=read.table(inputF, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")

ddist=datadist(rt)
options(datadist="ddist")

lrmModel=lrm(Type~ TNFAIP3+JUN+PPP1R15A+INHBB+DDIT4, data=rt, x=T, y=T, )
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.0001,0.1,0.3,0.5,0.7,0.9,0.99),
              lp=F, funlabel="Risk of Disease")

pdf("Nom.pdf", width=8, height=6)
plot(nomo)
dev.off()

cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=6, height=6)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()

rt$Type=ifelse(rt$Type=="con", 0, 1)
dc=decision_curve(Type ~ TNFAIP3+JUN+PPP1R15A+INHBB+DDIT4, data=rt, 
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)

pdf(file="DCA.pdf", width=6, height=6)
plot_decision_curve(dc,
                    curve.names="genes",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()

pdf(file="clinical_impact.pdf", width=6, height=6)
plot_clinical_impact(dc,
                     confidence.intervals=T,
                     col = c("red", "blue"))
dev.off()
