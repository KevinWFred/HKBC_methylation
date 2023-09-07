#!/usr/bin/env Rscript
library(data.table)
library("readxl")
setwd("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code")
#read DNA methylation
meth1=as.data.frame(fread("../data/HKBC.methylation.set1.after.QC.txt"))
rownames(meth1)=meth1$probeID
#the last column is probeID
meth1=meth1[,-ncol(meth1)]
sampletable1=read.csv("../data/dataMapping_first_HKBC_set.csv")
dim(sampletable1) #189 7
sum(colnames(meth1)%in% sampletable1$labID) #189
meth1=meth1[,match(sampletable1$labID,colnames(meth1))]
meth1_NCIID=sampletable1$subjectID

meth2=as.data.frame(fread("../data/HKBC.methylation.set2.after.QC.txt"))
rownames(meth2)=meth2$probeID
meth2=meth2[,-ncol(meth2)]
tmp=unlist(strsplit(colnames(meth2),"_"))
tmp1=tmp[seq(1,length(tmp),3)]
colnames(meth2)=tmp1
all(rownames(meth1)==rownames(meth2))
sampletable2=read_excel("../data/HongKongBreast_EPICmeth_ID_Mapping_16Aug2021.xls")
dim(sampletable2) #198 11
sum(sampletable2$`CGR SampleId` %in% colnames(meth2)) #198
meth2=meth2[,match(sampletable2$`CGR SampleId`,colnames(meth2))]
meth2_NCIID=sampletable2$NCI_ID

#Some of the samples have tree categories
wgs_normal=read_excel("../data/WGS-Normal database.xls",sheet = "data")
wgs_normal=wgs_normal[wgs_normal$`NCI ID` %in% c(meth1_NCIID,meth2_NCIID),]
wgs_normal=wgs_normal[!is.na(wgs_normal$`Tree categories`),]
wgs_tumor=read_excel("../data/WGS-Tumor database.xls",sheet = "data")
wgs_tumor=wgs_tumor[wgs_tumor$NCI_ID %in% c(meth1_NCIID,meth2_NCIID),]
wgs_tumor=wgs_tumor[!is.na(wgs_tumor$`Tree categories`),]
idx1=which(wgs_tumor$NCI_ID %in% sampletable1$subjectID[sampletable1$Tissue=="Tumor"]) #10
idx2=which(wgs_tumor$NCI_ID %in% sampletable2$NCI_ID[sampletable2$`Attribute Value`=="Tumor"]) #79
Tree_tumorsample=wgs_tumor$NCI_ID[c(idx1,idx2)]
idx1=which(wgs_normal$`NCI ID` %in% sampletable1$subjectID[sampletable1$Tissue=="Normal"]) #9
idx2=which(wgs_normal$`NCI ID` %in% sampletable2$NCI_ID[sampletable2$`Attribute Value`=="Normal"]) #33
Tree_normalsample=wgs_normal$`NCI ID`[c(idx1,idx2)]
length(intersect(Tree_tumorsample,Tree_normalsample)) #41
wgs_tumor=wgs_tumor[match(wgs_normal$`NCI ID`,wgs_tumor$NCI_ID),]

#combine sample table
allsampletable=rbind(data.frame(SampleID=colnames(meth1),NCIID=meth1_NCIID,Tissue=sampletable1$Tissue,Tree=NA,Batch=1),
                     data.frame(SampleID=colnames(meth2),NCIID=meth2_NCIID,Tissue=sampletable2$`Attribute Value`,Tree=NA,Batch=2))
idx=which(allsampletable$NCIID %in% wgs_normal$`NCI ID` & allsampletable$Tissue=="Tumor")
idx1=match(allsampletable$NCIID[idx],wgs_normal$`NCI ID`)
allsampletable$Tree[idx[idx1]]=wgs_normal$`Tree categories`
idx=which(allsampletable$NCIID %in% wgs_normal$`NCI ID` & allsampletable$Tissue=="Normal")
idx1=match(allsampletable$NCIID[idx],wgs_normal$`NCI ID`)
allsampletable$Tree[idx[idx1]]=wgs_normal$`Tree categories`

#combine methylation
allmeth=cbind(meth1,meth2)
#there are 3 0 beta values
idx=which(allmeth==0,arr.ind = T)
allmeth[idx]=1e-10
library(lumi)
allmvalue=beta2m(allmeth)

#to remove batch effect
library(sva)
#At same time to preserve the variation due to Tissue difference
modcombat<-model.matrix(~Tissue, data=allsampletable)
allmvalue.combat= ComBat(dat=as.matrix(allmvalue), batch=allsampletable$Batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
allmeth.combat=m2beta(allmvalue.combat)
save(allmeth.combat,file="../result/DNAmeth_combat.RData")

#pca plot for two datasets
plot_2pca=function(dat1=meth_batch1,dat2=meth_batch2)
{
  
  tmp=intersect(rownames(dat1),rownames(dat2))
  
  dat1=dat1[match(tmp,rownames(dat1)),]
  
  dat2=dat2[match(tmp,rownames(dat2)),]
  
  dat=cbind(dat1,dat2)
  
  pcadat=prcomp(t(dat),scale = F)
  
  #plot(pcadat$sdev^2/sum(pcadat$sdev^2)*100,ylab="Variance explained (%)",xlab="PC",cex.axis=1.3,cex.lab=1.3)
  
  pcadat=pcadat$x
  
  mycol=c(rep("blue",ncol(dat1)),rep("green",ncol(dat2)))
  
  plot(pcadat[,1],pcadat[,2],ylim=c(min(pcadat[,2]),2*max(pcadat[,2])),col=mycol,xlab="PC1",ylab="PC2",cex.axis=1.2,cex.lab=1.2)
  
  legend("topleft",legend=c("Batch1","Batch2"),col=c("blue","green"),cex=0.8,pch=1)
}

tumor_meth_batch1=allmeth[,allsampletable$Batch==1 & allsampletable$Tissue=="Tumor"]
tumor_meth_batch2=allmeth[,allsampletable$Batch==2 & allsampletable$Tissue=="Tumor"]
plot_2pca(dat1=tumor_meth_batch1,dat2=tumor_meth_batch2)
combat_tumor_meth_batch1=allmeth.combat[,allsampletable$Batch==1 & allsampletable$Tissue=="Tumor"]
combat_tumor_meth_batch2=allmeth.combat[,allsampletable$Batch==2 & allsampletable$Tissue=="Tumor"]
plot_2pca(dat1=combat_tumor_meth_batch1,dat2=combat_tumor_meth_batch2)
