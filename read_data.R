#!/usr/bin/env Rscript
library(data.table)
library("readxl")
setwd("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code")
meth1=as.data.frame(fread("../data/HKBC.methylation.set1.after.QC.txt"))
rownames(meth1)=meth1$probeID
meth1=meth1[,-ncol(meth1)]
meth2=as.data.frame(fread("../data/HKBC.methylation.set2.after.QC.txt"))
rownames(meth2)=meth2$probeID
meth2=meth2[,-ncol(meth2)]
all(rownames(meth1)==rownames(meth2))
tmp=unlist(strsplit(colnames(meth2),"_"))
tmp1=tmp[seq(1,length(tmp),3)]
which(duplicated(tmp1))
colnames(meth2)=tmp1
idx=duplicated(tmp1)
meth2=meth2[,!idx]
#meth=cbind(meth1,meth2)
pheno=read_excel("../data/WGS-Tumor database.xls",sheet = "data")
sampletable1=read.csv("../data/dataMapping_first_HKBC_set.csv")
sum(colnames(meth1)%in% sampletable1$labID) #189
meth1=meth1[,colnames(meth1)%in% sampletable1$labID]
idx=match(sampletable1$labID[sampletable1$Tissue=="Tumor"],colnames(meth1))
meth1_tumor=meth1[,idx]
colnames(meth1_tumor)=sampletable1$subjectID[sampletable1$Tissue=="Tumor"]
idx=match(sampletable1$labID[sampletable1$Tissue=="Normal"],colnames(meth1))
meth1_normal=meth1[,idx]
colnames(meth1_normal)=sampletable1$subjectID[sampletable1$Tissue=="Normal"]


sampletable2=read_excel("../data/HongKongBreast_EPICmeth_ID_Mapping_16Aug2021.xls")
sum(colnames(meth2) %in% sampletable2$`CGR SampleId`) #198
meth2=meth2[,colnames(meth2) %in% sampletable2$`CGR SampleId`]
idx=match(sampletable2$`CGR SampleId`[sampletable2$`Attribute Value`=="Tumor"],colnames(meth2))
meth2_tumor=meth2[,idx]
colnames(meth2_tumor)=sampletable2$NCI_ID[sampletable2$`Attribute Value`=="Tumor"]
meth2=meth2[,colnames(meth2) %in% sampletable2$`CGR SampleId`]
idx=match(sampletable2$`CGR SampleId`[sampletable2$`Attribute Value`=="Normal"],colnames(meth2))
meth2_normal=meth2[,idx]
colnames(meth2_normal)=sampletable2$NCI_ID[sampletable2$`Attribute Value`=="Normal"]

sampletable=read_excel("../data/HKBC_Trio_analysis.xlsx")
#2 samples removed (HKB1298 HKB2012) (SC335668 SC335672)

wgs_tumor=read_excel("../data/WGS-Tumor database.xls",sheet = "data")
sum(wgs_tumor$NCI_ID %in% sampletable1$subjectID) #10
sum(wgs_tumor$NCI_ID %in% sampletable2$NCI_ID) #79

Tree_samples_tumor=wgs_tumor$NCI_ID[!is.na(wgs_tumor$`Tree categories`)]
wgs_normal=read_excel("../data/WGS-Normal database.xls",sheet = "data")
Tree_samples_normal=wgs_normal$`NCI ID`[!is.na(wgs_normal$`Tree categories`)] #42 subjects have Tree categories
all(Tree_samples_tumor %in% Tree_samples_normal) 
sum(Tree_samples_tumor %in% colnames(meth1_tumor)) #9
sum(Tree_samples_tumor %in% colnames(meth2_tumor)) #32
Tree_samples_tumor[!Tree_samples_tumor %in% c(colnames(meth1_tumor),colnames(meth2_tumor))] #"2087" tumor is missing
sum(Tree_samples_normal %in% colnames(meth1_normal)) #9
sum(Tree_samples_normal %in% colnames(meth2_normal)) #32
Tree_samples_normal[!Tree_samples_normal %in% c(colnames(meth1_normal),colnames(meth2_normal))] #"2087" normal is missing

tumor_meth=cbind(meth1_tumor,meth2_tumor)
tumor_meth=tumor_meth[,colnames(tumor_meth) %in% Tree_samples_tumor]
idx=match(colnames(tumor_meth),wgs_tumor$NCI_ID)
wgs_tumor=wgs_tumor[idx,]
tumor_pheno=data.frame(NCI_ID=wgs_tumor$NCI_ID,AIMS=factor(wgs_tumor$AIMS),AGE=wgs_tumor$AGE,Tree=factor(wgs_tumor$`Tree categories`,ordered=T))
tumor_pheno$batch="batch1"
tumor_pheno$batch[which(tumor_pheno$NCI_ID %in% colnames(meth2_tumor))]="batch2"
tumor_pheno$batch=as.factor(tumor_pheno$batch)
sum(is.na(tumor_pheno$AIMS)) # 2 AIMS missing
normal_meth=cbind(meth1_normal,meth2_normal)
normal_meth=normal_meth[,colnames(normal_meth) %in% Tree_samples_normal]
idx=match(colnames(normal_meth),wgs_normal$`NCI ID`)
wgs_normal=wgs_normal[idx,]
normal_pheno=data.frame(NCI_ID=wgs_normal$`NCI ID`,AIMS=factor(wgs_normal$AIMS_T),AGE=wgs_normal$AGE,Tree=factor(wgs_normal$`Tree categories`,ordered=T))
normal_pheno$batch="batch1"
normal_pheno$batch[which(normal_pheno$NCI_ID %in% colnames(meth2_normal))]="batch2"
normal_pheno$batch=as.factor(normal_pheno$batch)

normal_meth=normal_meth[,match(colnames(tumor_meth),colnames(normal_meth))]
normal_pheno=normal_pheno[match(tumor_pheno$NCI_ID,normal_pheno$NCI_ID),]

library(lumi)
tumor_methM=beta2m(tumor_meth)
normal_methM=beta2m(normal_meth)


removeconstrows=function(dat)
{
  tmp=apply(dat,1,sd)
  idxconst=tmp==0
  dat=dat[!idxconst,]
  return(dat) 
}

pca=function(dat=tumor_meth,check=F,bigpca=F,numpc=NULL,seed=1000)
{
  if (bigpca==T)
  {
    #use random rows
    set.seed(seed)
    idx=sample(1:nrow(dat))[1:floor(nrow(dat)*0.2)]
    dat=as.data.frame(dat[idx,])
    
  }
  
  if (check==T)
  {
    #remove constant rows
    dat=removeconstrows(dat)
    
  }
  
  pc=prcomp(t(dat),scale=F,center=T)
  varprop=pc$sdev^2/sum(pc$sdev^2)
  plot(varprop)
  if (is.null(numpc))
  {
    numpc=sum(varprop>0.01) #
  }
  print(paste0("number of PCs: ",numpc))
  print(paste0("proportion of variance explained: ",round(sum(varprop[1:numpc]),digits = 3)))
  res=data.frame(matrix(NA,nrow=numpc,ncol=ncol(dat)))
  
  colnames(res)=c(colnames(dat))
  
  rownames(res)=paste0("pc",1:numpc)
  
  res[,1:ncol(res)]=t(pc$x)[1:numpc,]
  res=as.data.frame(t(res))
  return(res)
  
}
tumor_meth_pca=pca(check=T,numpc = 20)
normal_meth_pca=pca(dat=normal_meth,check=T,numpc = 20)

save(tumor_meth,tumor_meth_pca,tumor_pheno,normal_meth,normal_meth_pca,normal_pheno,file="../result/tumor_normal_meth.RData")
plot(tumor_meth_pca$pc1,tumor_meth_pca$pc2,col=tumor_pheno$batch)
plot(normal_meth_pca$pc1,normal_meth_pca$pc2,col=tumor_pheno$batch)


# check for batch effect
plot_2pca=function(dat1=meth_batch1,dat2=meth_batch2)
{
  
  tmp=intersect(rownames(dat1),rownames(dat2))
  
  dat1=dat1[match(tmp,rownames(dat1)),]
  
  dat2=dat2[match(tmp,rownames(dat2)),]
  
  dat=cbind(dat1,dat2)
  
  pcadat=prcomp(t(dat),scale = F)
  
  plot(pcadat$sdev^2/sum(pcadat$sdev^2)*100,ylab="Variance explained (%)",xlab="PC",cex.axis=1.3,cex.lab=1.3)
  
  pcadat=pcadat$x
  
  mycol=c(rep("blue",ncol(dat1)),rep("green",ncol(dat2)))
  
  plot(pcadat[,1],pcadat[,2],ylim=c(min(pcadat[,2]),1.4*max(pcadat[,2])),col=mycol,xlab="PC1",ylab="PC2",cex.axis=1.2,cex.lab=1.2)
  
  legend("topleft",legend=c("Batch1","Batch2"),col=c("blue","green"),pch=1,cex=1.2)
  
}
tumor_meth_batch1=tumor_meth[,colnames(tumor_meth) %in% colnames(meth1_tumor)]
tumor_meth_batch2=tumor_meth[,colnames(tumor_meth) %in% colnames(meth2_tumor)]
normal_meth_batch1=normal_meth[,colnames(normal_meth) %in% colnames(meth1_normal)]
normal_meth_batch2=normal_meth[,colnames(normal_meth) %in% colnames(meth2_normal)]
meth_batch1=cbind(tumor_meth_batch1,normal_meth_batch1)
meth_batch2=cbind(tumor_meth_batch2,normal_meth_batch2)
plot_2pca()
plot_2pca(dat1=tumor_meth_batch1,dat2=tumor_meth_batch2)
plot_2pca(dat1=normal_meth_batch1,dat2=normal_meth_batch2)

#batch effect correction
allpheno=rbind(tumor_pheno,normal_pheno)
allpheno$type="normal"
allpheno$type[1:nrow(tumor_pheno)]="tumor"
allpheno$type=as.factor(allpheno$type)
allpheno$NCI_ID[1:ncol(tumor_meth)]=paste0(allpheno$NCI_ID[1:ncol(tumor_meth)],"_T")
allpheno$Tree=factor(allpheno$Tree,levels = c(1,2,3),ordered = F)
allmeth=cbind(tumor_meth,normal_meth)
colnames(allmeth)[1:ncol(tumor_meth)]=paste0(colnames(allmeth)[1:ncol(tumor_meth)],"_T")
all(colnames(allmeth)==allpheno$NCI_ID)



allmvalue=beta2m(allmeth)
library(sva)

#you need to combine outcomes info from batch1 and batch2 samples into outcome dataframe.

#make sure sample orders in outcome and mvalue are the same 

modcombat<-model.matrix(~type+Tree, data=allpheno) 

#batch is the batch variable ("batch1" or "batch2")  

allmvalue.combat= ComBat(dat=allmvalue, batch=allpheno$batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

allmeth.combat=m2beta(allmvalue.combat)

allmeth.combat_batch1=allmeth.combat[,colnames(allmeth.combat) %in% allpheno$NCI_ID[allpheno$batch=="batch1"]]
allmeth.combat_batch2=allmeth.combat[,colnames(allmeth.combat) %in% allpheno$NCI_ID[allpheno$batch=="batch2"]]
plot_2pca(dat1=allmeth.combat_batch1,dat2=allmeth.combat_batch2)

allmeth.combat_tumor_batch1=allmeth.combat[,colnames(allmeth.combat) %in% allpheno$NCI_ID[allpheno$batch=="batch1" & allpheno$type=="tumor"]]
allmeth.combat_tumor_batch2=allmeth.combat[,colnames(allmeth.combat) %in% allpheno$NCI_ID[allpheno$batch=="batch2" & allpheno$type=="tumor"]]
plot_2pca(dat1=allmeth.combat_tumor_batch1,dat2=allmeth.combat_tumor_batch2)

allmeth.combat_normal_batch1=allmeth.combat[,colnames(allmeth.combat) %in% allpheno$NCI_ID[allpheno$batch=="batch1" & allpheno$type=="normal"]]
allmeth.combat_normal_batch2=allmeth.combat[,colnames(allmeth.combat) %in% allpheno$NCI_ID[allpheno$batch=="batch2" & allpheno$type=="normal"]]
plot_2pca(dat1=allmeth.combat_normal_batch1,dat2=allmeth.combat_normal_batch2)

library(Rtsne)
tsne_plot=function(dat=allmeth.combat,seed=1000,perplexity=5,pheno=allpheno)
{
  pcadat=prcomp(t(dat),scale = F)
  pcadat=pcadat$x
  set.seed(seed)
  tsnedat = Rtsne(pcadat, pca = F,perplexity = perplexity,theta = 0)
  tsnedat=tsnedat$Y
  rownames(tsnedat)=pheno$NCI_ID
  xmin=min(tsnedat[,1])
  xmax=max(tsnedat[,1])
  
  xmin=xmin-0.4*(xmax-xmin)
  plot(tsnedat[,1],tsnedat[,2],col=pheno$type,pch=as.numeric(pheno$batch),cex.axis=1.3,cex.lab=1.3,xlim=c(xmin,xmax),xlab="TSNE1",ylab="TSNE2")
  
  legend("topleft",legend=c("normal_1","normal_2","tumor_1","tumor_2"),col=c(1,1,2,2),pch=c(1,2,1,2),cex=1.2)
  
}
tsne_plot(dat=allmeth.combat,seed=2000,perplexity=3,pheno=allpheno)
tsne_plot(dat=allmeth,seed=1000,perplexity=3,pheno=allpheno)

pca_plot=function(dat=allmeth.combat,pheno=allpheno,optscale=F)
{
  pcadat=prcomp(t(dat),scale = optscale)
  pcadat=pcadat$x
  if (any(rownames(pcadat)!=pheno$NCI_ID)) stop("wrong order")
  xmin=min(pcadat[,1])
  xmax=max(pcadat[,1])
  
  xmin=xmin-0.4*(xmax-xmin)
  plot(pcadat[,1],pcadat[,2],col=pheno$type,pch=as.numeric(pheno$batch),cex.axis=1.3,cex.lab=1.3,xlim=c(xmin,xmax),xlab="PC1",ylab="PC2")
  
  legend("topleft",legend=c("normal_1","normal_2","tumor_1","tumor_2"),col=c(1,1,2,2),pch=c(1,2,1,2),cex=1.2)
  return(pcadat)
  
}

tmp=pca_plot(dat=allmeth,pheno=allpheno)
all(allpheno$NCI_ID==rownames(tmp))
allpheno=cbind(allpheno,tmp[,1:20])
allpheno$Tree=factor(allpheno$Tree,levels = c(1,2,3),ordered = T)

tmp=pca_plot(dat=allmeth.combat,pheno=allpheno)
all(allpheno$NCI_ID==rownames(tmp))
allpheno_combat=cbind(allpheno[,1:6],tmp[,1:20])

tumor_pheno=allpheno[allpheno$type=="tumor",]
tumor_meth=allmeth[,colnames(allmeth) %in% tumor_pheno$NCI_ID]
all(colnames(tumor_meth)==tumor_pheno$NCI_ID)
normal_pheno=allpheno[allpheno$type=="normal",]
normal_meth=allmeth[,colnames(allmeth) %in% normal_pheno$NCI_ID]
all(colnames(normal_meth)==normal_pheno$NCI_ID)
idx=match(gsub("_T","",colnames(tumor_meth)),colnames(normal_meth))
normal_meth=normal_meth[,idx]
all(colnames(normal_meth) == gsub("_T","",colnames(tumor_meth)))
all(rownames(normal_meth) %in%rownames(tumor_meth))
normal_pheno=normal_pheno[match(colnames(normal_meth),normal_pheno$NCI_ID),]
all(colnames(normal_meth)==normal_pheno$NCI_ID)
tmp=pca_plot(dat=tumor_meth,pheno=tumor_pheno)
tumor_pheno=cbind(tumor_pheno[,1:6],tmp[,1:20])
tmp=pca_plot(dat=normal_meth,pheno=normal_pheno)
normal_pheno=cbind(normal_pheno[,1:6],tmp[,1:20])
diff_meth=tumor_meth-normal_meth
tmp=pca_plot(dat=diff_meth,pheno=tumor_pheno)
diff_pheno=cbind(tumor_pheno[,1:6],tmp[,1:20])
tumor_meth_combat=allmeth.combat[,colnames(allmeth.combat) %in% tumor_pheno$NCI_ID]
normal_meth_combat=allmeth.combat[,colnames(allmeth.combat) %in% normal_pheno$NCI_ID]
all(colnames(normal_meth_combat)==normal_pheno$NCI_ID)
#colnames(normal_meth_combat)=paste0(colnames(normal_meth_combat),"_T")
idx=match(gsub("_T","",colnames(tumor_meth_combat)),colnames(normal_meth_combat))
normal_meth_combat=normal_meth_combat[,idx]
all(colnames(normal_meth_combat) == gsub("_T","",colnames(tumor_meth_combat)))
all(rownames(normal_meth_combat) %in%rownames(tumor_meth_combat))
tmp=pca_plot(dat=tumor_meth_combat,pheno=tumor_pheno)
tumor_combat_pheno=cbind(tumor_pheno[,1:6],tmp[,1:20])
tmp=pca_plot(dat=normal_meth_combat,pheno=normal_pheno)
normal_combat_pheno=cbind(normal_pheno[,1:6],tmp[,1:20])
diff_meth_combat=tumor_meth_combat-normal_meth_combat
tmp=pca_plot(dat=diff_meth_combat,pheno=tumor_pheno)
diff_pheno_combat=cbind(tumor_pheno,tmp[,1:20])
save(diff_pheno_combat,diff_meth_combat,diff_pheno,diff_meth,tumor_meth,tumor_pheno,tumor_combat_pheno,normal_meth,normal_pheno,normal_combat_pheno,allmeth,allmeth.combat,allpheno,allpheno_combat,file="../result/combat_tumor_normal.RData")

