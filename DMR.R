#!/usr/bin/env Rscript
setwd("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code")

#datasets for 41 tumor normal pairs before and after batch effect removal 
#normal_meth: normal data before
#normal_meth_combat: normal data after 
#normal_pheno: phenotype data for normal samples
#allmeth.combat: (normal + tumor) after
#allpheno: phenotype data for allmeth
load("../result/combat_tumor_normalData1.RData") 

#annotation for epic array
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
epicanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

library(lumi)
library(limma)

#m values
#allmeth.combat is the methylation after batch effect removal
allmethM=beta2m(allmeth.combat)
allpheno$id=allpheno$NCI_ID
allpheno$id=as.factor(gsub("_T","",allpheno$id))

#get the hyper-methylated CpGs
designmat <- model.matrix(~type+AGE+batch+PC1+PC2+PC3+PC4,data=allpheno)
fit <- lmFit(as.matrix(allmethM), designmat)
result=eBayes(fit)
colnames(result$p.value)[2]
out=as.data.frame(result$p.value)
out$fdr=p.adjust(out[,2],method = "BH")
out$fwer=p.adjust(out[,2],method = "bonferroni")
out$mean_normal=rowMeans(normal_meth)
out$mean_tumor=rowMeans(tumor_meth)
idx=match(rownames(out),epicanno$Name)
out$gene=epicanno$UCSC_RefGene_Name[idx]
hypercpgs=rownames(out)[which(out$fdr<0.05& out$mean_tumor>out$mean_normal)]
length(hypercpgs) #54670

#to detect DMR associate Tree==3 in normal
epicanno$chr <- as.character(epicanno$chr)
epicanno$chr <- gsub("chr","",epicanno$chr)  
epicanno <- epicanno[epicanno$chr!="X"&epicanno$chr!="Y",]
epicanno$chr <- as.numeric(epicanno$chr)
epicanno$pos <- as.numeric(epicanno$pos)
epicanno <- epicanno[row.names(epicanno) %in% row.names(allmeth),]
epicanno <- epicanno[order(epicanno$chr,epicanno$pos),]
epicanno=as.data.frame(epicanno)
bumphunterfun=function(designMat= designmat,dat,coef=2,B=100)
{
  library(bumphunter)
  library(doParallel)
  registerDoParallel(cores = 1)
  set.seed(39283)
  idx=match(rownames(dat),epicanno$Name)
  dmrs <- bumphunter(dat,designMat,chr=epicanno$chr[idx],pos=epicanno$pos[idx],nullMethod="bootstrap",cutoff=0.1,coef=coef,B=B,type="beta")
  #DMR results are saved in dmrs$table
  idx=which(dmrs$table$p.valueArea==0)
  #to get FDR
  dmrs$table$p.valueArea[idx]=1/nrow(dmrs$table)/B
  dmrs$table$fdrArea=p.adjust(dmrs$table$p.valueArea,method = "BH")
  return(dmrs)
}
normal_pheno$Tree3=normal_pheno$Tree==3
designmat <- model.matrix(~Tree3+AGE,data=normal_pheno)
dat=normal_meth_combat[match(c(hypercpgs),rownames(normal_meth_combat)),]
dmrs_hyper_normal_tree3=bumphunterfun(designMat = designmat,dat=dat)

#to annotate DMRs
processsemicolon=function(txt="5'UTR;1stExon;1stExon;1stExon;5'UTR;5'UTR;5'UTR;1stExon;5'UTR")
{
  txt=as.character(txt)
  res=txt
  for (i in 1:length(txt))
  {
    if (grepl(";",txt[i]))
    {
      tmp=unlist(strsplit(txt[i],";"))
      res[i]=paste0(unique(tmp),collapse="|")
    }
  }
  
  return(res)
}

#the following code has some functions used for bumphnter results annotation
source("./generate_Manhattanplot.R")
epicanno$UCSC_RefGene_Name=processsemicolon(epicanno$UCSC_RefGene_Name)

annotate_bumphunter=function(dat=dmrs_tumor_tree1,plotflag=F)
{
  idx0=which(dat$table$fwerArea<0.05)
  if (length(idx0)>0)
  {
    dat1=dat$table[idx0,]
    idx=match(rownames(dat$coef),epicanno$Name)
    anno=epicanno[idx,]
    dat1$genename=NA
    #cpg list
    dat4=NULL
    for (i in 1:nrow(dat1))
    {
      tmp=unique(anno$UCSC_RefGene_Name[dat1$indexStart[i]:dat1$indexEnd[i]])
      tmp=tmp[tmp!=""]
      tmp=paste0(tmp,collapse = ",")
      dat1$genename[i]=tmp
      
      #generate cpg list
      cpgs=anno$Name[dat1$indexStart[i]:dat1$indexEnd[i]]
      tmp=annotate_cpg(cpgs)
      tmp$DMR=paste0("DMR",i)
      tmp$pvalue=dat1$p.valueArea[i]
      tmp$fdr=dat1$fdrArea[i]
      dat4=rbind(dat4,tmp)
    }
    if (plotflag==T)
    {
      #udpate the pvalues and fwers
      dat2=data.frame(chr=rep(NA,length(dat$coef)),start=NA,pvalue=NA,fwer=NA,stringsAsFactors = F)
      rownames(dat2)=rownames(dat$coef)
      for (i in 1:nrow(dat$table))
      {
        cpgs=anno$Name[dat$table$indexStart[i]:dat$table$indexEnd[i]]
        idx=match(cpgs,anno$Name)
        idx1=match(cpgs,rownames(dat2))
        dat2$chr[idx1]=anno$chr[idx]
        dat2$start[idx1]=anno$pos[idx]
        dat2$pvalue[idx1]=dat$table$p.valueArea[i]
        dat2$fwer[idx1]=dat$table$fwerArea[i]
      }
      dat3=dat2[!is.na(dat2$chr),1:3]
      par(mar=c(4,4,2,1))
      draw_manhattan(pvalues=dat3,maxy=max(-log10(dat3$pvalue)))
      dat3=dat2[!is.na(dat2$chr),]
      idx=order(dat3$pvalue)
      dat3=dat3[idx,]
      idx=which(dat3$fwer<0.05)
      tmp=-log10(dat3$pvalue[length(idx)])
      abline(h=tmp,col="red")
    }
    #dat1:DMR results
    #dat4:CpGs within DMRs
    return(list(dat1=dat1,dat4=dat4))
  }
}
bumphunter_hyper_normal_tree3=annotate_bumphunter(dat=dmrs_hyper_normal_tree3) #

#get DMRs in supplementary table 2
tmp=bumphunter_hyper_normal_tree3$dat4
tmp$meanBeta_tree1=tmp$meanBeta_tree2=tmp$meanBeta_tree3=NA
for (i in 1:nrow(tmp))
{
  idx=which(rownames(normal_meth)==tmp$cpg[i])
  idx1=which(normal_pheno$Tree==1)
  tmp$meanBeta_tree1[i]=round(mean(unlist(normal_meth[idx,idx1]),na.rm=T),2)
  idx2=which(normal_pheno$Tree==2)
  tmp$meanBeta_tree2[i]=round(mean(unlist(normal_meth[idx,idx2]),na.rm=T),2)
  idx3=which(normal_pheno$Tree==3)
  tmp$meanBeta_tree3[i]=round(mean(unlist(normal_meth[idx,idx3]),na.rm=T),2)
}
#write.table(tmp,file="../result/bumphunter_normal_tree3_.csv",row.names=F,sep=",",quote=F)

#visualize DMR Figure 5c and supFig 11
source("./generate_DMRplot.R")
library(gtable)
plot_DMR=function(Data=normal_meth_combat,Data.group=normal_pheno$Tree,bhres=bumphunter_hyper_normal_tree3,ymax=1.1,outprefix="Normal_tree3")
{
  plotres=list()
  dmrs=unique(bhres$dat4$DMR)
  for (i in 1:length(dmrs))
  {
    cpgs=bhres$dat4$cpg[which(bhres$dat4$DMR==dmrs[i])]
    annotation=epicanno
    annotation=annotation[annotation$Name %in% rownames(Data),]
    idx=match(cpgs,annotation$Name)
    probe.annotation=annotation[idx,]
    #groups=as.character(unique(Data.group))
    groups=c("Tumour Only","Shared Tree","Multiple trees")
    levels(Data.group)=c("Tumour Only","Shared Tree","Multiple trees") #the legend position GIMP
    title=unique(bhres$dat4$genename[which(bhres$dat4$DMR==dmrs[i])])
    tmp=Plot.gene(Data, Data.group, probe.annotation,plot.title=title,groups = groups,ylim.max = ymax)
    print(tmp$Line)
    #ggsave(paste0("../result/",outprefix,"DMR",i,"_line.pdf"),width=12,height=8)
    print(tmp$Box)
    #ggsave(paste0("../result/",outprefix,"DMR",i,"_box.pdf"),width=11,height=6)
    plotres[i]=tmp
  }
  return(plotres)
}

allDMRplot=plot_DMR(Data=normal_meth)
