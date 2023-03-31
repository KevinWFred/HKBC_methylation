#!/usr/bin/env Rscript
setwd("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code")
qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL,title="")
  
{
  
  pvalue=pvalue[!is.na(pvalue)]
  n=length(pvalue)
  par(mar=c(5,5,2,1))
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log base 10)",
       
       ylab="Observed p-value (log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.2,cex.axis=1.2)
  title(main=title,cex=1.2)
  abline(0,1,lty=2)
  chisq <- qchisq(1-pvalue,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
  
}

load("../result/combat_tumor_normal.RData") #batch effect removal 
library(lumi)
library(limma)
library(sva)
allmethM=beta2m(allmeth.combat)
allpheno$id=allpheno$NCI_ID
allpheno$id=as.factor(gsub("_T","",allpheno$id))
designmat <- model.matrix(~type+AGE+batch+PC1+PC2+PC3+PC4,data=allpheno)
fit <- lmFit(as.matrix(allmethM), designmat)
# n.sv = num.sv(as.matrix(allmethM),designmat)
# designmat1 <- model.matrix(~1,data=allpheno)
# svobj = sva(as.matrix(allmethM),designmat,designmat1,n.sv=n.sv)
# modSv = as.data.frame(cbind(designmat,svobj$sv))
# fit <- lmFit(as.matrix(allmethM), modSv)
result=eBayes(fit)
colnames(result$p.value)[2]
#qqplot(result$p.value[,2],title="Compare tumor vs normal")
out=as.data.frame(result$p.value)
out$fdr=p.adjust(out[,2],method = "BH")
out$fwer=p.adjust(out[,2],method = "bonferroni")
all(row.names(out)==rownames(tumor_meth))
all(row.names(out)==rownames(normal_meth))
out$mean_normal=rowMeans(normal_meth)
out$mean_tumor=rowMeans(tumor_meth)
idx=match(rownames(out),epicanno$Name)
out$gene=epicanno$UCSC_RefGene_Name[idx]
save(out,file="../result/diff_tumor_normal.RData")
gene="REC8"
gene="DNM3"
idx=which(grepl(gene,out$gene))
out1=out[idx,]
out1$meandiff=out1$mean_normal-out1$mean_tumor
out1=out1[order(out1$meandiff),]
cpgs=rownames(out1)[out1$fdr<0.05]
par(mfrow=c(2,3))
for (i in 1:6)
{
  cpg=cpgs[i]
  idx=which(rownames(allmeth)==cpg)
  boxplot(unlist(allmeth[idx,])~allpheno$type,main=paste0(gene,":",cpg),ylab="Beta value",xlab="")
}
dim(out1)
out2=out1[out1$fdr<0.05,]
tmp=c(out2$mean_tumor,out2$mean_normal)
par(mfrow=c(1,1))
boxplot(tmp~c(rep("tumor",nrow(out2)),rep("normal",nrow(out2))),ylab="Mean beta value",xlab="",main=gene)
normal_cutoff1=0.2
normal_cutoff2=0.8
hypercpgs=rownames(out)[which(out$fdr<0.05& out$mean_tumor>out$mean_normal & out$mean_normal<normal_cutoff1)]
hypercpgs=rownames(out)[which(out$fdr<0.05& out$mean_tumor>out$mean_normal)]
length(hypercpgs) #54670
hypocpgs=rownames(out)[which(out$fdr<0.05& out$mean_tumor<out$mean_normal)]
length(hypocpgs) #70400

#bumphunter
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(minfi)
epicanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

epicanno$chr <- as.character(epicanno$chr)
epicanno$chr <- gsub("chr","",epicanno$chr)  
epicanno <- epicanno[epicanno$chr!="X"&epicanno$chr!="Y",]
epicanno$chr <- as.numeric(epicanno$chr)
epicanno$pos <- as.numeric(epicanno$pos)
epicanno <- epicanno[row.names(epicanno) %in% row.names(allmeth),]
epicanno <- epicanno[order(epicanno$chr,epicanno$pos),]
epicanno=as.data.frame(epicanno)
bumphunterfun=function(designMat= designmat,dat,coef=2)
{
  library(bumphunter)
  library(doParallel)
  registerDoParallel(cores = 1)
  set.seed(39283)
  idx=match(rownames(dat),epicanno$Name)
  dmrs <- bumphunter(dat,designMat,chr=epicanno$chr[idx],pos=epicanno$pos[idx],nullMethod="bootstrap",cutoff=0.1,coef=coef,B=100,type="beta")
  return(dmrs)
}

# dat=dat[match(epicanno$Name,rownames(dat)),]
# designmat <- model.matrix(~type+AGE,data=allpheno)
# dmrs1 <- bumphunter(as.matrix(allmeth),designmat,chr=epicanno$chr,pos=epicanno$pos,cutoff=0.1,B=500,type="beta")
# designmat <- model.matrix(~type+AGE+PC1+PC2+PC3+PC4+PC5,data=allpheno)
# dmrs1_5 <- bumphunter(as.matrix(allmeth),designmat,chr=epicanno$chr,pos=epicanno$pos,cutoff=0.1,B=500,type="beta")

#on tumor/normal
tumor_meth_combat=allmeth.combat[,colnames(allmeth.combat) %in% tumor_pheno$NCI_ID]
all(colnames(tumor_meth_combat)==tumor_pheno$NCI_ID)
normal_meth_combat=allmeth.combat[,colnames(allmeth.combat) %in% normal_pheno$NCI_ID]
all(colnames(normal_meth_combat)==normal_pheno$NCI_ID)
all(colnames(normal_meth_combat) == gsub("_T","",colnames(tumor_meth_combat)))
all(rownames(normal_meth_combat) == rownames(tumor_meth_combat))

#all cpgs
designmat <- model.matrix(~I(Tree==1)+AGE,data=tumor_pheno)
dmrs_tumor_tree1=bumphunterfun(designMat = designmat,dat=tumor_meth_combat)
designmat <- model.matrix(~I(Tree==1)+AGE,data=normal_pheno)
dmrs_normal_tree1=bumphunterfun(designMat = designmat,dat=normal_meth_combat)
designmat <- model.matrix(~I(Tree==2)+AGE,data=tumor_pheno)
dmrs_tumor_tree2=bumphunterfun(designMat = designmat,dat=tumor_meth_combat)
designmat <- model.matrix(~I(Tree==2)+AGE,data=normal_pheno)
dmrs_normal_tree2=bumphunterfun(designMat = designmat,dat=normal_meth_combat)
designmat <- model.matrix(~I(Tree==3)+AGE,data=tumor_pheno)
dmrs_tumor_tree3=bumphunterfun(designMat = designmat,dat=tumor_meth_combat)
designmat <- model.matrix(~I(Tree==3)+AGE,data=normal_pheno)
dmrs_normal_tree3=bumphunterfun(designMat = designmat,dat=normal_meth_combat)
designmat <- model.matrix(~I(Tree==3)+AGE+PC1+PC2,data=normal_pheno)
dmrs_normal_tree3PC2=bumphunterfun(designMat = designmat,dat=normal_meth_combat)
designmat <- model.matrix(~I(Tree==3)+AGE,data=tumor_pheno[tumor_pheno$Tree %in% c(1,3),])
dmrs_tumor_tree3_1=bumphunterfun(designMat = designmat,dat=tumor_meth_combat[,tumor_pheno$Tree %in% c(1,3)])
designmat <- model.matrix(~I(Tree==3)+AGE,data=tumor_pheno[tumor_pheno$Tree %in% c(2,3),])
dmrs_tumor_tree3_2=bumphunterfun(designMat = designmat,dat=tumor_meth_combat[,tumor_pheno$Tree %in% c(2,3)])
designmat <- model.matrix(~I(Tree==1)+AGE,data=tumor_pheno[tumor_pheno$Tree %in% c(1,2),])
dmrs_tumor_tree1_2=bumphunterfun(designMat = designmat,dat=tumor_meth_combat[,tumor_pheno$Tree %in% c(1,2)])
designmat <- model.matrix(~I(Tree==3)+AGE,data=normal_pheno[normal_pheno$Tree %in% c(1,3),])
dmrs_normal_tree3_1=bumphunterfun(designMat = designmat,dat=normal_meth_combat[,normal_pheno$Tree %in% c(1,3)])
designmat <- model.matrix(~I(Tree==3)+AGE,data=normal_pheno[normal_pheno$Tree %in% c(2,3),])
dmrs_normal_tree3_2=bumphunterfun(designMat = designmat,dat=normal_meth_combat[,normal_pheno$Tree %in% c(2,3)])
designmat <- model.matrix(~I(Tree==1)+AGE,data=normal_pheno[normal_pheno$Tree %in% c(1,2),])
dmrs_normal_tree1_2=bumphunterfun(designMat = designmat,dat=normal_meth_combat[,normal_pheno$Tree %in% c(1,2)])

#hyperhyop
designmat <- model.matrix(~I(Tree==1)+AGE,data=tumor_pheno)
dat=tumor_meth_combat[match(c(hypercpgs,hypocpgs),rownames(tumor_meth_combat)),]
dmrs_hyperhypo_tumor_tree1=bumphunterfun(designMat = designmat,dat=dat)
designmat <- model.matrix(~I(Tree==3)+AGE,data=tumor_pheno)
dat=tumor_meth_combat[match(c(hypercpgs,hypocpgs),rownames(tumor_meth_combat)),]
dmrs_hyperhypo_tumor_tree3=bumphunterfun(designMat = designmat,dat=dat)

designmat <- model.matrix(~I(Tree==1)+AGE,data=normal_pheno)
dat=normal_meth_combat[match(c(hypercpgs,hypocpgs),rownames(normal_meth_combat)),]
dmrs_hyperhypo_normal_tree1=bumphunterfun(designMat = designmat,dat=dat)
designmat <- model.matrix(~I(Tree==3)+AGE,data=normal_pheno)
dat=normal_meth_combat[match(c(hypercpgs,hypocpgs),rownames(normal_meth_combat)),]
dmrs_hyperhypo_normal_tree3=bumphunterfun(designMat = designmat,dat=dat)

#hyper
designmat <- model.matrix(~I(Tree==1)+AGE,data=tumor_pheno)
dat=tumor_meth_combat[match(c(hypercpgs),rownames(tumor_meth_combat)),]
dmrs_hyper_tumor_tree1=bumphunterfun(designMat = designmat,dat=dat)
designmat <- model.matrix(~I(Tree==3)+AGE,data=tumor_pheno)
dat=tumor_meth_combat[match(c(hypercpgs),rownames(tumor_meth_combat)),]
dmrs_hyper_tumor_tree3=bumphunterfun(designMat = designmat,dat=dat)

designmat <- model.matrix(~I(Tree==1)+AGE,data=normal_pheno)
dat=normal_meth_combat[match(c(hypercpgs),rownames(normal_meth_combat)),]
dmrs_hyper_normal_tree1=bumphunterfun(designMat = designmat,dat=dat)
designmat <- model.matrix(~I(Tree==3)+AGE,data=normal_pheno)
dat=normal_meth_combat[match(c(hypercpgs),rownames(normal_meth_combat)),]
dmrs_hyper_normal_tree3=bumphunterfun(designMat = designmat,dat=dat)

#hypo
designmat <- model.matrix(~I(Tree==1)+AGE,data=tumor_pheno)
dat=tumor_meth_combat[match(c(hypocpgs),rownames(tumor_meth_combat)),]
dmrs_hypo_tumor_tree1=bumphunterfun(designMat = designmat,dat=dat)
designmat <- model.matrix(~I(Tree==3)+AGE,data=tumor_pheno)
dat=tumor_meth_combat[match(c(hypocpgs),rownames(tumor_meth_combat)),]
dmrs_hypo_tumor_tree3=bumphunterfun(designMat = designmat,dat=dat)

designmat <- model.matrix(~I(Tree==1)+AGE,data=normal_pheno)
dat=normal_meth_combat[match(c(hypocpgs),rownames(normal_meth_combat)),]
dmrs_hypo_normal_tree1=bumphunterfun(designMat = designmat,dat=dat)
designmat <- model.matrix(~I(Tree==3)+AGE,data=normal_pheno)
dat=normal_meth_combat[match(c(hypocpgs),rownames(normal_meth_combat)),]
dmrs_hypo_normal_tree3=bumphunterfun(designMat = designmat,dat=dat)

#save.image(file="../result/bumphunterres.RData")

#read bumphunter res
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
      # cpgs=anno$Name[dat1$indexStart[i]:dat1$indexEnd[i]]
      # for (j in 1:length(cpgs))
      # {
      #   idx=which(rownames(epicmeth)==cpgs[j])
      #   boxplot(epicmeth[idx,]~sampletable1$type)
      # }
      tmp=unique(anno$UCSC_RefGene_Name[dat1$indexStart[i]:dat1$indexEnd[i]])
      tmp=tmp[tmp!=""]
      tmp=paste0(tmp,collapse = ",")
      dat1$genename[i]=tmp
      
      #generate cpg list
      cpgs=anno$Name[dat1$indexStart[i]:dat1$indexEnd[i]]
      tmp=annotate_cpg(cpgs)
      tmp$DMR=paste0("DMR",i)
      tmp$pvalue=dat1$p.valueArea[i]
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
    
    return(list(dat1=dat1,dat4=dat4))
  }
}
bumphunter_tumor_tree1=annotate_bumphunter(dat=dmrs_tumor_tree1)
bumphunter_tumor_tree2=annotate_bumphunter(dat=dmrs_tumor_tree2)#
bumphunter_tumor_tree3=annotate_bumphunter(dat=dmrs_tumor_tree3)
bumphunter_normal_tree1=annotate_bumphunter(dat=dmrs_normal_tree1)
bumphunter_normal_tree2=annotate_bumphunter(dat=dmrs_normal_tree2)
bumphunter_normal_tree3=annotate_bumphunter(dat=dmrs_normal_tree3,plotflag = T) #
bumphunter_tumor_tree3_1=annotate_bumphunter(dat=dmrs_tumor_tree3_1)
bumphunter_tumor_tree3_2=annotate_bumphunter(dat=dmrs_tumor_tree3_2)
bumphunter_tumor_tree1_2=annotate_bumphunter(dat=dmrs_tumor_tree1_2)
bumphunter_normal_tree3_1=annotate_bumphunter(dat=dmrs_normal_tree3_1)
bumphunter_normal_tree3_2=annotate_bumphunter(dat=dmrs_normal_tree3_2)
bumphunter_normal_tree1_2=annotate_bumphunter(dat=dmrs_normal_tree1_2)
# View(dmrs_hyper_normal_tree3$table)
# dat=bumphunter_normal_tree3$dat1[,c(1:3,9,11:15)]
# dat=cbind(DMR=paste0("DMR",1:nrow(dat)),dat)
# View(dat)
# 
# qqplot(dmrs_hyper_normal_tree3$table$p.valueArea)
bumphunter_hyperhypo_tumor_tree1=annotate_bumphunter(dat=dmrs_hyperhypo_tumor_tree1)
bumphunter_hyperhypo_tumor_tree3=annotate_bumphunter(dat=dmrs_hyperhypo_tumor_tree3)
bumphunter_hyperhypo_normal_tree1=annotate_bumphunter(dat=dmrs_hyperhypo_normal_tree1)
bumphunter_hyperhypo_normal_tree3=annotate_bumphunter(dat=dmrs_hyperhypo_normal_tree3) #

bumphunter_hyper_tumor_tree1=annotate_bumphunter(dat=dmrs_hyper_tumor_tree1)
bumphunter_hyper_tumor_tree3=annotate_bumphunter(dat=dmrs_hyper_tumor_tree3)
bumphunter_hyper_normal_tree1=annotate_bumphunter(dat=dmrs_hyper_normal_tree1)
bumphunter_hyper_normal_tree3=annotate_bumphunter(dat=dmrs_hyper_normal_tree3) #

bumphunter_hypo_tumor_tree1=annotate_bumphunter(dat=dmrs_hypo_tumor_tree1)
bumphunter_hypo_tumor_tree3=annotate_bumphunter(dat=dmrs_hypo_tumor_tree3)
bumphunter_hypo_normal_tree1=annotate_bumphunter(dat=dmrs_hypo_normal_tree1)
bumphunter_hypo_normal_tree3=annotate_bumphunter(dat=dmrs_hypo_normal_tree3) #

view_dmr=function(bhres=bumphunter_tumor_tree1_2)
{
  dat=bhres$dat1[,c(1:3,9,11:15)]
  dat=cbind(DMR=paste0("DMR",1:nrow(dat)),dat)
  View(dat)
  return(dat)
}
view_dmr()
view_dmr4=function(bhres=bumphunter_normal_tree1_3)
{
  dat=bhres$dat4
  View(dat)
  return(dat)
}
tmp=view_dmr4(bhres=bumphunter_normal_tree3)
tmp1=tmp[tmp$DMR=="DMR1",1:3]
write.table(tmp1,file="../result/bumphunter_normal_tree3_REC8.csv",row.names=F,sep=",",quote=F)
tmp1=tmp[tmp$DMR=="DMR3",1:3]
write.table(tmp1,file="../result/bumphunter_normal_tree3_DNM3.csv",row.names=F,sep=",",quote=F)
write.table(tmp,file="../result/bumphunter_normal_tree3.csv",row.names=F,sep=",",quote=F)

#visualize DMR
source("./generate_DMRplot.R")
library(gtable)
plot_DMR=function(Data=normal_meth_combat,Data.group=normal_pheno$Tree,bhres=bumphunter_normal_tree3,ymax=1.1,outprefix="Normal_tree3")
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
    groups=as.character(unique(Data.group))
    title=unique(bhres$dat4$genename[which(bhres$dat4$DMR==dmrs[i])])
    tmp=Plot.gene(Data, Data.group, probe.annotation,plot.title=title,groups = groups,ylim.max = ymax)
    print(tmp$Line)
    ggsave(paste0("../result/",outprefix,"DMR",i,"_line.pdf"),width=12,height=8)
    print(tmp$Box)
    ggsave(paste0("../result/",outprefix,"DMR",i,"_box.pdf"),width=12,height=8)
    plotres[i]=tmp
  }
  return(plotres)
}
allplotres1=plot_DMR(Data=normal_meth)

#DMR on another data
allplotres1=plot_DMR(Data=tumor_meth_combat,Data.group=tumor_pheno$Tree,outprefix="Normal_tree3_tumor")

plot_gene=function(cpgs,groups=c("1","2","3"),title="",ymax=0.8)
  
{
  
  annotation=epicanno
  annotation=annotation[annotation$Name %in% rownames(Data),]
  
  idx=match(cpgs,annotation$Name)
  probe.annotation=annotation[idx,]
  
  Plot.gene(Data, Data.group, probe.annotation,plot.title=title,groups = groups,ylim.max = ymax)
  
}

#sbatch --mem=6umphunter_normal_tree34g --cpus-per-task=12 --time=3-12:00:00 --gres=lscratch:64  /data/BB_Bioinformatics/Kevin/HKBC_methylation/code/diff_analysis.R