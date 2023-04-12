#!/usr/bin/env Rscript

setwd("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code")

get_pvalues=function(opt="tumor")
{
  res=NULL
  for (i in 1:78)
  {
    tmp=read.table(paste0("../result/split/",opt,"_",i,".txt"),header=T)
    res=rbind(res,tmp)
  }
  #res=res[!is.na(res$P),]
  return(res)
}

res=get_pvalues()
res_normal=get_pvalues(opt="normal")
res_diff=get_pvalues(opt="diff")
res_diffm=get_pvalues(opt="diffM")
res_lr=get_pvalues(opt="tumor_lr")
res_normal_lr=get_pvalues(opt="normal_lr")
res_lrPC=get_pvalues(opt="tumorPC_lr")
res_normal_lrPC=get_pvalues(opt="normalPC_lr")
res_lrtnoPC=get_pvalues(opt="lrtnoPC")
res_lrtPC=get_pvalues(opt="lrtPC")

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
png("../result/tumor_probes_qqplot.png",res=100)
qqplot(res$P,title="Tumor")
dev.off()
qqplot(res_normal$P,title="Normal")
qqplot(res_normal$P,ylim=c(0,10))
qqplot(res_lr$P,title="Logistic regression, Tree3 vs others, tumor")
qqplot(res_normal_lr$P,title="Logistic regression, Tree3 vs others, normal")
qqplot(res_lrPC$P,title="Logistic regression, Tree3 vs others, tumor, adjust for 3 PCs")
qqplot(res_normal_lrPC$P,title="Logistic regression, Tree3 vs others, normal, adjust for 3 PCs")
png("../result/diff_probes_qqplot.png",res=100)
qqplot(res_diff$P,title="Tumor - normal")
dev.off()
qqplot(res_diffm$P,title="Tumor - normal, M value")
qqplot(res_diff$P,ylim=c(0,5),title="Tumor - normal")
qqplot(res_lrtnoPC$P,title="Interaction")
qqplot(res_lrtnoPC$P,ylim=c(0,10))
sum(res_lrtnoPC$P<0.05/nrow(res_lrtnoPC),na.rm=T)
qqplot(res_lrtPC$P,title="Interaction, 4PCs")
qqplot(res_lrtPC$P,ylim=c(0,10))
plot(-log10(res_normal$P),-log10(res$P),xlab="Normal -log10(P)",ylab="Tumor -log10(P)",cex.lab=1.2,cex.axis=1.2)
plot(res_normal$Beta,res$Beta,xlab="Normal effect size",ylab="Tumor effect size", cex.lab=1.2,cex.axis=1.2)
cpg="cg08271447"
cpg="cg09187087"
cpg="cg24579224"
cpg="cg11613644"
dat=cbind.data.frame(allpheno,meth=unlist(allmeth[which(rownames(allmeth)==cpg),]))
dat1=cbind.data.frame(tumor_pheno,meth=unlist(tumor_meth[which(rownames(tumor_meth)==cpg),]))
dat2=cbind.data.frame(normal_pheno,meth=unlist(normal_meth[which(rownames(normal_meth)==cpg),]))
boxplot(dat$meth~dat$type)
boxplot(dat1$meth~dat1$Tree,main=paste0(cpg," tumor"),xlab="Tree",ylab="Beta")
boxplot(dat2$meth~dat2$Tree,main=paste0(cpg," normal"),xlab="Tree",ylab="Beta")

#manhattan plot
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
source("./generate_Manhattanplot.R")
idx=match(res$Probe,epicanno$Name)
png("../result/ordlr_tumor.png",width=960)
draw_manhattan(pvalues=data.frame(chr=epicanno$chr[idx],start=epicanno$pos[idx],P=res$P),maxy=max(-log10(res$P),na.rm=T))
dev.off()
idx=match(res_normal$Probe,epicanno$Name)
png("../result/ordlr_normal.png",width=960)
draw_manhattan(pvalues=data.frame(chr=epicanno$chr[idx],start=epicanno$pos[idx],P=res_normal$P),maxy=max(-log10(res_normal$P),na.rm=T))
dev.off()
#pathway analysis
library(missMethyl)
# #OPT can be "GO" or "KEGG"
# generich=function(dat=res,opt="GO")
# {
#   dat$P.adj=p.adjust(dat$P,method="BH")
#   selectprobes=dat$Probe[which(dat$P<0.005)]
# 
#   par(mfrow=c(1,1))
#   gst <- gometh(sig.cpg=selectprobes, collection=opt,sig.genes = T, array.type = c("EPIC"))
#   print(topGSA(gst, number=10))
#   return(gst)
# }

# 
# tumor_gst_GO=generich()
# tumor_gst_KEGG=generich(opt="KEGG")
# normal_gst_GO=generich(dat=res_normal)
# normal_gst_KEGG=generich(dat=res_normal,opt="KEGG")
# save(tumor_gst_GO,tumor_gst_KEGG,normal_gst_GO,normal_gst_KEGG,file="../result/genrichres.RData")

#for tumor-normal comparison
load("../result/diff_tumor_normal.RData")
generich=function(dat=out,opt="GO",opt1="all")
{
  dat$P.adj=p.adjust(dat$typetumor,method="bonferroni")
  if (opt1=="all")
  {
    selectprobes=rownames(dat)[which(dat$P.adj<0.05)]
  }
  if (opt1=="hyper")
  {
    selectprobes=rownames(dat)[which(dat$P.adj<0.05 & dat$mean_tumor-dat$mean_normal>0.1)]
  }
  if (opt1=="hypo")
  {
    selectprobes=rownames(dat)[which(dat$P.adj<0.05 & dat$mean_normal-dat$mean_tumor>0.1)]
  }
  
  par(mfrow=c(1,1))
  gst <- gometh(sig.cpg=selectprobes, collection=opt,sig.genes = T, array.type = c("EPIC"))
  print(topGSA(gst, number=10))
  return(gst)
}

diff_gst_GO=generich()
tmp=cbind.data.frame(GO=rownames(diff_gst_GO),diff_gst_GO)
tmp=tmp[order(tmp$P.DE),]
tmp=tmp[tmp$FDR<0.05,]

write.csv(tmp,file="../result/diff_tumor_normal_genrich_GO.csv",row.names = F)
diff_gst_KEGG=generich(opt="KEGG")
tmp=cbind.data.frame(KEGG=rownames(diff_gst_KEGG),diff_gst_KEGG)
tmp=tmp[order(tmp$P.DE),]
tmp=tmp[tmp$P.DE<0.05,]
write.csv(tmp,file="../result/diff_tumor_normal_genrich_KEGG.csv",row.names = F)
save(diff_gst_GO,diff_gst_KEGG,file="../result/diff_tumor_normal_genrichres.RData")

diff_gst_GO_hyper=generich(opt1="hyper")
tmp=cbind.data.frame(GO=rownames(diff_gst_GO_hyper),diff_gst_GO_hyper)
tmp=tmp[order(tmp$P.DE),]
tmp=tmp[tmp$FDR<0.05,]
write.csv(tmp,file="../result/diff_tumor_normal_genrich_GO_hyper.csv",row.names = F)

diff_gst_GO_hypo=generich(opt1="hypo")
diff_gst_KEGG_hyper=generich(opt="KEGG",opt1="hyper")
tmp=cbind.data.frame(KEGG=rownames(diff_gst_KEGG_hyper),diff_gst_KEGG_hyper)
tmp=tmp[order(tmp$P.DE),]
tmp=tmp[tmp$P.DE<0.05,]
write.csv(tmp,file="../result/diff_tumor_normal_genrich_KEGG_hyper.csv",row.names = F)

diff_gst_KEGG_hypo=generich(opt="KEGG",opt1="hypo")
tmp=cbind.data.frame(KEGG=rownames(diff_gst_KEGG_hypo),diff_gst_KEGG_hypo)
tmp=tmp[order(tmp$P.DE),]
tmp=tmp[tmp$P.DE<0.05,]
write.csv(tmp,file="../result/diff_tumor_normal_genrich_KEGG_hypo.csv",row.names = F)

save(diff_gst_GO,diff_gst_KEGG,diff_gst_GO_hyper,diff_gst_GO_hypo,
     diff_gst_KEGG_hyper,diff_gst_KEGG_hypo,file="../result/diff_tumor_normal_genrichres.RData")

