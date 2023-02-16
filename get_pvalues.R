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
res_lr=get_pvalues(opt="tumor_lr")
res_normal_lr=get_pvalues(opt="normal_lr")
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
qqplot(res_lr$P,title="")
qqplot(res_normal_lr$P,title="")

plot(-log10(res_normal$P),-log10(res$P),xlab="Normal -log10(P)",ylab="Tumor -log10(P)",cex.lab=1.2,cex.axis=1.2)
plot(res_normal$Beta,res$Beta,xlab="Normal effect size",ylab="Tumor effect size", cex.lab=1.2,cex.axis=1.2)
