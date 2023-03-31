#!/usr/bin/env Rscript

setwd("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code")
load("../result/tumor_normal_meth.RData")

idx=seq(1,nrow(meth1_tumor),10000)
if (max(idx)<nrow(meth1_tumor)) idx=c(idx,nrow(meth1_tumor)+1)
cmd=data.frame(rep("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code/ordinal_lr.R ",(length(idx)-1)),type="tumor",startidx=NA,endidx=NA,outfile=NA)
for (i in 1:(length(idx)-1))
{
  cmd$startidx[i]=idx[i]
  cmd$endidx[i]=idx[i+1]
  cmd$outfile[i]=paste0("/data/BB_Bioinformatics/Kevin/HKBC_methylation/result/split/tumor_",i,".txt")
}
write.table(cmd,file="ordinal_lr_tumor.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#swarm -f /data/BB_Bioinformatics/Kevin/HKBC_methylation/code/ordinal_lr_tumor.swarm -g 16 --module R/4.2.0 --time=4:00:00 --gres=lscratch:64

cmd=data.frame(rep("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code/ordinal_lr.R ",(length(idx)-1)),type="normal",startidx=NA,endidx=NA,outfile=NA)
for (i in 1:(length(idx)-1))
{
  cmd$startidx[i]=idx[i]
  cmd$endidx[i]=idx[i+1]
  cmd$outfile[i]=paste0("/data/BB_Bioinformatics/Kevin/HKBC_methylation/result/split/normal_",i,".txt")
}
write.table(cmd,file="ordinal_lr_normal.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#swarm -f /data/BB_Bioinformatics/Kevin/HKBC_methylation/code/ordinal_lr_normal.swarm -g 16 --module R/4.2.0 --time=4:00:00 --gres=lscratch:64
