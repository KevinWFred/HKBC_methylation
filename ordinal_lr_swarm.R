#!/usr/bin/env Rscript

load("../result/combat_tumor_normal.RData")

idx=seq(1,nrow(allmeth),10000)
if (max(idx)<nrow(allmeth)) idx=c(idx,nrow(allmeth)+1)
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

cmd=data.frame(rep("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code/ordinal_lr.R ",(length(idx)-1)),type="diff",startidx=NA,endidx=NA,outfile=NA)
for (i in 1:(length(idx)-1))
{
  cmd$startidx[i]=idx[i]
  cmd$endidx[i]=idx[i+1]
  cmd$outfile[i]=paste0("/data/BB_Bioinformatics/Kevin/HKBC_methylation/result/split/diff_",i,".txt")
}
write.table(cmd,file="ordinal_lr_diff.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#swarm -f /data/BB_Bioinformatics/Kevin/HKBC_methylation/code/ordinal_lr_diff.swarm -g 16 --module R/4.2.0 --time=4:00:00 --gres=lscratch:64

cmd=data.frame(rep("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code/ordinal_lr.R ",(length(idx)-1)),type="diffM",startidx=NA,endidx=NA,outfile=NA)
for (i in 1:(length(idx)-1))
{
  cmd$startidx[i]=idx[i]
  cmd$endidx[i]=idx[i+1]
  cmd$outfile[i]=paste0("/data/BB_Bioinformatics/Kevin/HKBC_methylation/result/split/diffM_",i,".txt")
}
write.table(cmd,file="ordinal_lr_diffm.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#swarm -f /data/BB_Bioinformatics/Kevin/HKBC_methylation/code/ordinal_lr_diffm.swarm -g 16 --module R/4.2.0 --time=4:00:00 --gres=lscratch:64
