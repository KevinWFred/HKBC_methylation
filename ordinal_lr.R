#!/usr/bin/env Rscript
library(MASS)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) stop("Four inputs are required: start indices of meth, tumor or normal, and output file name")

setwd("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code")

rm_outlier=function(dat)
{
  for (j in 1:3)
  {
    idx=which(dat$Tree==j)
    #idx1=which(abs(dat$meth[idx]-median(dat$meth[idx],na.rm=T))>6*mad(dat$meth[idx],na.rm = T))
    #dat$meth[idx[idx1]]=NA
    tmp=boxplot.stats(dat$meth[idx])$out
    dat$meth[idx][dat$meth[idx] %in% tmp]=NA
  }
  return(dat)
}
rm_dat=function(dat)
{
  tmp1=max(abs(dat$meth[dat$Tree==1]),na.rm=T)
  tmp2=max(abs(dat$meth[dat$Tree==3]),na.rm=T)
  tmp3=max(abs(dat$meth[dat$Tree==3]),na.rm=T)
  if(tmp1<0.01 | tmp2<0.02|tmp3<0.01)
  {
    dat=NULL
  }
  return(dat)
}
#load("../result/tumor_normal_meth.RData")
load("../result/combat_tumor_normal.RData") #batch effect removal 
library(lumi)
allmethM=beta2m(allmeth.combat)
tumor_methM=allmethM[,colnames(allmethM) %in% allpheno$NCI_ID[allpheno$type=="tumor"]]
#tumor_pheno=allpheno[allpheno$type=="tumor",]
all(colnames(tumor_methM)==tumor_pheno$NCI_ID )
normal_methM=allmethM[,colnames(allmethM) %in% allpheno$NCI_ID[allpheno$type=="normal"]]
#normal_pheno=allpheno[allpheno$type=="normal",]
all(colnames(normal_methM)==normal_pheno$NCI_ID )
tumor_meth_combat=allmeth.combat[,colnames(allmeth.combat) %in% tumor_pheno$NCI_ID]
all(colnames(tumor_meth_combat)==tumor_pheno$NCI_ID)
normal_meth_combat=allmeth.combat[,colnames(allmeth.combat) %in% normal_pheno$NCI_ID]
all(colnames(normal_meth_combat)==normal_pheno$NCI_ID)
all(colnames(normal_meth_combat) == gsub("_T","",colnames(tumor_meth_combat)))
all(rownames(normal_meth_combat) == rownames(tumor_meth_combat))
opt=as.character(args[1]) #tumor or normal
if (opt=="tumor")
{
  meth=tumor_methM
  pheno=tumor_combat_pheno
}
if (opt=="normal")
{
  meth=normal_methM
  pheno=normal_combat_pheno
}
if (opt == "diff")
{
  meth=diff_meth_combat
  pheno=diff_pheno_combat
}
if (opt == "diffM")
{
  meth=tumor_methM-normal_methM
  pheno=diff_pheno_combat
}
idxstart=as.integer(args[2])
idxend=as.integer(args[3])
meth=meth[idxstart:(idxend-1),]
out=data.frame(Probe=rownames(meth),P=rep(NA,nrow(meth)),Beta=NA,SE=NA,sd=NA,rm=F)
rownames(out)=rownames(meth)
#i=which(rownames(meth)=="cg13911857")
#i=which(rownames(meth)=="cg06491116")
#i=which(rownames(meth)=="cg26391777")
#i=which(rownames(meth)=="cg03094935")

for (i in 1:nrow(meth))
{
  if (i %%1000==0) cat(i,"..")
  
  dat=cbind.data.frame(pheno,meth=unlist(meth[i,]))
  #dat=rm_outlier(dat=dat)
  out$sd[i]=sd(dat$meth,na.rm=T)
  #dat=rm_dat(dat=dat)
  if(!is.null(dat))
  {
    
    if (sd(dat$meth,na.rm=T)>0.01)
    {
      dat$meth=scale(dat$meth)
      #boxplot(dat$meth~dat$Tree)
      #boxplot(dat$meth~dat$Tree,xlab="Tree",ylab="M-value",main=rownames(meth)[i])
      if (!opt %in% c("diff","diffM"))
      {
        m1 = tryCatch(
          expr = {
            polr(formula = Tree ~ meth+AGE+batch+PC1+PC2, data = dat, Hess = TRUE)
          },
          error = function(e){ 
            return(NULL)
          },
          warning=function(w) {
            return(NULL)
          }
        )
      }
      else
      {
        m1 = tryCatch(
          expr = {
            polr(formula = Tree ~ meth+AGE+batch, data = dat, Hess = TRUE)
          },
          error = function(e){ 
            return(NULL)
          },
          warning=function(w) {
            return(NULL)
          }
        )
      }
      if(!is.null(m1))
      {
        ctable = tryCatch(
          expr= {coef(summary(m1))
          },
          error = function(e){ 
            return(NULL)
          }
        )
        if (!is.null(ctable))
        {
          ctable = cbind(ctable, "p_value" = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2)
          out[i, "SE"] = ctable["meth", "Std. Error"]
          out[i, "P"] = ctable["meth", "p_value"]
          out[i, "Beta"] = ctable["meth", "Value"]
        }
      }
    }
  }else
  {
    out$rm[i]=T
  }
  
}
# qqplot(out$P)
outfile=as.character(args[4])
write.table(out,file=outfile,row.names = F,sep="\t",quote=F)
print("done")