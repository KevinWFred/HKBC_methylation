#!/usr/bin/env Rscript
library(MASS)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) stop("Four inputs are required: start indices of meth, tumor or normal, and output file name")

setwd("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code")
#load("../result/tumor_normal_meth.RData")
load("../result/combat_tumor_normal.RData") #batch effect removal 
library(lumi)
allmethM=beta2m(allmeth.combat)
tumor_methM=allmethM[,colnames(allmethM) %in% allpheno$NCI_ID[allpheno$type=="tumor"]]
tumor_pheno=allpheno[allpheno$type=="tumor",]
all(colnames(tumor_methM)==tumor_pheno$NCI_ID )
normal_methM=allmethM[,colnames(allmethM) %in% allpheno$NCI_ID[allpheno$type=="normal"]]
normal_pheno=allpheno[allpheno$type=="normal",]
all(colnames(normal_methM)==normal_pheno$NCI_ID )

opt=as.character(args[1]) #tumor or normal
if (opt=="tumor")
{
  meth=tumor_methM
  pheno=tumor_pheno
}else
{
  meth=normal_methM
  pheno=normal_pheno
}
idxstart=as.integer(args[2])
idxend=as.integer(args[3])
meth=meth[idxstart:(idxend-1),]
out=data.frame(Probe=rownames(meth),P=rep(NA,nrow(meth)),Beta=NA,SE=NA)
rownames(out)=rownames(meth)
#i=which(rownames(meth)=="cg13911857")
#i=which(rownames(meth)=="cg06491116")
#i=which(rownames(meth)=="cg26391777")
#i=which(rownames(meth)=="cg03094935")

for (i in 1:nrow(meth))
{
  if (i %%1000==0) cat(i,"..")
  
  dat=cbind.data.frame(pheno,meth=unlist(meth[i,]))
  #boxplot(dat$meth~dat$Tree)
  #boxplot(dat$meth~dat$Tree,xlab="Tree",ylab="M-value",main=rownames(meth)[i])
  # dat1=dat[dat$batch=="batch1",]
  # boxplot(dat1$meth~dat1$Tree,xlab="Tree",ylab="M-value",main=rownames(meth)[i])
  # dat2=dat[dat$batch=="batch2",]
  # boxplot(dat2$meth~dat2$Tree,xlab="Tree",ylab="M-value",main=rownames(meth)[i])
  # idx=which(abs(dat$meth-median(dat$meth))>12*mad(dat$meth))
  # if (length(idx)<3)
  #  dat$meth[idx]=NA
  m1 = tryCatch(
    expr = {
      polr(formula = Tree ~ meth+AGE+PC1+PC2, data = dat, Hess = TRUE)
    },
    error = function(e){ 
      return(NULL)
    }
  )
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
outfile=as.character(args[4])
write.table(out,file=outfile,row.names = F,sep="\t",quote=F)
print("done")