#!/usr/bin/env Rscript
library(MASS)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) stop("Four inputs are required: start indices of meth, tumor or normal, and output file name")

setwd("/data/BB_Bioinformatics/Kevin/HKBC_methylation/code")
load("../result/tumor_normal_meth.RData")
library(lumi)
tumor_methM=beta2m(tumor_meth)
normal_methM=beta2m(normal_meth)
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
for (i in 1:nrow(meth))
{
  if (i %%1000==0) cat(i,"..")
  dat=cbind.data.frame(pheno,meth=unlist(meth[i,]))
  #boxplot(dat$meth~dat$Tree)
  # idx=which(abs(dat$meth-median(dat$meth))>12*mad(dat$meth))
  # if (length(idx)<3)
  #  dat$meth[idx]=NA
  m1 = tryCatch(
    expr = {
      polr(formula = Tree ~ meth+AGE+batch, data = dat, Hess = TRUE)
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
