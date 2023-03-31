#!/usr/bin/env Rscript

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

annotate_cpg=function(cpg="cg25305703")
{
  anno1=epicanno
  anno1$chr=factor(anno1$chr,levels = c(1:22,"X","Y"))
  idx=order(anno1$chr,anno1$pos)
  anno1=anno1[idx,]
  idx=match(cpg,anno1$Name)
  res=data.frame(cpg=cpg,chr=anno1$chr[idx],position=anno1$pos[idx],genegroup=as.character(anno1$UCSC_RefGene_Group[idx]),
                 promoter=F,genename=NA,stringsAsFactors = F)
  for (i in 1:length(idx))
  {
    if (grepl("Promoter",anno1$Regulatory_Feature_Group[idx[i]])) res$promoter[i]=T
    # if (!is.na(anno1$Enhancer[idx[i]]))
    # {
    #   if(anno1$Enhancer[idx[i]]==T) res$enhancer[i]=T
    # }
    
    if (anno1$UCSC_RefGene_Name[idx[i]]!="" & !is.na(anno1$UCSC_RefGene_Name[idx[i]]))
    {
      res$genename[i]=anno1$UCSC_RefGene_Name[idx[i]]
    }else #intergenic
    {
      tmp=unique(anno1$UCSC_RefGene_Name[(idx[i]-min(500,idx[i]-1)):idx[i]])
      tmp=tmp[tmp!=""]
      tmp=unlist(strsplit(tmp,"|",fixed=T))
      tmp1=unique(anno1$UCSC_RefGene_Name[(idx[i]+1):(idx[i]+500)])
      tmp1=tmp1[tmp1!=""]
      tmp1=unlist(strsplit(tmp1,"|",fixed=T))
      res$genename[i]=paste0(tmp[length(tmp)],"/",tmp1[1])
    }
    #it has gencode annoation
    if (anno1$UCSC_RefGene_Name[idx[i]]=="" & anno1$GencodeBasicV12_NAME[idx[i]]!="") res$genename[i]=anno1$GencodeBasicV12_NAME[idx[i]]
  }
  res$genename=processsemicolon(res$genename)
  res$genegroup=processsemicolon(res$genegroup)
  return(res)
}

#sort table based on chromosome
sortgenetable=function(genetable)
{
  genetable$chr=as.character(genetable$chr)
  genetable$chr=gsub("chr","",genetable$chr)
  genetable$chr=gsub("23","X",genetable$chr)
  genetable$chr=gsub("24","Y",genetable$chr)
  if (class(genetable$start[1])=="factor")
  {
    genetable$start=as.numeric(as.character(genetable$start))
  }
  if (class(genetable$start[1])=="character")
  {
    genetable$start=as.numeric(genetable$start)
  }
  
  chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
  res=data.frame(matrix(NA,nrow=0,ncol=ncol(genetable)))
  for (chr in chrs)
  {
    tmptable=genetable[which(genetable$chr==chr),]
    tmptable=tmptable[order(tmptable$start),]
    res=rbind(res,tmptable)
  }
  return(res)
}

#data have 3 columns: chr, pos, pvalue
mhtplot2=function (data, control = mht.control(),  ...) 
{
  for (p in c("grid")) {
    if (length(grep(paste("^package:", p, "$", sep = ""), 
                    search())) == 0) {
      if (!require(p, quietly = TRUE, character.only = TRUE)) 
        warning(paste("mhtplot needs package `", p, "' to be fully functional; please install", 
                      sep = ""))
    }
  }
  #length of each chromosome
  chrs=c(1:22,"X","Y")
  dictfile="/data/BB_Bioinformatics/Kevin/tools/database/ucsc.hg19.compact.dict"
  chrlen=read.table(file=dictfile,sep="\t",skip=1,stringsAsFactors = F)
  chrlen=chrlen$V3
  chrlen=gsub("LN:","",chrlen,fixed=T)
  chrlen=as.numeric(chrlen)
  names(chrlen)=chrs
  
  data2 <- data[!apply(is.na(data), 1, any), ]
  idxnoNA <- !apply(is.na(data), 1, any)
  #the index of nonNA values. NAs will not be plotted
  idxnoNA <- which(idxnoNA==T)                  
  n2 <- dim(data2[1])
  chr <- as.character(data2[, 1])
  uniqchr=unique(chr)
  chrs1=chrs[chrs %in% uniqchr]
  idx=match(chrs1,uniqchr)
  uniqchr=uniqchr[idx]
  #if only draw chrom with data
  chrlen1=chrlen[names(chrlen) %in% chr]
  chrlen2=c(0,cumsum(chrlen1)) #start point of each chr
  names(chrlen2)=c(names(chrlen1),"nextchr")
  pos <- newpos <- data2[, 2]
  for (i in 1:length(uniqchr))
  {
    idx=which(data2[,1]==uniqchr[i])
    startp=chrlen2[which(names(chrlen2)==uniqchr[i])]
    # print(uniqchr[i])
    # print(range(idx))
    # print(startp)
    #the overall cooridinate start from 0
    pos[idx] = newpos[idx] =data2[idx,2]+startp
  }
  if (class(data2[1,3])=="character") data2[,3]=as.numeric(data2[,3])
  p <- data2[, 3]
  
  #this part has problem, the orders is chr1,chr10,chr2,...; should change to chr1,chr2,...,chr10,... 
  #allchr <- as.vector(tablechr)
  
  allchr=NULL
  for (i in 1:length(uniqchr))
  {
    allchr=c(allchr,sum(chr==uniqchr[i]))
  }
  n.chr <- length(allchr)
  type <- control$type
  #usepos <- control$usepos
  logscale <- control$logscale
  base <- control$base
  cutoffs <- control$cutoffs
  colors <- control$colors
  labels <- control$labels
  srt <- control$srt
  gap <- control$gap
  pcex <- control$cex
  yline <- control$yline
  xline <- control$xline
  colorlist <- colors()
  if (is.null(colors)) 
    colors <- sample(colorlist, n.chr)
  
  if (is.null(labels)) 
    labels <- uniqchr
  if (is.null(gap)) 
    gap <- 0
  
  CMindex <- cumsum(allchr)
  #the overall cooridinate start from 0
  CM=newpos
  args <- list(...)
  if ("ylim" %in% names(args)) 
    dp <- seq(args$ylim[1], args$ylim[2], length = sum(allchr))
  else dp <- seq(min(p), max(p), length = sum(allchr))
  if (logscale) 
    y <- -log(dp, base)
  else y <- dp
  y1 <- min(y)
  y2 = max(y)
  par(xaxt = "n", yaxt = "n")
  #   xy <- xy.coords(CM, y)
  #   plot(xy$x, xy$y, type = "n", ann = FALSE, axes = FALSE, ...)
  plot(x=chrlen2[c(1,length(chrlen2))],y=c(floor(y1),ceiling(y2)),type="n",ann = FALSE, axes = FALSE)
  #axis(1)
  #axis(2)
  par(xaxt = "s", yaxt = "s")
  for (i in 1:n.chr) {
    u <- CMindex[i]
    l <- CMindex[i] - allchr[i] + 1
    idx <- l:u
    #cat("Plotting points ", l, "-", u, "\n")
    if (logscale) 
      y <- -log(p[idx], base)
    else y <- p[idx]
    col.chr <- colors[i]
    if (type == "l") 
      lines(CM[idx], y, col = col.chr, cex = pcex, ...)
    else 
    {
      points(CM[idx], y, col = col.chr, cex = pcex, ...)
    }
    #text(ifelse(i == 1, CM[1], CM[l]), y1, pos = 1, offset = 1, 
    #     labels[i], srt = srt, ...)
    # if (i<=15 | (i %% 2==1 & i>15))
    # {
    #   text(1/2*(chrlen2[i]+chrlen2[i+1]), y1, pos = 1, offset = 1, labels[i], srt = srt, cex=1.1)
    # }
    text(1/2*(chrlen2[i]+chrlen2[i+1]), y1, pos = 1, offset = 1, labels[i], srt = srt, cex=0.8)
  }
  
  if (!is.null(cutoffs)) 
    segments(0, cutoffs, n2 + gap * n.chr, cutoffs)
  if ("ylab" %in% names(args)) 
    mtext(args$ylab, 2, line = yline, las = 0, cex=1.3)
  else mtext(ifelse(logscale, paste("-log", base, "(Observed value)", 
                                    sep = ""), "Observed value"), 2, line = yline, las = 0, ...)
  if ("xlab" %in% names(args)) 
    xlabel <- args$xlab
  else xlabel <- ifelse(is.null(names(chr)), "Chromosome", names(chr))
  mtext(xlabel, 1, line = xline, las = 0, cex=1.3)
  #CM:coordinate of points, for idxnoNA rows in the original data which don't have NAs
  return(res=list(CM=CM,CMindex=CMindex,chrlen2=chrlen2,y1=y1,y2=y2,idxnoNA=idxnoNA)) # x coordinate and number of points in each chr
}

library(gap)
draw_manhattan=function(pvalues=NULL,fdrs=NULL,fdrthreshold=0.05,maxy=6,chrs=NULL,keepres=F,logscale=T,main=NULL,ylab=NULL)
{
  #if fdrs is a dataframe, the first row is pvalue, second row is fdr
  if (is.null(pvalues))
  {
    pvalues=fdrs[,1]
    if(class(pvalues)=="character") pvalues=as.numeric(pvalues)
    names(pvalues)=rownames(fdrs)
  }else
  {
    if (class(pvalues)=="data.frame")
    {
      if (! "chr" %in% colnames(pvalues) | ! "start" %in% colnames(pvalues)) #if the dataframe have start and chr columns, it already have coordinates
      {
        tmp=rownames(pvalues)
        pvalues=pvalues[,ncol(pvalues)] #the last column is pvalues
        if(class(pvalues)=="character") pvalues=as.numeric(pvalues)
        names(pvalues)=tmp
      }
    }
  }
  if (class(pvalues)=="numeric") #without gene coordinate, needs to figure out
  {
    if (is.null(chrs))
    {
      chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
      #chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
    }
    geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)
    #form a dataframe with gene's position and p value
    df_pvalues=data.frame(gene=names(pvalues),pvalue=pvalues)
    genetable=merge(df_pvalues,geneposition,all.x=T)
    idxkeep=which(genetable$chr %in% chrs)
    genetable=genetable[idxkeep,]
    genetable=sortgenetable(genetable)
    genenames=genetable$gene
    #form the genetable for mahanttan plot
    genetable=data.frame(chr=genetable$chr,start=genetable$start,pvalue=genetable$pvalue,stringsAsFactors=F)
    rownames(genetable)=genenames
  }else
  {
    if (class(pvalues)=="data.frame")
    {
      warning("the columns should be chr, start and pvalue")
      genetable=pvalues
      genetable=sortgenetable(genetable)
    }
  }
  
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  par(las=1, xpd=TRUE,cex.axis=1.2,cex=1,cex.lab=1.5)
  ops<-mht.control(colors=color2,yline=1,xline=2,usepos=T,srt=0,logscale=logscale)
  if (is.null(ylab))
  {
    if (logscale==T)
    {
      ylab=expression(Observed~~-log[10](italic(p)))
    }else
    {
      ylab=expression(Observed~~italic(p))
    }
  }
  
  res=mhtplot2(genetable,ops,pch=10,bg=color2,ylab=ylab)
  #res=mhtplot2(genetable,ops,pch=10,bg=color2,cex.axis=1.2,cex.lab=1.2)
  #allpos keeps the overall cooridinate
  genetable=cbind(genetable,allpos=rep(NA,nrow(genetable)))
  genetable$allpos[res$idxnoNA]=res$CM
  #lines(c(0,max(res$CM)),c(0,0),lwd=2)
  if (logscale==T)
  {
    axis(2,at=c(floor(res$y1):(ceiling(res$y2))),pos=0, cex.axis=1.2)
  }else
  {
    axis(2,pos=0, cex.axis=1.2)
  }
  xat=c(res$chrlen2)
  axis(1,pos=floor(res$y1),labels=FALSE,tick=T,at=xat)
  if (! is.null(main))
    title(main=main)
  
  #to add fdr line
  #draw the fdr cutoff line
  if (!is.null(fdrs))
  {
    dictfile="/data/BB_Bioinformatics/Kevin/tools/database/ucsc.hg19.compact.dict"
    chrlen=read.table(file=dictfile,sep="\t",skip=1,stringsAsFactors = F)
    chrlen=chrlen$V3
    chrlen=gsub("LN:","",chrlen,fixed=T)
    chrlen=as.numeric(chrlen)
    names(chrlen)=c(1:22,"X","Y")
    #available chrs length
    chrlen1=chrlen[names(chrlen) %in% genetable$chr]
    
    #change fdrs,pvalues to vector
    fdrs1=fdrs
    if (class(fdrs)=="data.frame")
    {
      tmp=rownames(fdrs)
      fdrs1=fdrs[,ncol(fdrs)] #assume fdr is in the last column
      names(fdrs1)=tmp
    }
    pvalues1=pvalues
    if (class(pvalues)=="data.frame")
    {
      tmp=rownames(pvalues)
      pvalues1=pvalues[,ncol(pvalues)] #assume pvalue is in the last column
      names(pvalues1)=tmp
    }
    
    smallfdrs=fdrs1[fdrs1<=fdrthreshold]
    if (length(smallfdrs)>0) #if have fdr<fdrthreshold
    {
      idxs=which(fdrs1 %in% smallfdrs)
      pvaluethreshold=max(pvalues1[idxs])
      segments(0,-log10(pvaluethreshold),par('usr')[2],-log10(pvaluethreshold),col="red")
      #text(sum(chrlen1)/2,-log10(pvaluethreshold)+0.5,paste0("FDR=",fdrthreshold))
      #text(sum(chrlen1)*0.9,-log10(pvaluethreshold)+0.25,paste0("FDR=",fdrthreshold),cex=1.1)
      text(sum(chrlen1)*0.9,-log10(pvaluethreshold)+0.25,paste0("FDR=",round(max(smallfdrs),digits = 2)),cex=1.1)
    }else
    {
      warning("no genes with small fdr were found")
      #find the gene with min fdr
      
      minfdr=round(min(fdrs1,na.rm=T),digits = 2)
      idxs=which(fdrs1==minfdr)
      
      #use the one with max p-value
      thepvalue=max(pvalues1[idxs])
      #abline(h=-log10(thepvalue),col="red")
      segments(0,-log10(thepvalue),par('usr')[2],-log10(thepvalue),col="red")
      if (class(fdrs)=="data.frame")
      {
        #text(sum(chrlen1)/2,-log10(thepvalue)+0.5,paste0("FDR=",minfdr))
        text(sum(chrlen1)*0.9,-log10(thepvalue)+0.25,paste0("FDR=",minfdr),cex=1.1)
      }else
      {
        #text(sum(chrlen1)/2,-log10(thepvalue)+0.5,paste0("FDR=",minfdr))
        text(sum(chrlen1)*0.9,-log10(thepvalue)+0.25,paste0("FDR=",minfdr),cex=1.1)
      }
    }
  }
  par(cex.axis=1,cex=1,cex.lab=1)
  if (keepres==T) return(genetable)
}
