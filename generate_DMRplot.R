library(ggplot2)

qqplot <- function(pvalue = NULL,main = "qqplot", xlim=NULL, ylim=NULL)
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (-log10)",
       ylab="Observed p-value (-log10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.2,cex.axis=1.2)
  abline(0,1,lty=2)
}

plot.probe <- function(Data.j, Data.group, annotation.j,
                       types = c("Beta"),
                       title.color = "black", 
                       p.value = NULL, qvalue = NULL){
  
  plot.df <- data.frame(Beta = as.numeric(Data.j), 
                        Group = Data.group)
  
  plot.probe.list <- vector("list", length(types))
  for (i in c(1:length(types))) {
    p.tmp <- NULL
    if(types[i] == "Beta"){
      p.tmp <- ggplot(plot.df, aes(Group, Beta))
    }else if(types[i] == "M"){
      p.tmp <- ggplot(plot.df, aes(Group, M))
    }else if(types[i] == "Abs.median.diff"){
      p.tmp <- ggplot(plot.df, aes(Group, Median.diff.abs))
    }
    
    p.title <- paste(rownames(annotation.j), " - ", 
                     annotation.j$genename, " : ", annotation.j$chr, " - ", annotation.j$pos,
                     '\n', annotation.j$Regulatory_Feature_Group, sep='')
    if(!is.null(p.value)){
      p.title <- paste(p.title, '\n', "p-value ", signif(p.value, digits = 3), sep='')
    }
    if(!is.null(qvalue)){
      p.title <- paste(p.title, '\n', "qvalue ", signif(qvalue, digits = 3), sep='')
    }
    
    p.tmp <- p.tmp +
      geom_boxplot(colour = "red", outlier.shape=NA) +
      geom_jitter() +
      ggtitle(p.title) +
      theme(
        plot.title = element_text(color = title.color)
      )
    
    plot.probe.list[[i]] <- p.tmp
  }
  return(plot.probe.list)
}

#process UCSC_RefGene_Group info, make it unique for a probe
process_group=function(myannotation=probe.annotation)
{
  idx=which(grepl(";",myannotation$UCSC_RefGene_Group))
  if (length(idx)>0)
  {
    myannotation0=myannotation
    idx=order(myannotation$chr,myannotation$pos)
    myannotation=myannotation[idx,]
    idx=which(grepl(";",myannotation$UCSC_RefGene_Group))
    idx0=which(!grepl(";",myannotation$UCSC_RefGene_Group))
    allgroups=c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR")
    for (i in idx)
    {
      tmp=unique(unlist(strsplit(myannotation$UCSC_RefGene_Group[i],';')))
      if (length(tmp)>0)
      {
        idx1=match(tmp,allgroups)
        myannotation$UCSC_RefGene_Group[i]=paste0(tmp[order(idx1)],collapse = ";")
      }
    }
  }
  idx=which(grepl(";",myannotation$UCSC_RefGene_Group))
  idx0=which(!grepl(";",myannotation$UCSC_RefGene_Group))
  if (length(idx)>0)
  {
    for (i in idx)
    {
      tmp=unique(unlist(strsplit(myannotation$UCSC_RefGene_Group[i],';')))
      idx1=match(tmp,allgroups)
      groupstart=tmp[order(idx1)]
      groupend=groupstart[length(groupstart)]
      groupstart=groupstart[1]
      idxbefore=which(!is.na(myannotation$UCSC_RefGene_Group))
      idxbefore=idxbefore[idxbefore<i]
      idxbefore=idxbefore[length(idxbefore)]
      groupbefore=unique(myannotation$UCSC_RefGene_Group[idxbefore])
      idx2=match(groupbefore,allgroups)
      idx3=idx1[idx1>=idx2]
      #groupbefore=groupbefore[order(idx1)]
      # groupbefore=groupbefore[length(groupbefore)]
      # idxafter=idx0[idx0>i]
      # groupafter=unique(myannotation$UCSC_RefGene_Group[idxafter])
      # idx1=match(groupafter,allgroups)
      # groupafter=groupafter[order(idx1)]
      # groupafter=groupafter[1]
      # if (match(groupend,allgroups)>match(groupafter,allgroups) & 
      #     match(groupstart,allgroups)<=match(groupafter,allgroups))
      # {
      #   myannotation$UCSC_RefGene_Group[i]=groupstart
      # }
      # if (match(groupstart,allgroups)<match(groupbefore,allgroups) & 
      #     match(groupend,allgroups)<=match(groupafter,allgroups))
      # {
      #   myannotation$UCSC_RefGene_Group[i]=groupend
      # }
      myannotation$UCSC_RefGene_Group[i]=allgroups[idx3[1]]
    }
    idx=match(rownames(myannotation0),rownames(myannotation))
    myannotation=myannotation[idx,]
  }
  
  return(myannotation)
}
# Plot probes of a gene.
# Default known variables:
#   'Data', 'Data.group'
# Parameters:
#   'probe.annotation': probe annotation
#   'genome': hg18, hg19, hg38, mm9, mm10
#   'additional.probe.list': list of types of DMRs, each type is a list of DMRs
#   'additional.track.level': track in the plot for the types of DMRs
#                             numerical values starting from 1
#                             same length as 'additional.probe.list'
Plot.gene <- function(Data, Data.group, probe.annotation, groups=NULL,ylim.max=1.1,
                      additional.probe.list = NULL, 
                      additional.track.level = NULL,
                      additional.line = NULL,
                      plot.title = NULL){
  ylim.min <- if(is.null(additional.probe.list)) -0.2 else -0.25 - length(additional.track.level) * 0.05
  #ylim.min <- if(is.null(additional.probe.list)) -0.09 else -0.11 - length(additional.track.level) * 0.05
  
  if(!is.matrix(Data)) Data <- data.matrix(Data)
  
  if(!is.null(additional.probe.list)){
    additional.probe.list <- lapply(additional.probe.list, function(x){
      if(is.null(names(x))) names(x) <- c(1:length(x))
      return(x)
    })
  }
  
  p.output.list <- list(Line = NULL, Box = NULL)
  
  probe.annotation=process_group(probe.annotation)
  annotation.gene.list.i <- probe.annotation
  annotation.gene.list.i <- annotation.gene.list.i[rownames(annotation.gene.list.i) %in% rownames(Data),]
  annotation.gene.list.i <- annotation.gene.list.i[order(annotation.gene.list.i$pos),]
  
  
  # Check promoter region
  promoter.crd <- unique(annotation.gene.list.i[annotation.gene.list.i$Regulatory_Feature_Group %in% c("Promoter_Associated", "Promoter_Associated_Cell_type_specific"),]$Regulatory_Feature_Name)
  promoter.crd <- promoter.crd[promoter.crd != ""]
  promoter.plot.df <- NULL
  if(!is.null(promoter.crd) && length(promoter.crd) != 0){
    promoter.crd.tmp <- lapply(promoter.crd, function(x){
      return(as.numeric(strsplit(x, split = ":|-", fixed = FALSE)[[1]]))
    })
    
    promoter.plot.df <- do.call("rbind", lapply(promoter.crd.tmp, function(x){
      data.frame(xmin = x[2], xmax = x[3], ymin = -0.1, ymax = -0.05, Region.type = "Promoter region", label = "")
      #data.frame(xmin = x[2], xmax = x[3], ymin = -0.05, ymax = -0.03, Region.type = "Promoter region", label = "")
    }))
  }
  
  if (!is.null(promoter.plot.df))
  {
    promoter.plot.df$xmax[nrow(promoter.plot.df)]=min(max(annotation.gene.list.i$pos)+500,promoter.plot.df$xmax[nrow(promoter.plot.df)])
  }
  ####
  
  # Check CpG island region
  island.crd <- unique(annotation.gene.list.i$Islands_Name)
  island.crd <- island.crd[island.crd != ""]
  island.plot.df <- NULL
  if(!is.null(island.crd) && length(island.crd) != 0){
    island.crd.tmp <- lapply(island.crd, function(x){
      return(as.numeric(strsplit(x, split = ":|-", fixed = FALSE)[[1]][2:3]))
    })
    
    island.plot.df <- do.call("rbind", lapply(island.crd.tmp, function(x){
      data.frame(xmin = x[1], xmax = x[2], ymin = -0.15, ymax = -0.1, Region.type = "CpG island", label = "")
      #data.frame(xmin = x[1], xmax = x[2], ymin = -0.07, ymax = -0.05, Region.type = "CpG island", label = "")
    }))
  }
  if (!is.null(island.plot.df))
  {
    island.plot.df$xmax[nrow(island.plot.df)]=min(max(annotation.gene.list.i$pos)+500,island.plot.df$xmax[nrow(island.plot.df)])
  }
  ####
  
  # # Gene regions from annotation
  # gene.df <- NULL
  # annotation.gene.list.i.genes <- unique(annotation.gene.list.i$genename)
  # annotation.gene.list.i.genes <- annotation.gene.list.i.genes[annotation.gene.list.i.genes != ""]
  # if(length(annotation.gene.list.i.genes) != 0){
  #   gene.df <- do.call("rbind", lapply(annotation.gene.list.i.genes, function(j){
  #     annotation.gene.list.i.j <- annotation.gene.list.i[annotation.gene.list.i$genename == j,]
  #     gene.df.tmp <- data.frame(xmin = min(annotation.gene.list.i.j$pos), 
  #                               xmax = max(annotation.gene.list.i.j$pos), 
  #                               ymin = -0.2, 
  #                               ymax = -0.15,
  #                               Region.type = "Gene overlapping",
  #                               label = "")
  #     return(gene.df.tmp)
  #   }))
  # }
  # ####
  
  # Group from annotation
  group.df <- NULL
  annotation.gene.list.i.groups <- unique(annotation.gene.list.i$UCSC_RefGene_Group)
  annotation.gene.list.i.groups <- annotation.gene.list.i.groups[!annotation.gene.list.i.groups %in% c("",NA)]
  if(length(annotation.gene.list.i.groups) != 0){
    group.df <- do.call("rbind", lapply(annotation.gene.list.i.groups, function(j){
      annotation.gene.list.i.j <- annotation.gene.list.i[which(annotation.gene.list.i$UCSC_RefGene_Group == j),]
      group.df.tmp <- data.frame(xmin = min(annotation.gene.list.i.j$pos),
                                 xmax = max(annotation.gene.list.i.j$pos),
                                 ymin = -0.2,
                                 ymax = -0.15,
                                 Region.type = j,
                                 label = "")
      # group.df.tmp <- data.frame(xmin = min(annotation.gene.list.i.j$pos), 
      #                            xmax = max(annotation.gene.list.i.j$pos), 
      #                            ymin = -0.09, 
      #                            ymax = -0.07,
      #                            Region.type = j,
      #                            label = "")
      return(group.df.tmp)
    }))
  }
  ####
  
  # Additional regions
  additional.df <- NULL
  if(!is.null(additional.probe.list)){
    additional.df <- do.call("rbind", lapply(c(1:length(additional.probe.list)), function(j){
      if(length(additional.probe.list[[j]]) == 0) return(NULL)
      do.call("rbind", lapply(c(1:length(additional.probe.list[[j]])), function(l){
        additional.df.tmp <- annotation.gene.list.i[additional.probe.list[[j]][[l]],]$pos
        additional.df.tmp <- data.frame(xmin = min(additional.df.tmp), 
                                        xmax = max(additional.df.tmp), 
                                        ymin = -0.2 - additional.track.level[j] * 0.05, 
                                        ymax = -0.2 - (additional.track.level[j]-1) * 0.05,
                                        Region.type = names(additional.probe.list)[j],
                                        label = names(additional.probe.list[[j]])[l])
      }))
    }))
    additional.df$Region.type <- factor(additional.df$Region.type, levels = names(additional.probe.list))
  }
  ####
  
  # Combine the regions into one dataframe
  #Add.region.df <- rbind(promoter.plot.df, island.plot.df, group.df, additional.df)
  Add.region.df <- rbind(promoter.plot.df, island.plot.df, group.df, additional.df)
  if(!is.null(Add.region.df)){
    Add.region.df$Region.type <- factor(Add.region.df$Region.type,
                                        levels = c("Promoter region", 
                                                   "CpG island", 
                                                   "TSS1500",
                                                   "TSS200",
                                                   "5'UTR",
                                                   "1stExon",
                                                   "Body",
                                                   "3'UTR",
                                                   levels(additional.df$Region.type)))
    Add.region.df$Region.type=droplevels(Add.region.df$Region.type)
    
    Add.region.df.missing.type <- levels(Add.region.df$Region.type)[! levels(Add.region.df$Region.type) %in% Add.region.df$Region.type]
    if(length(Add.region.df.missing.type) != 0){
      Add.region.df <- rbind(Add.region.df,
                             data.frame(xmin = Add.region.df$xmin[1],
                                        xmax = Add.region.df$xmax[1],
                                        ymin = 100,
                                        ymax = 101,
                                        Region.type = Add.region.df.missing.type, 
                                        label = ""))
    }
  }
  
  ####
  
  Data.i <- Data[rownames(annotation.gene.list.i),]
  
  plot.df <- data.frame(Beta = c(Data.i),
                        Sample = rep(colnames(Data.i), each = nrow(Data.i)),
                        Probe = rep(rownames(Data.i), ncol(Data.i)),
                        Pos = rep(as.numeric(annotation.gene.list.i$pos), ncol(Data.i)),
                        Strand = rep(as.character(annotation.gene.list.i$strand), ncol(Data.i)),
                        Group = rep(Data.group, each = nrow(Data.i)))
  #CpG islands are not overlaped with CpGs, remove them
  idx=which(Add.region.df$Region.type=="CpG island")
  if (length(idx)>0)
  {
    rmidx=NULL
    for (j in idx)
    if (Add.region.df$xmin[j]>max(as.integer(as.character(plot.df$Pos))) | Add.region.df$xmax[j]< min(as.integer(as.character(plot.df$Pos))))
    {
      rmidx=c(rmidx,j)
    }
    if (!is.null(rmidx))
      Add.region.df=Add.region.df[-rmidx,]
  }
  
  if (!is.null(groups))
  {
    plot.df$Group=factor(plot.df$Group,levels = groups)
  }
  
  plot.df$Probe <- factor(plot.df$Probe, levels = rownames(Data.i))
  
  Data.group.unique=unique(Data.group)
  plot.generegion=data.frame(Beta=rep(NA,nrow(Data.i)*length(Data.group.unique)),
                             Probe=rep(rownames(Data.i),length(Data.group.unique)),
                             Pos=rep(as.numeric(annotation.gene.list.i$pos),length(Data.group.unique)),
                             Group=rep(Data.group.unique,each=nrow(Data.i)))
  if (!is.null(groups))
  {
    plot.generegion$Group=factor(plot.generegion$Group,levels = groups)
  }
  for (i in 1:nrow(plot.generegion))
  {
    plot.generegion$Beta[i]=mean(Data.i[which(rownames(Data.i)==plot.generegion$Probe[i]),which(Data.group==plot.generegion$Group[i])],na.rm=T)
  }
  #fit <- aov(Beta ~ Group + Pos, data = plot.df)
  plot.generegion$Group=factor(plot.generegion$Group,levels=as.character(unique(plot.generegion$Group)))
  # Methylation scatterplot
  allcolors=c("blue","red","green","purple","brown","black")
  p1 <- ggplot() + 
    theme_bw() +
    theme(text = element_text(size=14),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.direction="horizontal",
          legend.position = "top",
          legend.title = element_blank(),legend.text = element_text(size=10),axis.title.x = element_text(margin = margin(t = 15))) +
    geom_point(data = plot.generegion, aes(x=Pos, y=Beta, color=Group),size=1.5) + labs(x="Chromosome coordinate",y="Mean of beta") +
    scale_color_manual(values = allcolors[1:length(Data.group.unique)])+
    coord_cartesian(ylim = c(ylim.min, ylim.max))
  for (y in rev(levels(plot.generegion$Group))) {
    p1 <- p1 +
      geom_line(data = plot.generegion[plot.generegion$Group == y,], aes(x=Pos, y=Beta, color=Group, group = Group), alpha = 0.4,size=1.5)
  }
  leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box") 
  if(!is.null(Add.region.df)){
    p2=geom_rect(data = Add.region.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Region.type), 
                 color = "white")
    #leg2 <- gtable_filter(ggplot_gtable(ggplot_build(p2)), "guide-box") 
    # p1 <- p1 + geom_rect(data = Add.region.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Region.type), 
    #                      color = "white")
  }
  p1 <- p1+
    scale_y_continuous(breaks=seq(0,1,0.1),labels=c(seq(0,1,0.1)))+
    #geom_text(data = plot.df, aes(x=Pos, label = Strand), y=1, size = 10) +
    scale_shape(solid = FALSE) +
    #ylim(ylim.min, 1) +
    #geom_hline(yintercept = 0) +
    # geom_hline(yintercept = -0.1, linetype = "dashed") +
    ggtitle(paste(plot.title,", ", "chr",probe.annotation$chr[1],", ", length(unique(plot.df$Probe)), " probes", sep='')) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(), 
          legend.position = "top")
  
  x1=min(as.numeric(as.character(plot.df$Pos)))
  x2=max(as.numeric(as.character(plot.df$Pos)))
  
  #top legend for groups
  plotNew <- p1 +
    scale_x_continuous(breaks=seq(x1,x2,length.out=3)) +
    annotation_custom(grob = leg1, xmin = x1, xmax = x2, ymin = 0.85*ylim.max, ymax = 0.95*ylim.max)
  p=plotNew+p2+theme(text = element_text(size=14), axis.text = element_text(face="bold"), plot.title = element_text(hjust = 0.5),
                     legend.title = element_blank(), 
                     legend.position = "bottom")+guides(color = FALSE)+guides(fill = guide_legend(nrow = 1))
  
  # plotNew <- arrangeGrob(leg2, plotNew,
  #                        heights = unit.c(leg1$height, unit(1, "npc") -  leg1$height), ncol = 1)
  
  # Change coordinates into factors
  plot.df$Pos <- factor(plot.df$Pos, levels = sort(unique(plot.df$Pos)))
  
  Add.region.factor.df <- NULL
  if(!is.null(Add.region.df)){
    Add.region.factor.df <- do.call("rbind", lapply(c(1:nrow(Add.region.df)), function(x){
      pos.levels <- as.numeric(levels(plot.df$Pos))
      position.tmp <- which(pos.levels >= Add.region.df$xmin[x] & pos.levels <= Add.region.df$xmax[x])
      if(length(position.tmp) == 0){
        return(NULL)
      }else{
        return(data.frame(xmin = position.tmp[1]-0.35, 
                          xmax = tail(position.tmp, 1)+0.35, 
                          ymin = Add.region.df$ymin[x], 
                          ymax = Add.region.df$ymax[x],
                          Region.type = Add.region.df$Region.type[x],
                          label = Add.region.df$label[x]))
      }
    }))
    
    Add.region.factor.df.missing.type <- levels(Add.region.factor.df$Region.type)[! levels(Add.region.factor.df$Region.type) %in% Add.region.factor.df$Region.type]
    if(length(Add.region.factor.df.missing.type) != 0){
      Add.region.factor.df <- rbind(Add.region.factor.df,
                                    data.frame(xmin = Add.region.factor.df$xmin[1],
                                               xmax = Add.region.factor.df$xmax[1],
                                               ymin = 100,
                                               ymax = 101,
                                               Region.type = Add.region.factor.df.missing.type,
                                               label = ""))
    }
  }
  
  additional.labels <- NULL
  if(!all(Add.region.factor.df$label == "")){
    additional.labels <- Add.region.factor.df[Add.region.factor.df$label != "",]
    additional.labels$xmid <- (additional.labels$xmin + additional.labels$xmax) / 2
    additional.labels$ymid <- (additional.labels$ymin + additional.labels$ymax) / 2
  }
  
  additional.line.df <- NULL
  if(!is.null(additional.line)){
    additional.line.df <- data.frame(Pos = levels(plot.df$Pos),
                                     Mean = additional.line)
  }
  
  # Methylation boxplot
  p.3 <- ggplot() + theme_bw() +
    theme(text = element_text(size=14),
          axis.text = element_text(face="bold"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "top",
          legend.title = element_blank()) 
  if(!is.null(additional.line.df)){
    p.3 <- p.3 + geom_line(data = additional.line.df, aes(x=Pos, y=Mean, group=1), color = "black", linetype = "dashed")
  }
  p.3 <- p.3 +
    geom_boxplot(data = plot.df, aes(x=Probe, y=Beta, color=Group), alpha = 0.5, width = 0.6)+scale_color_manual(values = allcolors[1:length(Data.group.unique)])+
    labs(x="CpG",y="Beta")+
    coord_cartesian(ylim = c(ylim.min, ylim.max))
  
  leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p.3)), "guide-box") 
  if(!is.null(Add.region.factor.df)){
    p2=geom_rect(data = Add.region.factor.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Region.type), 
                 color = "white")
    
    p.3 <- p.3 + geom_rect(data = Add.region.factor.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Region.type), 
                           color = "white")
  }
  if(!is.null(additional.labels)){
    p.3 <- p.3 + geom_text(aes(x = xmid, y = ymid, label = label), data = additional.labels, size = 4)
  }
  p.3 <- p.3 +
    #geom_text(data = plot.df, aes(x=Pos, label = Strand), y=1, size = 10) +
    scale_shape(solid = FALSE) +
    scale_y_continuous(breaks=seq(0,1,0.2),labels=c(seq(0,1,0.2)))+
    #ylim(ylim.min, 1) +
    #geom_hline(yintercept = 0) +
    #geom_vline(xintercept = 0.4)+
    #geom_segment(aes(x = 0.4, y = 0, xend = 0.4, yend = 1))+
    #geom_hline(yintercept = -0.1, linetype = "dashed") +
    ggtitle(paste(plot.title,", ", "chr",probe.annotation$chr[1],", ", length(unique(plot.df$Probe)), " probes", sep='')) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(), 
          legend.position = "top")
  
  
  x1=factor(plot.df$Probe[1],levels = rownames(Data.i))
  n=1+as.integer(nrow(Data.i)/2)
  x2=factor(plot.df$Probe[n],levels = rownames(Data.i))
  plotNew <- p.3 + 
    annotation_custom(grob = leg1,xmin = x1, xmax = x2, ymin = 0.85*ylim.max, ymax = 0.95*ylim.max)
  p.3=plotNew+p2+theme(plot.title = element_text(hjust = 0.5),
                       legend.title = element_blank(), legend.text=element_text(size=10),
                       legend.direction="horizontal",
                       legend.position = "bottom")+guides(color = "none")+guides(fill = guide_legend(nrow = 1))
  
  if (nrow(Data.i)>50)
  {
    p.3=p.3+theme(axis.text.x=element_blank())
  }
  
  p.output.list$Line <- p
  p.output.list$Box <- p.3
  
  return(p.output.list)
}
