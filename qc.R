# QC
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/qc.R")

plotNumCpgs <- function(props, names.arg) {
  x <- c(props[1:2],NA,props[3:4],NA,props[5:6],NA,props[7])
  names.arg <- c(names.arg[1:2],NA,names.arg[3:4],NA,names.arg[5:6],NA,names.arg[7])
  platform <- c("ERRBS","SSMethylSeq","CpGiant", "WGBS")
  cols <- c("lightblue","lightblue","white","salmon","salmon","white","plum2","plum2","white","lightgreen")
  yaxis <- pretty(c(0,max(pretty(range(x,na.rm=T)))),n=10)
  b <- barplot(x,
    ylim=range(yaxis),
    main="Number of C's in CpG context covered",
    ylab=expression(paste("number of C's (",10^6,")")),
    xlab="platform",
    sub="coverage >10x",
    names.arg=names.arg,
    col=cols,
    space=0,
    cex.names=0.9,
    axes=F)
  axis(2, at=yaxis, labels=yaxis/1000000,las=1)
  at <- seq(1,length(b),by=3)
  at[length(at)] <- b[length(b)]
  mtext(platform, side=1, at=at, line=2) # platform labels
  # text(b,x,x,srt=45,pos=1,offset=2)
  text(b,x,x,pos=1, cex=0.7) # cpg counts in bars
}


plotMeanCpgCov <- function(props, names.arg) {
  x <- c(props[1:2],NA,props[3:4],NA,props[5:6],NA,props[7])
  names.arg <- c(names.arg[1:2],NA,names.arg[3:4],NA,names.arg[5:6],NA,names.arg[7])
  platform <- c("ERRBS","SSMethylSeq","CpGiant", "WGBS")
  cols <- c("lightblue","lightblue","white","salmon","salmon","white","plum2","plum2","white","lightgreen")
  yaxis <- pretty(c(0,max(pretty(range(x,na.rm=T)))),n=10)
  b <- barplot(x,
    ylim=range(yaxis),
    main="Mean coverage per C in CpG context",
    ylab="mean coverage (count)",
    xlab="platform",
    sub="coverage >10x",
    names.arg=names.arg,
    col=cols,
    space=0,
    cex.names=0.9,
    axes=F)
  axis(2, at=yaxis, labels=yaxis, las=1)
  at <- seq(1,length(b),by=3)
  at[length(at)] <- b[length(b)]
  mtext(platform, side=1, at=at, line=2)
  # text(b,x,x,srt=45,pos=1,offset=2)
  text(b,x,x,pos=1)
}


plotMappingEfficiency <- function(props, names.arg) {
  x <- c(props[1:2],NA,props[3:4],NA,props[5:6],NA,props[7])
  names.arg <- c(names.arg[1:2],NA,names.arg[3:4],NA,names.arg[5:6],NA,names.arg[7])
  platform <- c("ERRBS","SSMethylSeq","CpGiant", "WGBS")
  cols <- c("lightblue","lightblue","white","salmon","salmon","white","plum2","plum2","white","lightgreen")
  yaxis <- pretty(c(0,max(pretty(range(x,na.rm=T)))),n=10)
  b <- barplot(x,
    ylim=range(yaxis),
    main="Mapping efficiency",
    ylab="percent (%)",
    xlab="platform",
    names.arg=names.arg,
    col=cols,
    space=0,
    cex.names=0.9,
    axes=F)
  axis(2, at=yaxis, labels=yaxis, las=1)
  at <- seq(1,length(b),by=3)
  at[length(at)] <- b[length(b)]
  mtext(platform, side=1, at=at, line=2) # platform labels
  # text(b,x,x,srt=45,pos=1,offset=2)
  text(b,x,x,pos=1)
}


plotAlignStats <- function(props, names.arg) {
  tmp <- rbind(t(props[,8:10]),props[,7]-apply(props[,8:10],1,sum))
  tmp <- cbind(tmp[,1:2], NA, tmp[,3:4], NA, tmp[,5:6], NA, tmp[,7])
  names.arg <- c(names.arg[1:2],NA,names.arg[3:4],NA,names.arg[5:6],NA,names.arg[7])
  platform <- c("ERRBS","SSMethylSeq","CpGiant", "WGBS")
  rownames(tmp)[4] <- "Rejected"
  legendlabels <- c("Unique alignments","No alignments","Ambiguous mapping","Sequences rejected")
  bcols <- c("green","pink","gray90","blue")
  yaxis <- pretty(c(0,max(props[,7])),n=10)
  b <- barplot(tmp,
    space=0,
    col=bcols,
    ylim=range(yaxis),
    main="Read alignment counts",
    xlab="platform",
    ylab=expression(paste("counts (",10^6,")")),
    names.arg=names.arg,
    cex.names=0.9,
    axes=F)
  axis(2, at=yaxis, labels=yaxis/1000000,las=1)
  at <- seq(1,length(b),by=3)
  at[length(at)] <- b[length(b)]
  mtext(platform, side=1, at=at, line=2) # platform labels
  # legend(b[3]-0.28,yaxis[10],legend=legendlabels,fill=bcols,bg=rgb(1,1,1,alpha=0.9),box.col=rgb(0.7,0.7,0.7,alpha=0.4),cex=0.8)
  legend("topleft",legend=legendlabels,fill=bcols,bg=rgb(1,1,1,alpha=0.9),box.col=rgb(0.7,0.7,0.7,alpha=0.4),cex=1)
}


plotAlignStatsPerc <- function(props, names.arg) {
  tmp <- rbind(t(props[,8:10]),props[,7]-apply(props[,8:10],1,sum))
  tmp <- t(t(tmp) / props[,7] * 100)
  tmp <- cbind(tmp[,1:2], NA, tmp[,3:4], NA, tmp[,5:6], NA, tmp[,7])
  names.arg <- c(names.arg[1:2],NA,names.arg[3:4],NA,names.arg[5:6],NA,names.arg[7])
  platform <- c("ERRBS","SSMethylSeq","CpGiant", "WGBS")
  rownames(tmp)[4] <- "Rejected"
  legendlabels <- c("Unique alignments","No alignments","Ambiguous mapping","Sequences rejected")
  bcols <- c("green","pink","gray90","blue")
  b <- barplot(tmp,
    space=0,
    col=bcols,
    #ylim=range(yaxis),
    main="Read alignment proportions",
    xlab="platform",
    ylab="percent (%)",
    names.arg=names.arg,
    cex.names=0.9)
  at <- seq(1,length(b),by=3)
  at[length(at)] <- b[length(b)]
  mtext(platform, side=1, at=at, line=2) # platform labels  
  #legend(b[5]-0.28,20,legend=legendlabels,fill=bcols,bg=rgb(1,1,1,alpha=0.9),box.col=rgb(0.7,0.7,0.7,alpha=0.4),cex=0.8)
  legend(b[5]-b[1],20,legend=legendlabels,fill=bcols,bg=rgb(1,1,1,alpha=0.9),box.col=rgb(0.7,0.7,0.7,alpha=0.4),cex=1)
}


plotNumCpgvsOtherC <- function(props) {
  cols <- c("blue","blue","red","red","green","green")
  xaxis <- pretty(range(props[,1]))
  yaxis <- pretty(range(props[,4]))
  xlab <- colnames(props)[1]
  par(mar=c(5,5,4,2))
  plot(props[,1],props[,4],
    xlim=range(xaxis),
    ylim=range(yaxis),
    main="CpG's vs other C considered (CHH,CHG)",
    xlab=expression(paste("number of CpG's covered (", 10^6, ")")),
    ylab=expression(paste("total other C considered (", 10^6, ")")),
    sub="coverage >10x",
    col=cols,
    pch=21,
    cex=1.5,
    bg="gray",
    axes=F)
  points(props[,1],props[,4],
    cex=0.8,
    col=cols)
  text(props[,1],props[,4],rownames(props),pos=3)
  axis(1, at=xaxis, labels=xaxis/1000000)
  axis(2, at=yaxis, labels=yaxis/1000000, las=2)
  par(mar=c(5.1,4.1,4.1,2.1))
}

# end QC