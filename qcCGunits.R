# QC
if(FALSE){
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/qcCGunits.R")
}

library(GenomicRanges)
library(beeswarm)

plotNumCGunits <- function(props, names.arg, main="") {
  x <- c(props[1:2],NA,props[3:4],NA,props[5:6],NA,props[7])
  names.arg <- c(names.arg[1:2],NA,names.arg[3:4],NA,names.arg[5:6],NA,names.arg[7])
  platform <- c("ERRBS","SSMethylSeq","CpGiant", "WGBS")
  cols <- c("lightblue","lightblue","white","salmon","salmon","white","plum2","plum2","white","lightgreen")
  yaxis <- pretty(c(0,max(pretty(range(x,na.rm=T)))),n=10)
  b <- barplot(x,
    ylim=range(yaxis),
    main=main,
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


plotMeanCGunitCoverage <- function(props, names.arg, main="") {
  x <- c(props[1:2],NA,props[3:4],NA,props[5:6],NA,props[7])
  names.arg <- c(names.arg[1:2],NA,names.arg[3:4],NA,names.arg[5:6],NA,names.arg[7])
  platform <- c("ERRBS","SSMethylSeq","CpGiant", "WGBS")
  cols <- c("lightblue","lightblue","white","salmon","salmon","white","plum2","plum2","white","lightgreen")
  yaxis <- pretty(c(0,max(pretty(range(x,na.rm=T)))),n=10)
  b <- barplot(x,
    ylim=range(yaxis),
    main=main,
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

# plotCGcovDistBox <- function(mcgrlist,main="") {
  # cgcov <- list()
  # for (i in 1:length(mcgrlist)) {
    # cgcov[[i]] <- mcgrlist[[i]]$coverage
  # }
  # meancov <- sapply(cgcov, mean)
  # names(cgcov) <- names(mcgrlist)
  # bcols <- c("blue1","blue3","salmon1","salmon3","plum","plum4","lightgreen")
  # omar <- par()$mar
  # par(mar=c(9.1, 4.1, 4.1, 2.1))
  # b <- boxplot(cgcov,
   # main=main,
   # xlab="",
   # ylab="counts",
   # outline=F,
   # col=bcols,
   # las=2)
  # mtext("outliers removed",adj=1)
  # points(seq(1:length(meancov)), meancov, pch=3)
  # legend("topright",legend="mean",bty="n",pch=3)
  # par(mar=omar)
# }

plotCGcovDistBox <- function(mcgrlist,main="") {
  cgcov <- list()
  for (i in 1:length(mcgrlist)) {
    cgcov[[i]] <- mcgrlist[[i]]$coverage
  }
  meancov <- sapply(cgcov, mean)
  names(cgcov) <- names(mcgrlist)
  platform <- c("ERRBS","SSMethylSeq","CpGiant", "WGBS")
  bcols <- c("blue1","blue3","salmon1","salmon3","plum","plum4","lightgreen")
  b <- boxplot(cgcov,
    main=main,
    xlab="",
    ylab="counts",
    outline=F,
    col=bcols,
    las=2,
    axes=F)
  box()
  axis(2)
  axis(1, at=seq(1,length(cgcov)), labels=abtext,tick=F,line=-0.5)
  at=c(seq(1.5,length(cgcov),by=2),length(cgcov))
  mtext(platform, side=1, at=c(1.5,3.5,5.5,7), line=2) # platform labels
  mtext("outliers removed",adj=1)
  points(seq(1:length(meancov)), meancov, pch=3)
  legend("topright",legend="mean",bty="n",pch=3)
}






plotRegionCoverageSwarm <- function(percRegCovered, percRegCGcovered, main="") {
  cols <- c("blue1","blue3","salmon1","salmon3","plum","plum4")
  plot(0, 
    type="n", 
    main=main, 
    xlab="", 
    ylab="percent (%)",
    xlim=c(0.5,2.5),
    ylim=c(60,100),
    axes=F)
  y <- seq(60,100,by=5)
  segments(0, y, 3, y, col=rgb(0.5,0.5,0.5,alpha=0.5),lwd=0.5)
  axis(1, at=c(1,2),
    labels=c("Regions covered","CpGs covered in regions"),
    tick=F)
  axis(2, at=y, las=1)
  
  beeswarm(list(percRegCovered,percRegCGcovered),
    pwcol=rep(cols,2), 
    pch=21, 
    pwbg=rep(cols,2),
    cex=1.5,
    las=1,
    add=T)
  legend("bottomleft",
    names(percRegCovered),
    bg="snow",
    fill=cols,
    border=cols,
    box.col=rgb(0.6,0.6,0.6,alpha=0.5))
  box()
}


plotRegionCoverage <- function(percRegCovered, percRegCGcovered, main="") {
  cols <- c("blue1","blue4","salmon1","salmon3","plum","plum4")

  plot(percRegCovered, percRegCGcovered,
    main=main, 
    xlab="percent Regions covered (%)", 
    ylab="percent Region CG's covered (%)",
    pch=19,
    col=cols,
    las=2)
  x  <- range(pretty(range(percRegCovered)))
  xx <- seq(x[1],x[2],by=5)
  y  <- range(pretty(range(percRegCGcovered)))
  yy <- seq(y[1],y[2],by=5)
  segments(x[1], yy, x[2], yy, col=rgb(0.5,0.5,0.5,alpha=0.5),lwd=0.5)
  segments(xx, y[1], xx, y[2], col=rgb(0.5,0.5,0.5,alpha=0.5),lwd=0.5)
  
  legend("bottomright",
    names(percRegCovered),
    bg="snow",
    fill=cols,
    border=cols,
    box.col=rgb(0.6,0.6,0.6,alpha=0.5))
  box()
}



getRegionCoverage <- function(regionsList, mcgrlist) {
  mat <- cbind(c(3,3,1,1,2,2), c(1,2,3,4,5,6))  
  perccovered <- apply(mat, 1, function(i) {
    x <- findOverlaps(regionsList[[ i[1] ]], mcgrlist[[ i[2] ]])
    y <- length(unique(queryHits(x)))
    perccovered <- y / length(regionsList[[ i[1] ]]) * 100
  })
  names(perccovered) <- names(mcgrlist)[ mat[,2] ]  
  return(perccovered)
}

getRegionCGCoverage <- function(regionsList, mcgrlist) {
  mat <- cbind(c(3,3,1,1,2,2), c(1,2,3,4,5,6))  
  perccovered <- apply(mat, 1, function(i) {
    x <- findOverlaps(regionsList[[ i[1] ]], mcgrlist[[ i[2] ]])
    y <- length(unique(subjectHits(x)))
    perccovered <- y / sum(regionsList[[ i[1] ]]$numCGss) * 100
  })
  names(perccovered) <- names(mcgrlist)[ mat[,2] ]  
  return(perccovered)
}


getRegionCoverageStats <- function(regionsList, mcgrlist) {
  mat <- cbind(c(3,3,1,1,2,2), c(1,2,3,4,5,6))  
  statsmat <- numeric()
  statsmat <- t(apply(mat, 1, function(i) {
      x <- findOverlaps(regionsList[[ i[1] ]], mcgrlist[[ i[2] ]])
      statsmat <- rbind(statsmat,
       c(length(regionsList[[ i[1] ]]), sum(regionsList[[ i[1] ]]$numCGss), length(mcgrlist[[ i[2] ]]), length(unique(queryHits(x))), length(unique(subjectHits(x)))))
    }))
    
  colnames(statsmat) <- c("totalRegions", "totalCGregions", "totalCGunits","overlapRegions","overlapCGunits")
  rownames(statsmat) <- names(mcgrlist)[mat[,2]]
  return(statsmat)
}

#
# MAIN
#
  
if (!exists("regionsList")) {
  cat("loading regionsList from regionsList.rda\n")
  load("regionsList.rda")
}
if (!exists("mcgrlist")) {
  cat("loading mcgrlist from mcgrlist_CGunits10x.rda\n")
  load("mcgrlist_CGunits10x.rda")
}


makepdf <- F

atime  <- format(Sys.time(), "%Y%m%d")
outdir <- paste("qcCGunits",atime, sep="-")
if (!file.exists(outdir)) {
  dir.create(outdir)
}
cat("output directory is [", outdir, "]\n")


abtext <- c("A","B","A_opt","B_min","A_opt","B_min","A")

run <- function() {
  
  
  if (makepdf) { pdf(paste(outdir, "numberCGunitsCovered.pdf", sep="/"),pointsize=11) }
    props <- sapply(mcgrlist, length)
    plotNumCGunits(props, abtext, main="Number of CG units covered")
  if (makepdf) { dev.off() } else { dev.new() }
  
  if (makepdf) { pdf(paste(outdir, "meanCGunitsCoverage.pdf", sep="/"),pointsize=11) }
    props <- sapply(mcgrlist, function(i) { x <- round(mean(i$coverage),digits=2) })
    plotMeanCGunitCoverage(props, abtext, main="Mean CG unit coverage")
  if (makepdf) { dev.off() } else { dev.new() }
  
  if (makepdf) { pdf(paste(outdir, "medianCGunitsCoverage.pdf", sep="/"),pointsize=11) }
    props <- sapply(mcgrlist, function(i) { x <- round(median(i$coverage),digits=2) })
    plotMeanCGunitCoverage(props, abtext, main="Median CG unit coverage")
  if (makepdf) { dev.off() } else { dev.new() }
  
  if (makepdf) { pdf(paste(outdir, "cgUnitCoverageBox.pdf", sep="/"),pointsize=11) }
   #cat("plotting box plot of cg units \n")
   plotCGcovDistBox(mcgrlist, main="CG unit coverage >=10x")
  if (makepdf) { dev.off() }  

  

  
}


saveMat <- function(props, file="stats") {
  afile <- paste(file, "csv", sep=".")
  cat(",",file=afile)
  ow <- options("warn")
  options(warn = -1)
  write.table(props,file=afile,append=T,quote=F,sep=",")
  options(ow) # reset
}


plotNumberCGunitCoverage <- function(statsmatRegionCoverage, percRegCGcovered) {
  
  abtext <- c("A","B","A_opt","B_min","A_opt","B_min","A") 
  
  x <- c(as.vector(rbind(statsmatRegionCoverage[,2],0)),length(mcgrlist[["WGBS"]]))
  x <- c(x[1:4],NA,x[5:8],NA,x[9:12],NA,x[13])
  names.arg <- c(abtext[1:2],NA,abtext[3:4],NA,abtext[5:6],NA,abtext[7])
  platform <- c("ERRBS","SSMethylSeq","CpGiant", "WGBS")
  yaxis <- pretty(c(0,max(pretty(range(x,na.rm=T)))),n=10)
  bcols <- c("gray50","lightblue","salmon","orange1","red")
  
  omar <- par()$mar
  par(mar=c(5.1,4.5,4.1,5.1))
  b <- barplot(x,
        main="Number of CG units",
        ylab=expression(paste("number of CG units (",10^6,")")),
        xlab="platform",  
        ylim=range(yaxis),
        sub="coverage >10x",
        col=bcols[1], border=NA,
        lwd=2,
        axes=F)
  axis(2, at=yaxis, labels=yaxis/1000000,las=1)
  at <- c(rowMeans(cbind(b[1:15],b[2:16]))[c(1,3,6,8,11,13)],b[length(b)])
  axis(1,at=at,labels=abtext,tick=F)
  at <- c(rowMeans(cbind(b[1:15],b[2:16]))[c(2,7,12)],b[length(b)])
  mtext(platform, side=1, at=at, line=2) # platform labels
  y <- c(as.vector(rbind(statsmatRegionCoverage[,5], statsmatRegionCoverage[,3] - statsmatRegionCoverage[,5])),length(mcgrlist[["WGBS"]]))
  y <- c(y[1:4],NA,y[5:8],NA,y[9:12],NA,y[13])
  barplot(y, add=T,col=c(bcols[2:3],bcols[2:3],NA,bcols[2:3],bcols[2:3],NA,bcols[2:3],bcols[2:3],NA,bcols[4]), border=NA, axes=F)
  
  segments(b[1],tail(seq(0,10000000, by=1000000)),b[15],tail(seq(0,10000000, by=1000000)),col=rgb(1,0,0,alpha=0.1))    
  
  axis(4, at=seq(0,10000000, by=1000000),
    labels=seq(0,10000000, by=1000000) / 100000,
    col=bcols[5],
    las=2)
  mtext("percent CG units in regions covered", side=4, line=2.2, adj=0.2)
  points(b[ c(1, 3, 6, 8, 11, 13) ], percRegCGcovered * 100000,
    col=bcols[5],
    pch=15,
    lwd=2)

    
    
  legend("topleft",bty="n",
    legend=c("total CG units in predicted regions","covered CG units in regions","covered CG units outside regions","covered CG units","percent CG units in regions covered"),
    fill=bcols)
  
  par(mar=omar)
  
}



runRegionRecovery <- function() {
  
  if (!exists("statsmatRegionCoverage")) {
    rda <- "statsmatRegionCoverage.rda"
    if (!file.exists(rda)){
      statsmatRegionCoverage <- getRegionCoverageStats(regionsList, mcgrlist)
      save(statsmatRegionCoverage, file="statsmatRegionCoverage.rda")
      saveMat(statsmatRegionCoverage, "statsmatRegionCoverage")
    } else {
      cat("loading",rda,"\n")
      load(rda)
    }      
  }
  
  if (!exists("percRegCovered")) {
    rda <- "percRegCovered.rda"
    if (!file.exists(rda)){
      cat("computing percRegCovered\n")
      percRegCovered <- getRegionCoverage(regionsList, mcgrlist[1:6])
      save(percRegCovered,file=rda)
    } else {
      cat("loading",rda,"\n")
      load(rda)
    }
  }
  if (!exists("percRegCGcovered")) {
    rda <- "percRegCGcovered.rda"
    if (!file.exists(rda)){
      cat("computing percRegCGcovered\n")
      percRegCGcovered <- getRegionCGCoverage(regionsList, mcgrlist[1:6])
      save(percRegCGcovered,file=rda)
    } else {
      cat("loading",rda,"\n")
      load(rda)
    }
  }


  if (makepdf) { pdf(paste(outdir,"numberCGunitCoverage.pdf",sep="/")) } else { dev.new() }
    #plotNumberCGunitCoverage(statsmatRegionCoverage)
    plotNumberCGunitCoverage(statsmatRegionCoverage, percRegCGcovered)
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"targetRegAndCGcoverageSwarm.pdf",sep="/")) } else { dev.new() }
    plotRegionCoverageSwarm(percRegCovered, percRegCGcovered, main="Percent coverage of regions and C's in regions")
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"targetRegAndCGcoverageScatter.pdf",sep="/")) } else { dev.new() }
    plotRegionCoverage(percRegCovered, percRegCGcovered, main="Percent coverage of regions and C's in regions")
  if (makepdf) { dev.off() }



}


runCoverageDistribution <- function() {

  mat <- cbind(c(3,3,1,1,2,2),c(1,2,3,4,5,6))
  pcov <- list()
  pcov <- apply(mat,1,function(i) {
    print(i)
    overlaps <- findOverlaps(regionsList[[ i[1] ]], mcgrlist[[ i[2] ]])  # mspi
    x <- table(queryHits(overlaps))
    y <- cbind(regionsList[[ i[1] ]]$numCGss, 0)
    y[as.numeric(names(x)),2] <- x
    x <- y[,2] / y[,1] * 100
    if(length(which(x > 100)) != 0) {
      stop("something wrong in [",i,"]")
    }
    pcov[[ i[2] ]] <- x[!is.nan(x)]
  })
  names(pcov) <- names(mcgrlist)[mat[,2]]
  save(pcov, file="pcov.rda")  

  cols <- c("blue1","blue3","salmon1","salmon3","plum","plum4")
  
  if (makepdf) { pdf(paste(outdir,"targetRegCGunitCoverageDensity.pdf",sep="/")) } else { dev.new() }

  plot(0,type="n",xlim=c(-5,105), ylim=c(0,0.4),
    main="Target regions CG unit coverage",
    xlab="percent coverage of region (%)",
    ylab="density")
    
  for (i in 1:length(pcov)) {
    lines(density(pcov[[i]],bw=1),col=cols[i])
  }
  
  legend("topleft", names(pcov), lty=1, col=cols, bty="n") 
  
  if (makepdf) { dev.off() } 
}

runCoverageDistributionHeatmap <- function() {
  
  makepdf <- TRUE
  outdir <- "qcCGunits-20160606"
  load("pcov.rda")  

  mat <- numeric()
  for (i in 1:length(pcov)) {
    x <- hist(pcov[[i]], breaks=100, plot=FALSE)
    mat <- cbind(mat, x$counts)
  }
  colnames(mat) <- names(pcov)
  
  omar <- par()$mar
  
  rgb.palette <- colorRampPalette(c("white","lightblue","blue","black"), space = "rgb")
  breaks <- c(0,500,1000,5000,10000,50000,100000,150000,250000,350000)
  inset <- -0.16
  if (makepdf) { 
    pdf(paste(outdir,"targetRegCGunitCoverageHeatmap.pdf",sep="/"))
    inset <- -0.28
  } else { dev.new() }
  
  par(mar=c(5.1,9,4.1,6.1))
  
  image(mat,   
   col=rgb.palette(length(breaks)-1),   
   breaks=breaks,
   main="Target regions CG unit coverage",
   xlab="percent coverage CG units in region (%)",
   axes=F)
   segys <- seq(-0.1,1.1,by=0.2)
   segments(-0.1, segys, 100.1, segys, lwd=4, col="white", lend=2)
   axis(1,at=seq(0,1,by=0.1), labels=seq(0,100,by=10))
   axis(2, at=seq(0,1, by=0.2),labels=colnames(mat), las=2, tick=F)
  legend("right",
    legend=rev(breaks),
    col=rev(rgb.palette(length(breaks))),
    title="Density",
    box.col="gray90",
    pch=22,
    pt.cex=3,
    pt.bg=rev(rgb.palette(length(breaks))),
    inset=inset,
    xpd=T)
    
  if (makepdf) { dev.off() }

  par(mar=omar)
  
 
  
  
  
  
  
  
  
  percmat <- numeric()
  for (i in 1:length(pcov)) {
    x <- hist(pcov[[i]], breaks=100, plot=FALSE)
    percmat <- cbind(percmat, round(x$counts/sum(x$counts) * 100,digits=3))
  }
  colnames(percmat) <- names(pcov)
  
  omar <- par()$mar
  
  rgb.palette <- colorRampPalette(c("white","lightblue","orange","blue","black"), space = "rgb")
  breaks <- c(0,0.5,1,2,5,10,20,50,75,90,100)
  inset <- -0.16
  if (makepdf) {
    pdf(paste(outdir,"targetRegCGunitCoverageHeatmapPerc.pdf",sep="/"))
    inset <- -0.28
  } else { dev.new() }
  
  par(mar=c(5.1,9,4.1,6.1))
  
  image(percmat,   
   col=rgb.palette(length(breaks)-1),   
   breaks=breaks,
   main="Percent of target regions with percent CG unit coverage",
   xlab="percent coverage CG units in region (%)",
   axes=F)
   segys <- seq(-0.1,1.1,by=0.2)
   segments(-0.1, segys, 100.1, segys, lwd=4, col="white", lend=2)
   axis(1,at=seq(0,1,by=0.1), labels=seq(0,100,by=10))
   axis(2, at=seq(0,1, by=0.2),labels=colnames(percmat), las=2, tick=F)
    
  legend("right",
    legend=rev(breaks),
    col=rev(rgb.palette(length(breaks))),
    title=paste(expression("Percent\nof\ndataset\nregions\n(%)")),
    box.col="white",
    pch=22,
    pt.cex=3,
    pt.bg=rev(rgb.palette(length(breaks))),
    inset=inset,
    xpd=T)

  if (makepdf) { dev.off() }

  par(mar=omar)
  
  
  
  
# pdf("targetRegCGunitCoverage_ecdf.pdf")  
# plot(ecdf(pcov[[1]]),col=rgb(0,0,1,alpha=0.5),verticals=T,pch="",lwd=3,
  # main="Target regions CG unit coverage",
  # xlab="percent CpG coverage in region (%)")
# plot(ecdf(pcov[[2]]),add=T,col=rgb(0,0,0.8,alpha=0.5),verticals=T,pch="",lwd=3)
# plot(ecdf(pcov[[3]]),add=T,col=rgb(0,1,0,alpha=0.5),verticals=T,pch="",lwd=2)
# plot(ecdf(pcov[[4]]),add=T,col=rgb(0,0.8,0,alpha=0.5),verticals=T,pch="",lwd=2)
# plot(ecdf(pcov[[5]]),add=T, col=rgb(1,0,0,alpha=0.5),verticals=T,pch="",lwd=1)
# plot(ecdf(pcov[[6]]),add=T, col=rgb(0.8,0,0,alpha=0.5),verticals=T,pch="",lwd=1)
# legend("topleft", legend=c("errbs","ssms","cpgiant"), col=c("blue","green","red"),lwd=2)
# dev.off()
  
  
  
  
}




if(F){
library(VennDiagram)

x <- union(mcgrlist[[1]], mcgrlist[[2]])
y <- union(mcgrlist[[3]], mcgrlist[[4]])
z <- union(mcgrlist[[5]], mcgrlist[[6]])

xx <- paste(seqnames(x), start(x), sep=".")
yy <- paste(seqnames(y), start(y), sep=".")
zz <- paste(seqnames(z), start(z), sep=".")



venn <- c(length(x), length(y), length(z),
                 length(intersect(xx,yy)),
                 length(intersect(yy,zz)),
                 length(intersect(xx,zz)),
                 length(intersect(zz, intersect(xx,yy))))
a <- round(venn/1000000, digits=2)
                 
draw.triple.venn(a[1],a[2],a[3],a[4],a[5],a[6],a[7])




cols <- c("lightblue","pink","plum")

tv <- draw.triple.venn(a[1], a[2], a[3], a[4], a[5], a[6], a[7],fill=c("lightblue","salmon","plum2"),category=c("ERRBS", "SSMethylSeq", "CpGiant"),margin=c(0.04,0.04,0.04,0.04), cat.cex = rep(2, 3),cex = rep(2, 7))

pdf("vennOverlap.pdf",width=10,height=10)

par(mar=c(1,1,4,1))
plot(0, xlim=c(-10,10), ylim=c(-10,10), type="n", xlab="",ylab="", axes=F)
grid.draw(tv)
title(main="Overlap of CpG sites (10^6)",cex.main=2)
par(mar=c(5.1,4.1,4.1,2.1))

dev.off()
}






# end QC

