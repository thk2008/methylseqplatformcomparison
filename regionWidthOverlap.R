# methylome sequencing region analysis
#
# Figure 1
#
#
if(F){
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/regionWidthOverlap.R")
}
require(GenomicRanges)
require(vioplot)
require("BSgenome.Hsapiens.UCSC.hg19")
# library(beeswarm)
# library(VennDiagram)


getDesignRegionsTotalCpGs <- function(regiongr) {
  require("BSgenome.Hsapiens.UCSC.hg19")
  gr         <- trim(regiongr) # sometimes chrM coordinates cause trouble
  # sometimes there is a C on the end of a region that may be a CG
  # which we see in the methylcalls but not in the designed regions
  tmp <- gr
  start(tmp) <- end(tmp)
  x <- getSeq(Hsapiens, tmp)
  y <- vcountPattern("C", x)
  end(tmp[which(y == 1)]) <- end(tmp[which(y == 1)]) + 1
  end(gr) <- end(tmp)
  regionSeqs <- getSeq(Hsapiens, gr)
  regiongr$numCGss <- vcountPattern("CG", regionSeqs)
  regiongr$numCGds <- regiongr$numCGss * 2
  return(regiongr)
}


plotRegionOverlap <- function(x, ymax) {
  totals <- apply(x, 2, sum)
  yaxis  <- pretty(c(0, ymax))
  b <- barplot(x,
              main="",
              xlab="platform",
              ylab=expression(paste("counts (",10^5,")")),
              ylim=range(yaxis),
              col=c("salmon","lightblue"),
              cex.names=1,
              axes=F)
  axis(2, at=yaxis,labels=yaxis/100000, las=2, cex.axis=1)
  text(b, x[1,],paste(signif(c(x["overlap",1] / totals[1] * 100, x["overlap",2] / totals[2] * 100),digits=3), "%"), pos=1,cex=1.0)
  title(main="target regions overlap",cex.main=1.0)
}

plotRegionOverlapPDF <- function(x, ymax) {
  omar <- par()$mar
  par(mar=c(5.1,6.1,4.1,2.1))
  totals <- apply(x, 2, sum)
  yaxis  <- pretty(c(0, ymax))
  b <- barplot(x,
              main="",
              xlab="",
              ylab="",
              ylim=range(yaxis),
              col=c("salmon","lightblue"),
              cex.names=2.5,
              axes=F)
  axis(2, at=yaxis,labels=yaxis/100000, las=2, cex.axis=2)
  text(b, x[1,],paste(signif(c(x["overlap",1] / totals[1] * 100, x["overlap",2] / totals[2] * 100),digits=3), "%"), pos=1,cex=2.0)
  title(main="Target regions overlap", cex.main=2.0,
        xlab="platform",        cex.lab=2,
        ylab=expression(paste("counts (",10^5,")")))
  par(mar=omar)
}


plotRegionWidthBoxplot <- function(region.widths) {
  b <- boxplot(region.widths,
    main="Target region length",
    xlab="platform",
    ylab="length (bases)",
    outline=F,
    col="gray90",
    ylim=c(50,550),
    axes=F)
  box()
  axis(2,at=seq(50,550,by=50),las=2)
  axis(1,at=1:3,labels=names(region.widths))
  mtext("outliers removed",adj=1)
}

plotRegionWidthBoxplotPDF <- function(region.widths) {
  omar <- par()$mar
  par(mar=c(5.1,6.1,4.1,2.1))
  b <- boxplot(region.widths,
    main="",
    xlab="",
    ylab="",
    outline=F,
    col="gray90",
    ylim=c(50,550),
    axes=F)
  box()
  axis(2,at=seq(50,550,by=50),las=2,cex.axis=2)
  axis(1,at=1:3,labels=names(region.widths),cex.axis=2)
  title(main="Target region length", cex.main=2,
        xlab="platform",      cex.lab=2)
  title(ylab="length (bases)", line=4, cex.lab=2)
  legend("topright","outliers removed",cex=1.5,bty="n")
  par(mar=omar)
}




plotRegionWidthVioplot <- function(region.widths) {
  violengths <- list()
  for (i in 1:length(region.widths)) {
    y <- region.widths[[i]]
    violengths[[i]] <- y[y <= 1000]
  }
  vioplot(violengths[[1]],violengths[[2]],violengths[[3]],
    h=25,
    rectCol="gray50",
    names=names(region.widths),
    col="gray90")
  title(main="Target region length",xlab="plaform",ylab="length (bases)")
  mtext("lengths > 1000 removed",adj=1)
}



getOverlapMatList <- function(regionsList) {
    mat <- cbind(c(1,1,2), c(2,3,3))
    overlapList <- apply(mat, 1, function(i) {
            labelA      <- names(regionsList)[ i[1] ]
            labelB      <- names(regionsList)[ i[2] ]
            x           <- regionsList[[ i[1] ]]
            y           <- regionsList[[ i[2] ]]
            o           <- findOverlaps(x, y, minoverlap=10)
            totalA      <- length(x)
            totalB      <- length(y)
            numoverlapA <- length(unique(queryHits(o)))
            numoverlapB <- length(unique(subjectHits(o)))

            overlaps <- c(totalA, totalB, numoverlapA, numoverlapB,length(o))
            names(overlaps) <- c("totalA", "totalB","numberOverlapA","numberOverlapB","numberOverlap")
            
            list(overlaps=overlaps, A=labelA, B=labelB)

          })
          
  return(overlapList)
}


plotOverlapTape <- function(overlapMatList) {
  tmp  <- sapply(overlapMatList, "[[", "overlaps")
  ymax <- max(pretty(range(tmp[1:2,]),n=10))
  omar <- par()$mar
  par(mar=c(5.1,6.6,4.1,6.6))
  b <- barplot(c(0,0,0),
    horiz=T,
    xlim=c(-ymax,ymax),
    xlab=expression(paste("number of regions (",10^3,")")),
    names.arg=sapply(overlapMatList, "[[", "B"),
    ylim=c(0,4),
    las=2,
    axes=F)
  title(main="Target region overlap",cex.main=2.0)
  axis(4, at=b, labels=sapply(overlapMatList, "[[", "A"), las=2, tick=F)
  yaxis <- pretty(c(-ymax,ymax),n=21)
  axis(1, at=yaxis,labels=abs(yaxis / 1000))
  segments(yaxis,0,yaxis,3.8,lty=2,col="gray90")
  barplot(tmp[1,],  names.arg="",horiz=T,add=T,axes=F,col="gray90") # totalA
  barplot(-tmp[2,], names.arg="",horiz=T,add=T,axes=F,col="gray90") # totalB
  barplot(tmp[3,],  names.arg="",horiz=T,add=T,axes=F,col="green")  # overlapA
  barplot(-tmp[4,], names.arg="",horiz=T,add=T,axes=F,col="green")  # overlapB
  segments(0,0,0,3.8,col="gray80")
  segments(0,0,0,3.8,lty=2)
  axis(2, at=b-0.1, labels=paste(round(tmp[4,] / tmp[2,] * 100, digits=1),"%",sep=""),las=2,tick=F,pos=-ymax)
  axis(4, at=b-0.1, labels=paste(round(tmp[3,] / tmp[1,] * 100, digits=1),"%",sep=""),las=2,tick=F,pos=ymax) 
  par(mar=omar)
}



##########
#        #
#  MAIN  #
#        #
##########
datadir <- "data"
makepdf <- F

atime  <- format(Sys.time(), "%Y%m%d")
outdir <- paste("regionWidthOverlap",atime, sep="-")
if (!file.exists(outdir)) {
  dir.create(outdir)
}
cat("output directory is [", outdir, "]\n")

if (!exists("regionsList")) {
  rda <- "regionsList.rda"
  if (!file.exists(rda)){
    cat("creating regionsList.rda\n")
    source("/home/thk2008/bin/methylseqplatformcomparison/scripts/getDataUtils.R")
    regionFiles <- paste(datadir, c("MethylSeq.84Mb.bed",
                                    "130912_HG19_CpGiant_4M_EPI_noRandom.bed", # was 130912_HG19_CpGiant_4M_EPI.bed
                                    "MspI-fragments-hg19-20140919.bed"),
                   sep="/")
    names(regionFiles) <- c("SSMethylSeq","CpGiant","MspI_84-334")
    regionsList <- getRegionsList(regionFiles)
    regionsList[[1]] <- getDesignRegionsTotalCpGs(regionsList[[1]])
    regionsList[[2]] <- getDesignRegionsTotalCpGs(regionsList[[2]])
    regionsList[[3]] <- getDesignRegionsTotalCpGs(regionsList[[3]])
    save(regionsList, file=rda)
  } else {
    cat("loading",rda,"\n")
    load(rda)
  }
}

region.widths <- lapply(regionsList, width)

if (!exists("overlapMatList")) {
  rda <- "overlapMatList.rda"
  if (!file.exists(rda)) {
    cat("creating overlapMatList.rda\n")
    overlapMatList <- getOverlapMatList(regionsList)
    save(overlapMatList, file="overlapMatList.rda")
  } else {
    cat("loading",rda,"\n")
    load(rda)
  }
}

run <- function() {

  if (makepdf) { 
    pdf(paste(outdir, "regionWidthBoxplot.pdf", sep="/"),pointsize=11)
     plotRegionWidthBoxplotPDF(region.widths)
    dev.off()
  } else {
    dev.new()
    plotRegionWidthBoxplot(region.widths)
  }

  # if (makepdf) { pdf(paste(outdir, "regionWidthVioplot.pdf", sep="/"),pointsize=11) } else { dev.new() }
    # plotRegionWidthVioplot(region.widths)
  # if (makepdf) { dev.off() }

  ymax <- max(sapply(regionsList, length))

  for (i in 1:length(overlapMatList)) {
    if (makepdf) {
      pdf(paste(outdir, paste("regionOverlap_",paste(colnames(overlapMatList[[i]][[1]]),collapse="-"),".pdf",sep=""), sep="/"),pointsize=11)
      plotRegionOverlapPDF(overlapMatList[[i]][[1]], ymax)
      dev.off()
      pdf(paste(outdir, paste("regionOverlap_",paste(colnames(overlapMatList[[i]][[1]]),collapse="-"),".pdf",sep=""), sep="/"),pointsize=11)
      plotOverlapTape(overlapMatList)
      dev.off()
    } else {
      dev.new()
      plotRegionOverlap(overlapMatList[[i]][[1]], ymax)
    }
  }
  
  
  if (makepdf) {
    pdf(paste(outdir, paste("regionOverlapTape.pdf",sep=""), sep="/"),pointsize=11)
      plotOverlapTape(overlapMatList)
    dev.off()
  } else {
    dev.new()
    plotOverlapTape(overlapMatList)
  }

}




if(F) {
     # draw.pairwise.venn(area1, area2, cross.area, category = rep("", 2),
         # euler.d = TRUE, scaled = TRUE, inverted = FALSE,
         # ext.text = TRUE, ext.percent = rep(0.05, 3), lwd =
         # rep(2, 2), lty = rep("solid", 2), col = rep("black",
         # 2), fill = NULL, alpha = rep(0.5, 2), label.col =
         # rep("black", 3), cex = rep(1, 3), fontface =
         # rep("plain", 3), fontfamily = rep("serif", 3), 
         # cat.pos = c(-50, 50), cat.dist = rep(0.025, 2), cat.cex =
         # rep(1, 2), cat.col = rep("black", 2), cat.fontface =
         # rep("plain", 2), cat.fontfamily = rep("serif", 2),
         # cat.just = rep(list(c(0.5, 0.5)), 2), cat.default.pos = "outer",
         # cat.prompts = FALSE, ext.pos = rep(0, 2),
         # ext.dist = rep(0, 2), ext.line.lty = "solid",
         # ext.length = rep(0.95, 2), ext.line.lwd = 1,
         # rotation.degree = 0, rotation.centre = c(0.5, 0.5),
         # ind = TRUE, sep.dist = 0.05, offset = 0, cex.prop =
         # NULL, print.mode = "raw", sigdigs = 3, ...)

	library(VennDiagram)
	load("overlapMatList.rda")

	makepdf <- F

  atime  <- format(Sys.time(), "%Y%m%d")
  outdir <- paste("regionWidthOverlap",atime, sep="-")
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  cat("output directory is [", outdir, "]\n")

	
# x <- overlapMatList[[1]]
# dev.off()
# vn <- draw.pairwise.venn( x$overlaps[1], x$overlaps[2], mean(x$overlaps[3:4]),
	# label.col=rgb(1,1,1,alpha=0),
	# cat.dist=rep(0.05, 2),
	# cat.cex=rep(1.5,2),
	# fill=c("lightblue","salmon"),
	# margin=c(0.04,0.04,0.04,0.04),
	# category=c(paste(x$A,"\n",round(x$overlaps[3]/x$overlaps[1]*100,digits=2),"%",sep=""),
	           # paste(x$B,"\n",round(x$overlaps[4]/x$overlaps[2]*100,digits=2),"%",sep="")))
# dev.off()
# plot(0, xlim=c(-10,10), ylim=c(-10,10), type="n", xlab="",ylab="", axes=F)
# grid.draw(vn)
# title(main="Region overlap",cex.main=2)

	
	x <- overlapMatList[[1]]
dev.off()
vn <- draw.pairwise.venn( x$overlaps[1], x$overlaps[2], mean(x$overlaps[3:4]),
	label.col=rgb(1,1,1,alpha=0),
	cat.dist=rep(0.05, 2),
	cat.cex=rep(1.5,2),
	fill=c("lightblue","salmon"),
	margin=c(0.04,0.04,0.04,0.04),
	category=c(paste(x$A,"\n",round(x$overlaps[3]/x$overlaps[1]*100,digits=2),"%",sep=""),
	           paste(x$B,"\n",round(x$overlaps[4]/x$overlaps[2]*100,digits=2),"%",sep="")))
dev.off()
plot(0, xlim=c(-10,10), ylim=c(-10,10), type="n", xlab="",ylab="", axes=F)
grid.draw(vn)
title(main="Region overlap",cex.main=2)


for (i in 
x <- overlapMatList[[1]]
dev.off()
vn <- draw.pairwise.venn( x$overlaps[1], x$overlaps[2], mean(x$overlaps[3:4]),
	label.col=rgb(1,1,1,alpha=0),
	cat.dist=rep(0.05, 2),
	cat.cex=rep(1.5,2),
	fill=c("lightblue","salmon"),
	margin=c(0.04,0.04,0.04,0.04),
	category=c(paste(x$A,"\n",round(x$overlaps[3]/x$overlaps[1]*100,digits=2),"%",sep=""),
	           paste(x$B,"\n",round(x$overlaps[4]/x$overlaps[2]*100,digits=2),"%",sep="")))
dev.off()
plot(0, xlim=c(-10,10), ylim=c(-10,10), type="n", xlab="",ylab="", axes=F)
grid.draw(vn)
title(main="Region overlap",cex.main=2)







	
}

  for (i in 1:length(overlapMatList)) {
    x <- overlapMatList[[i]]
    #dev.off()
    cpos <- c(0,0)
    if (i == 1) { cpos <- c(210,150) }
    vn <- draw.pairwise.venn( x$overlaps[1], x$overlaps[2], mean(x$overlaps[3:4]),
      label.col=rgb(1,1,1,alpha=0),
      cat.dist=rep(0.05, 2),cat.pos=cpos,
      cat.cex=rep(1.5,2),
      fill=c("lightblue","salmon"),
      margin=c(0.04,0.04,0.04,0.04),
      category=c(paste(x$A,"\n",round(x$overlaps[3]/x$overlaps[1]*100,digits=2),"%",sep=""),
                 paste(x$B,"\n",round(x$overlaps[4]/x$overlaps[2]*100,digits=2),"%",sep="")))
    #dev.off()
    
    if (makepdf) { pdf(paste(outdir, paste("regionOverlapVenn_", x$A, "-", x$B,".pdf",sep=""), sep="/"), pointsize=11) }
    
    plot(0, xlim=c(-10,10), ylim=c(-10,10), type="n", xlab="",ylab="", axes=F)
    grid.draw(vn)
    title(main="Region overlap",cex.main=2)
    
    if (makepdf) { dev.off() }
  }
}

