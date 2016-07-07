

if (FALSE) {

source("/scratch/thk2008/epi3bin/methylseqplatformcomparison/scripts/platformComparisonUtils.R")
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/platformComparisonUtils.R")

library(corrplot)
library(GenomicRanges)
library(methylKit)
library(VennDiagram)
library(Rsamtools)
library("BSgenome.Hsapiens.UCSC.hg19")

}


if (!exists("datadir")) {
  if (Sys.info()["nodename"] == "localhost.localdomain") {
    datadir <- "/scratch/thk2008/methylomecapture/data"
  }
  if (Sys.info()["nodename"] == "pc142472.med.cornell.edu") {
    datadir <- "/scratch001/thk2008/methylomecapture/data"
  }
  if (Sys.info()["nodename"] == "epicore03.pbtech") {
    datadir <- "/scratch001/thk2008_dat/methylomecapture/data"
  }
}

slabels <- c("ERRBS_A","ERRBS_B","SSMethylSeq_A_opt","SSMethylSeq_B_min","CpGiant_A_opt","CpGiant_B_min","WGBS")

summaryfiles <- paste( datadir,
  c("bisseq_pipeline_summary.IMR_90.txt",
    "bisseq_pipeline_summary.IMR90_A_TruSeq.txt",
    "bisseq_pipeline_summary.IMR90_1.txt",
    "bisseq_pipeline_summary.IMR90-_2.txt",
    "bisseq_pipeline_summary.IMR90_Roche_1ug.txt",
    "bisseq_pipeline_summary.IMR90_Roche_0_25ug.txt",
    "bisseq_pipeline_summary.R_WGBS_EpiG_IMR90_100_12_I12.txt"),
  sep="/")
names(summaryfiles) <- slabels

cpgfiles <- paste( datadir,
  c("cpg.IMR_90.mincov10.txt",
    "cpg.IMR90_A_TruSeq.mincov10.txt",
    "cpg.IMR90_1.mincov10.txt",
    "cpg.IMR90-_2.mincov10.txt",
    "cpg.IMR90_Roche_1ug.mincov10.txt",
    "cpg.IMR90_Roche_0_25ug.mincov10.txt",
    "cpg.R_WGBS_EpiG_IMR90_100_12_I12.mincov10.txt"),
  sep="/")

names(cpgfiles) <- slabels


cpgfiles0x <- paste( datadir,
  c("methylcall.CpG.IMR_90.mincov0.txt",
    "methylcall.CpG.IMR90_A_TruSeq.mincov0.txt",
    "methylcall.CpG.IMR90_1.mincov0.txt",
    "methylcall.CpG.IMR90-_2.mincov0.txt",
    "methylcall.CpG.IMR90_Roche_1ug.mincov0.txt",
    "methylcall.CpG.IMR90_Roche_0_25ug.mincov0.txt",
    "methylcall.CpG.R_WGBS_EpiG_IMR90_100_12_I12.mincov0.txt"),
  sep="/")
names(cpgfiles0x) <- slabels


bamfiles <- paste( datadir,
  c("IMR_90.sorted.bam",
  "IMR90_A_TruSeq.sorted.bam",
  "IMR90_1.sorted.bam",
  "IMR90-_2.sorted.bam",
  "IMR90_Roche_1ug.sorted.bam",
  "IMR90_Roche_0_25ug.sorted.bam",
  "R_WGBS_EpiG_IMR90_100_12_I12.sorted.bam"),
  sep="/")
names(bamfiles) <- slabels


regionFiles <- paste(datadir,
  c("MethylSeq.84Mb.bed",
    "130912_HG19_CpGiant_4M_EPI_noRandom.bed", # was 130912_HG19_CpGiant_4M_EPI.bed
    "MspI-fragments-hg19-20140919.bed"),
  sep="/")
names(regionFiles) <- c("SSMethylSeq","CpGiant","MspI_84-334")


# regionFiles <- paste(datadir,
  # c("130912_HG19_CpGiant_4M_EPI_noRandom.bed", # was 130912_HG19_CpGiant_4M_EPI.bed
    # "MethylSeq.84Mb.bed",
    # "MspI-fragments-hg19-20140919.bed",
    # "imr90frags.bed",
    # "imr90Afrags.bed",
    # "R_WGBS_EpiG_IMR90Afrags.bed",
    # "R_WGBS_EpiG_IMR90Bfrags.bed"),
  # sep="/")
# names(regionFiles) <- c("CpGiant","SSMethylSeq","MspI_84-334","ERRBS_A", "ERRBS_B", "WGBS_A", "WGBS_B")


########################
#                      #
#  On-target coverage  #
#                      #
########################



plotBinnedRegionsRecovered <- function(mcgr, regiongr) {
  # require(GenomicRanges)
  label <- paste("Recovered", names(regiongr), "regions for", names(mcgr))
  x <- findOverlaps(regiongr[[1]], mcgr[[1]])
  xx <- regiongr[[1]][unique(queryHits(x))]
  x <- width(xx)
  recovbins <- c(
                length(x[x>0 & x<50]),
                length(x[x>=50 & x<80]),
                length(x[x>=80 & x<100]),
                length(x[x>=100 & x<200]),
                length(x[x>=200 & x<400]),
                length(x[x>=400 & x<800]),
                length(x[x>=800 & x<1000]),
                length(x[x>=1000 & x<2000]),
                length(x[x>=2000 & x<4000]),
                length(x[x>=4000 & x<5000]),
                length(x[x>5000])
                )
  names(recovbins) <- c("1-49","50-79","80-99","100-199","200-399",
                        "400-799","800-999","1K-2K","2K-4K","4K-5K",
                        ">5K")

  x <- width(regiongr[[1]])
  allbins <- c(
                length(x[x>0 & x<50]),
                length(x[x>=50 & x<80]),
                length(x[x>=80 & x<100]),
                length(x[x>=100 & x<200]),
                length(x[x>=200 & x<400]),
                length(x[x>=400 & x<800]),
                length(x[x>=800 & x<1000]),
                length(x[x>=1000 & x<2000]),
                length(x[x>=2000 & x<4000]),
                length(x[x>=4000 & x<5000]),
                length(x[x>5000])
                )
  names(allbins) <- c("1-49","50-79","80-99","100-199","200-399",
                        "400-799","800-999","1K-2K","2K-4K","4K-5K",
                        ">5K")

  bins <- rbind(recovbins,allbins-recovbins)
  rownames(bins) <- c("recovered", "uncovered")

  acol <- c("blue","gray90")
  yaxis <- pretty(range(c(0,apply(bins,2,sum))))
  #pbin  <- round(bins/length(x)*100,digits=1)
  cex   <- 1
  par(mar=c(6,4,4,4))
  tmp <- barplot(bins,
                col=acol,
                ylim=c(0,tail(yaxis,n=1)),
                cex.axis=cex,
                cex.names=cex,
                xlab="",
                ylab="counts",
                main=label,
                las=2,
                axes=F)
  axis(2,at=yaxis,labels=paste(yaxis/1000,"K",sep=""),las=1,cex.axis=cex)
  # text(tmp,max(yaxis)-(max(yaxis)*0.05),labels=bins,cex=cex)
  # text(tmp,max(yaxis)-(max(yaxis)*0.09),labels=paste(pbin,"%"),cex=cex)
  # title(sub=paste("number of ", i, " regions above bar"),cex.sub=cex)
  # mtext(paste("total ", i, " regions",length(x)),adj=0,cex=cex)
  title(xlab="bin range (bp)", line=5)
  legend("topright",rownames(bins),fill=acol)
  par(mar=c(5.1,4.1,4.1,2.1))
}

plotBinnedRegionsRecovered2 <- function(mcgr, regiongr) {
  # require(GenomicRanges)
  label <- paste("Recovered", names(regiongr), "regions for", names(mcgr))
  x <- findOverlaps(regiongr[[1]], mcgr[[1]])
  xx <- regiongr[[1]][unique(queryHits(x))]
  x <- width(xx)
  recovbins <- c(length(x[x>0 & x<50]),
                length(x[x>=50 & x<84]),
                length(x[x>=84 & x<=334]),
                length(x[x>334 & x<=500]),
                length(x[x>500 & x<=1000]),
                length(x[x>1000 & x<=10000]),
                length(x[x>10000]))
  names(recovbins) <- c("1-49","50-83","84-334","335-500","501-1000","1k-10K",">10K")

  x <- width(regiongr[[1]])
  allbins <- c(length(x[x>0 & x<50]),
              length(x[x>=50 & x<84]),
              length(x[x>=84 & x<=334]),
              length(x[x>334 & x<=500]),
              length(x[x>500 & x<=1000]),
              length(x[x>1000 & x<=10000]),
              length(x[x>10000]))
  names(allbins) <- c("1-49","50-83","84-334","335-500","501-1000","1k-10K",">10K")

  bins <- rbind(recovbins,allbins-recovbins)
  rownames(bins) <- c("recovered", "uncovered")

  acol <- c("gray30","gray90")
  yaxis <- pretty(range(apply(bins,2,sum)))
  #pbin  <- round(bins/length(x)*100,digits=1)
  cex   <- 0.6
  par(mar=c(7,4,4,4))
  tmp <- barplot(bins,
                col=acol,
                ylim=c(0,tail(yaxis,n=1)),
                cex.axis=cex,
                cex.names=cex,
                xlab="bin range (bp)",
                ylab="counts",
                main=label,
                las=1,
                axes=F)
  axis(2,at=yaxis,labels=paste(yaxis/1000,"K",sep=""),las=1,cex.axis=cex)
  # text(tmp,max(yaxis)-(max(yaxis)*0.05),labels=bins,cex=cex)
  # text(tmp,max(yaxis)-(max(yaxis)*0.09),labels=paste(pbin,"%"),cex=cex)
  # title(sub=paste("number of ", i, " regions above bar"),cex.sub=cex)
  # mtext(paste("total ", i, " regions",length(x)),adj=0,cex=cex)
  legend("topleft",rownames(bins),fill=acol)
  par(mar=c(5.1,4.1,4.1,2.1))
}



plotRegionsRecovered <- function(mcgrlist, regiongrlist) {
  # require(GenomicRanges)
  regrecovered <- numeric()
  w <- width(regiongrlist[[1]])

  x <- findOverlaps(regiongrlist[[1]], mcgrlist[[1]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[1]]) * 100)

  ssregions <- regiongrlist[[1]][which(w>=50 & w<=500)]
  x <- findOverlaps(ssregions, mcgrlist[[1]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(ssregions) * 100)

  ssregions <- regiongrlist[[1]][which(w>=84 & w<=334)]
  x <- findOverlaps(ssregions, mcgrlist[[1]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(ssregions) * 100)

  x <- findOverlaps(regiongrlist[[1]], mcgrlist[[2]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[1]]) * 100)

  ssregions <- regiongrlist[[1]][which(w>=50 & w<=500)]
  x <- findOverlaps(ssregions, mcgrlist[[2]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(ssregions) * 100)

  ssregions <- regiongrlist[[1]][which(w>=84 & w<=334)]
  x <- findOverlaps(ssregions, mcgrlist[[2]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(ssregions) * 100)

  x <- findOverlaps(regiongrlist[[2]], mcgrlist[[3]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[2]]) * 100)

  x <- findOverlaps(regiongrlist[[2]], mcgrlist[[4]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[2]]) * 100)

  x <- findOverlaps(regiongrlist[[3]], mcgrlist[[5]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[3]]) * 100)

  x <- findOverlaps(regiongrlist[[3]], mcgrlist[[6]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[3]]) * 100)

  names(regrecovered) <- c(paste(names(mcgrlist)[1],"full", names(regiongrlist)[1]),
                           paste(names(mcgrlist)[1], names(regiongrlist)[1], "50-500bp"),
                           paste(names(mcgrlist)[1], names(regiongrlist)[1], "84-334bp"),
                           paste(names(mcgrlist)[2],"full", names(regiongrlist)[1]),
                           paste(names(mcgrlist)[2], names(regiongrlist)[1], "50-500bp"),
                           paste(names(mcgrlist)[2], names(regiongrlist)[1], "84-334bp"),                           
                           paste(names(mcgrlist)[3], names(regiongrlist)[2]),
                           paste(names(mcgrlist)[4], names(regiongrlist)[2]),
                           paste(names(mcgrlist)[5], names(regiongrlist)[3]),
                           paste(names(mcgrlist)[6], names(regiongrlist)[3]))

  par(mar=c(12,4,4,2))
  b <- barplot(regrecovered,
  main="Region recovery",
  ylab="percent (%)",
  xlab="",
  ylim=c(0,100),
  las=2)
  #title(xlab="samples and region datasets",line=10)
  text(b,regrecovered,round(regrecovered,digits=1),pos=1)
  par(mar=c(5.1,4.1,4.1,2.1))
}

# end On-target coverage

########################
#                      #
#  Fragment size       #
#                      #
########################

getBamFragmentLengths <- function(afile) {
  #require(Rsamtools)
  param <- ScanBamParam(what=c('isize'))
  isize <- scanBam(afile, param=param, use.name=T)
  return (isize)
}

plotFragmentSizeHist <- function(isize, aname) {
  plot(density(isize),
    main=paste(aname,"fragment size density"),
    xlab="fragment size (bp)")
}
plotFragmentSizeDensity <- function(isize, aname) {
  hist(isize,
    main=paste(aname,"fragment size histogram"),
    xlab="fragment size (bp)")
}

# end Fragment size



tryme <- function() {

# source("/home/thk2008/bin/correlation2.EC-AA-2045.utils.R")
  source("/home/thk2008/bin/corErrbsXt.utils.R")
  datadir <- "/scratch001/thk2008_dat/methylomecapture/data"
  slabels <- c("ERRBS_A","ERRBS_B","SSMethylSeq_A_opt","SSMethylSeq_B_min","CpGiant_A_opt","CpGiant_B_min","WGBS")

  cpgfiles <- paste( datadir,
    c("cpg.IMR_90.mincov10.txt",
      "cpg.IMR90_A_TruSeq.mincov10.txt",
      "cpg.IMR90_1.mincov10.txt",
      "cpg.IMR90-_2.mincov10.txt",
      "cpg.IMR90_Roche_1ug.mincov10.txt",
      "cpg.IMR90_Roche_0_25ug.mincov10.txt",
      "cpg.R_WGBS_EpiG_IMR90_100_12_I12.mincov10.txt"),
    sep="/")
  names(cpgfiles) <- slabels
  
  makepdf <- T
  atime  <- format(Sys.time(), "%Y%m%d")
  if (makepdf) {
    outdir <- paste("results",atime, sep="-")
    if (!file.exists(outdir)) dir.create(outdir)
    cat("output directory is [", outdir, "]\n")
  }

  colClasses <- c("character","NULL","NULL","NULL","integer","numeric","NULL")
  mclist <- list()
  for (i in names(cpgfiles)) {
    #add to test: nrows=10000)
    mclist[[i]] <- read.table(cpgfiles[i], header=T, colClasses=colClasses) 
  }
 
  pdf(paste(outdir, "pwcorrelation.pdf", sep="/"), width=10, height=10)
   plotPWcor(mclist)
  dev.off()

  
  
  mat <- numeric()
  for (i in 1:(length(mclist)-1)) {
    for (j in (i+1):length(mclist)) {
      mat <- rbind(mat,c(i,j))
    }
  }

  apply(mat,1,function(i) {
    x <- mclist[[i[1]]][,3]
    names(x) <- mclist[[i[1]]][,1]
    y <- mclist[[i[2]]][,3]
    names(y) <- mclist[[i[2]]][,1]
    z <- intersect(names(x),names(y))
    cat(paste("r =",round(cor(x[z],y[z]),digits=2)),"\n")
  })

  pdf(paste(outdir, "allIntersectCor.pdf", sep="/"), width=10, height=10)
    plotAllItersectCor(mclist)
  dev.off()
  
}





















