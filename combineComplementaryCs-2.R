
# combines complementary CpG's into single value

if (F) {
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/combineComplementaryCs-2.R")

}

require(GenomicRanges)
require("BSgenome.Hsapiens.UCSC.hg19")

createmcgr     <- T
createmcgrlist <- F


if(createmcgr) {
  
  # datafiles <- list.files("data",pattern="methylcall.CpG.*.mincov0.txt",full=T)

  slabels <- c("ERRBS_A","ERRBS_B","SSMethylSeq_A_opt","SSMethylSeq_B_min","CpGiant_A_opt","CpGiant_B_min","WGBS")
  datadir <- "/scratch001/thk2008_dat/methylomecapture/data"
  datafiles <- paste( datadir,
    c("methylcall.CpG.IMR_90.mincov0.txt",
      "methylcall.CpG.IMR90_A_TruSeq.mincov0.txt",
      "methylcall.CpG.IMR90_1.mincov0.txt",
      "methylcall.CpG.IMR90-_2.mincov0.txt",
      "methylcall.CpG.IMR90_Roche_1ug.mincov0.txt",
      "methylcall.CpG.IMR90_Roche_0_25ug.mincov0.txt",
      "methylcall.CpG.R_WGBS_EpiG_IMR90_100_12_I12.mincov0.txt"),
    sep="/")
  names(datafiles) <- slabels
  
  outdir <- "cgUnits10x"
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  
  colClasses <- c("character","character","numeric","character","numeric","numeric","numeric")
  
  for (afile in datafiles) {
    sname <- unlist(strsplit(afile,"\\."))[3]
    cat("reading file", afile, ", ")
    cpg.frame <- read.table(afile, colClasses=colClasses, header=T)
    rownames(cpg.frame) <- NULL
  
    cat("making granges, ")
    tmp <- rep("+", nrow(cpg.frame))
    tmp[cpg.frame[,4] == "R"] <- "-"
    cpg.frame[,4] <- as.factor(tmp) # convert all F/R to +/-
    cpg.frame[,6] <- round(cpg.frame[,5]*cpg.frame[,6]/100) # numCs
    cpg.frame[,7] <- round(cpg.frame[,5]*cpg.frame[,7]/100) # numTs
    tmp.gr <- GRanges(seqnames=cpg.frame[,2],
                      ranges=IRanges(cpg.frame[,3],end=cpg.frame[,3]),
                      strand=cpg.frame[,4],
                      coverage=cpg.frame[,5],
                      numC=cpg.frame[,6],
                      numT=cpg.frame[,7])
    chr.len <- seqlengths(Hsapiens)
    chr.len <- chr.len[grep("_", names(chr.len), invert = T)]
    seqlengths(tmp.gr) <- chr.len[names(seqlengths(tmp.gr))]
  
    cat("combining duplicates, ")
    dupgr <- which(duplicated(tmp.gr))
    o <- findOverlaps(tmp.gr[dupgr], tmp.gr)
    cdups <- GRanges()
    n <- length(unique(queryHits(o)))
    for (i in unique(queryHits(o))) {
      if ( ((i / n) %% 0.01) == 0) {
        cat(i/n*100,"%\n")
      }
      x <- tmp.gr[subjectHits(o[queryHits(o) == i])]
      y <- unique(x)
      values(y) <- lapply(values(x),sum)
      cdups <- c(cdups, y)
    }
    gr2 <- c(tmp.gr[-subjectHits(o)], cdups)
    gr2 <- sort(gr2)
    
    cat("combining sites, ")
    mcgr      <- gr2  
    mcgrplus  <- mcgr[strand(mcgr) == "+"]
    mcgrminus <- mcgr[strand(mcgr) == "-"]
    start(mcgrminus)  <- start(mcgrminus)-1
    end(mcgrminus)    <- end(mcgrminus)-1
    strand(mcgrminus) <- "+"
    x <- findOverlaps(mcgrplus,mcgrminus)
    mcgrplus <- c(mcgrplus, mcgrminus[-subjectHits(x)])
    mcgrplus[queryHits(x)]$coverage <- mcgrplus[queryHits(x)]$coverage + mcgrminus[subjectHits(x)]$coverage
    mcgrplus[queryHits(x)]$numC     <- mcgrplus[queryHits(x)]$numC + mcgrminus[subjectHits(x)]$numC
    mcgrplus[queryHits(x)]$numT     <- mcgrplus[queryHits(x)]$numT + mcgrminus[subjectHits(x)]$numT
    mcgr <- mcgrplus
    mcgr10x <- mcgr[which(mcgr$coverage >= 10)]
    strand(mcgr10x) <- "*"
    mcgr10x$freqC <- round(mcgr10x$numC / mcgr10x$coverage * 100, digits=3)
    mcgr10x$freqT <- round(mcgr10x$numT / mcgr10x$coverage * 100, digits=3)
    
    save(mcgr10x, file=paste(outdir, paste(sname,"cgUnits10xgr.rda",sep="_"), sep="/"))
  
    cat("done. saved",paste(sname,"cgUnits10xgr.rda",sep="_"),"\n")
  }

}

if (createmcgrlist) {
  
  require(GenomicRanges)
  # datafiles <- list.files("cgUnits10x", pattern="cgUnits10xgr.rda",full=T)
  slabels <- c("ERRBS_A","ERRBS_B","SSMethylSeq_A_opt","SSMethylSeq_B_min","CpGiant_A_opt","CpGiant_B_min","WGBS")
  datadir <- "/scratch001/thk2008_dat/methylomecapture/cgUnits10x"
  datafiles <- paste( datadir,
    c("IMR_90_cgUnits10xgr.rda", 
      "IMR90_A_TruSeq_cgUnits10xgr.rda",
      "IMR90_1_cgUnits10xgr.rda",
      "IMR90-_2_cgUnits10xgr.rda",
      "IMR90_Roche_1ug_cgUnits10xgr.rda",
      "IMR90_Roche_0_25ug_cgUnits10xgr.rda",
      "R_WGBS_EpiG_IMR90_100_12_I12_cgUnits10xgr.rda"),
    sep="/")
  names(datafiles) <- slabels

  for (i in 1:length(datafiles)) {
    if (file.exists(datafiles[i])) { 
      cat(basename(datafiles[i]), "exists\n")
    } else {
      stop("file does not exist [",datafiles[i],"]") 
    }
  }  
  
  mcgrlist <- list()
  for (i in 1:length(datafiles)) {
    cat(i,"\n")
    load(datafiles[i])
    mcgrlist[[i]] <- mcgr10x
    rm(mcgr10x)
  }
  names(mcgrlist) <- slabels
  save(mcgrlist, file="mcgrlist_CGunits10x.rda")
  
}




if(F){
  for (i in 1:length(mcgrlist)) {
    cat(i,"\n")
    grA <- mcgrlist[[i]]  
    grA$freqC <- round(grA$numC / grA$coverage * 100,digits=3)
    grA$freqT <- round(grA$numT / grA$coverage * 100,digits=3)
    mcgrlist[[i]] <- grA
  }
  save(mcgrlist, file="mcgrlist_CGunits10x.rda")
}



