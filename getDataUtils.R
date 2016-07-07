# get data functions

require(GenomicRanges)
require("BSgenome.Hsapiens.UCSC.hg19")

# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/getDataUtils.R")

# input  : summary files and names
# output : matrix of stats from files
getPropertiesMatrix <- function(sfiles) {
  propsmat <- numeric()
  for (afile in sfiles) {
  #afile <- summaryfiles[1]
    x         <- readLines(afile,n=40)
    n         <- grep("property\\s+value", x, perl=T)
    x         <- x[(n+2):(n+20)]
    props     <- sapply(strsplit(x," +"),function(x) {x[1]})
    propnames <- props[!is.na(props)] # these should all be identical
    vals      <- sapply(strsplit(x," +"),function(x) {x[2]})
    vals      <- as.numeric(vals[!is.na(vals)])
    propsmat  <- rbind(propsmat, vals)
  }
  colnames(propsmat) <- propnames
  rownames(propsmat) <- names(sfiles)
  return(propsmat)
}

savePropsCsv <- function(props, outdir=".") {
  afile <- paste(outdir, "propertysummary.csv",sep="/")
  cat(",",file=afile)
  ow <- options("warn")
  options(warn = -1)
  write.table(props,file=afile,append=T,quote=F,sep=",")
  options(ow) # reset
}


# input   : named cpg files character vector
# returns : named list of genomic ranges objects
getCpGgrList <- function(mcfiles) {
  mcgrlist <- list()
  for (i in names(mcfiles)) {
    cat("getting " ,i, "...\n")
    cpg.frame <- read.table(mcfiles[i], header=T)
    tmp <- rep("+", nrow(cpg.frame))
    tmp[cpg.frame[,4] == "R"] <- "-"
    cpg.frame[,4] <- as.factor(tmp)
    cpg.frame[,6] <- round(cpg.frame[,5]*cpg.frame[,6]/100) # numCs
    cpg.frame[,7] <- round(cpg.frame[,5]*cpg.frame[,7]/100) # numTs
    tmp.gr   <- GRanges(seqnames=cpg.frame[,2],ranges=IRanges(cpg.frame[,3],end=cpg.frame[,3]),strand=cpg.frame[,4],coverage=cpg.frame[,5],numC=cpg.frame[,6],numT=cpg.frame[,7])
    mcgr     <- unique(tmp.gr)
    chr.len  <- seqlengths(Hsapiens)  # get chromosome lengths
    chr.len  <- chr.len[grep("_", names(chr.len), invert = T)] # remove random chromosomes
    seqlengths(mcgr) <- chr.len[names(seqlengths(mcgr))]
    mcgrlist[[i]] <- mcgr
  }
  return (mcgrlist)
}


getCpgGr <- function(mcfile) {

  cpg.frame <- read.table(mcfile, header=T)
  tmp <- rep("+", nrow(cpg.frame))
  tmp[cpg.frame[,4] == "R"] <- "-"
  cpg.frame[,4] <- as.factor(tmp)
  # cpg.frame[,6] <- round(cpg.frame[,5]*cpg.frame[,6]/100) # numCs
  # cpg.frame[,7] <- round(cpg.frame[,5]*cpg.frame[,7]/100) # numTs
  tmp.gr   <- GRanges(seqnames=cpg.frame[,2],ranges=IRanges(cpg.frame[,3],end=cpg.frame[,3]),strand=cpg.frame[,4],coverage=cpg.frame[,5],numC=cpg.frame[,6],numT=cpg.frame[,7])
  mcgr     <- unique(tmp.gr)
  chr.len  <- seqlengths(Hsapiens)  # get chromosome lengths
  chr.len  <- chr.len[grep("_", names(chr.len), invert = T)] # remove random chromosomes
  seqlengths(mcgr) <- chr.len[names(seqlengths(mcgr))]
  
  return (mcgr)
}


getGrRegion <- function(regionfile, header=F) {
  bedframe <- read.table(regionfile, header=header, fill=T)
  regiongr <- GRanges(seqnames=bedframe[,1],ranges=IRanges(bedframe[,2],end=bedframe[,3]))
  chr.len  <- seqlengths(Hsapiens)  # get chromosome lengths
  chr.len  <- chr.len[grep("_", names(chr.len), invert = T)] # remove random chromosomes
  seqlengths(regiongr) <- chr.len[names(seqlengths(regiongr))]
  return (regiongr)
}


getRegionsList <- function(regionFiles) {
  slabels <- names(regionFiles)
  regionsList <- list()
  for (i in 1:length(regionFiles)) {    
    if (i == 1 || i == 2) {
      cat("getting", slabels[i], "\n")
      regionsList[ slabels[i] ] <- getGrRegion(regionFiles[i])
    }    
    if (i == 3) {
      cat("getting", slabels[i], "\n")
      regionsList[ slabels[i] ] <- getGrRegion(regionFiles[i])
      x <- width(regionsList[[i]])
      x <- which( (x >= 84) & (x <= 334) )
      mspi84.334 <- regionsList[[i]][x]
      regionsList[[i]] <- mspi84.334
      rm(x, mspi84.334)
    }
    if (i == 4) {
      cat("getting", slabels[i], "\n")
      x1 <- getGrRegion(regionFiles[i], header=T)
      #regionsList[ slabels[i] ] <- NULL
    }
    if (i == 5) {
      cat("getting", slabels[i], "\n")
      x2 <- getGrRegion(regionFiles[i], header=T)
      cat("union", slabels[i], "\n")
      regionsList[ slabels[i] ] <- union(x1,x2)
      rm(x1,x2)
    }
    if (i == 6) {
      cat("getting", slabels[i], "\n")
      x1 <- getGrRegion(regionFiles[i], header=T)
      # regionsList[ slabels[i] ] <- getGrRegion(regionFiles[i], header=T)
    }
    if (i == 7) {
      cat("getting", slabels[i], "\n")
      x2 <- getGrRegion(regionFiles[i], header=T)
      cat("union", slabels[i], "\n")
      regionsList[ slabels[i] ] <- union(x1,x2)
      rm(x1,x2)
    }
  }
  
  chr.len  <- seqlengths(Hsapiens)  # get chromosome lengths
  chr.len  <- chr.len[grep("_", names(chr.len), invert = T)] # remove random chromosomes
  for (i in 1:length(regionsList)) {
    seqlengths(regionsList[[i]]) <- chr.len[names(seqlengths(regionsList[[i]]))]
  }

  return(regionsList)
}
