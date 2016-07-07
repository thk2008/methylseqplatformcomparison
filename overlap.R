# overlap

# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/overlap.R")

require(GenomicRanges)
#require(corrplot)



getOverlapMat <- function(mcgrlist) {
  snames         <- names(mcgrlist)
  overlapMat     <- numeric()
  rnames         <- character()
  overlapMatPerc <- matrix(NA,length(snames),length(snames),dimnames=list(snames,snames))
  diag(overlapMatPerc) <- 100
  for (i in 1:(length(mcgrlist)-1)) {
    for (j in (i+1):length(mcgrlist)) {
      cat(names(mcgrlist)[i],"::",names(mcgrlist)[j],"\n")
      rnames <- c(rnames, paste(names(mcgrlist)[i],"::",names(mcgrlist)[j],sep=""))
      numi   <- length(mcgrlist[[i]])
      numj   <- length(mcgrlist[[j]])
      overlapCounts       <- sum(countOverlaps(mcgrlist[[i]], mcgrlist[[j]]))
      overlapMat          <- rbind(overlapMat, c(numi, numj, overlapCounts))
      overlapMatPerc[i,j] <- overlapCounts / numi * 100
      overlapMatPerc[j,i] <- overlapCounts / numj * 100
    }
  }
  rownames(overlapMat) <- rnames
  colnames(overlapMat) <- c("x::","::x","overlap")

  return (list(overlapmat=overlapMat,percoverlap=overlapMatPerc))
}


# looks like this:
# x <- overlaps$overlapmat[1,]
#     x::     ::x overlap
# 6629668 6696652 5955415
plotCpGsiteOverlap <- function(sdat, slabels) {
  y <- rbind(sdat[3], sdat[1:2] - sdat[3])
  colnames(y) <- unlist( strsplit( slabels,"::" ))
  perc <- (sdat[3] / sdat[1:2]) * 100
  bcols <- c("slateblue","salmon")
  b <- barplot(y,
    main=paste("CpG site overlap", slabels),
    ylab="number of CpG sites",
    sub=paste("common sites: ",sdat[3]),
    col=bcols)
  text(b,y[1,],paste(round(perc,digits=1),"%"),pos=3)
  mtext(paste("total sites: ",sdat[1:2]), side=1, at=b, line=2)
}




# end Overlap
