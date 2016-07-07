
if(F){
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/annotations.R")
}

require(GenomicRanges)
require(methylKit)
require("BSgenome.Hsapiens.UCSC.hg19")
require(vioplot)
library(beeswarm)

getAnnotsObj <- function(cpgfiles) {

  annotsobj <- read (as.list(cpgfiles), sample.id=as.list(names(cpgfiles)),assembly="hg19",treatment=rep(1,length(cpgfiles)))

  annotsobj <- lapply(annotsobj, as, Class="GRanges")
  names(annotsobj) <- names(cpgfiles)
  chr.len  <- seqlengths(Hsapiens)  # get chromosome lengths
  chr.len  <- chr.len[grep("_", names(chr.len), invert = T)] # remove random chromosomes
  for (i in 1:length(annotsobj)) {
    seqlengths(annotsobj[[i]]) <- chr.len[names(seqlengths(annotsobj[[i]]))]
  }
  return(annotsobj)
}


getAnnotStats <- function(annotsobj, gene.obj, cpg.obj) {

  annotStats <- list()

  for (i in names(annotsobj)) {
    cat(i)
    tmp <- list()
    gAnnots.i         <- annotate.WithGenicParts(annotsobj[[i]], gene.obj)
    fAnnots.i         <- annotate.WithFeature.Flank(annotsobj[[i]],cpg.obj$CpGi,cpg.obj$shores, feature.name="CpGi",flank.name="shores")
    # plotTargetAnnotation(gAnnots.i, precedence=T, main=paste(mclabels[i], "methylation annotation"))
    # plotTargetAnnotation(fAnnots.i,col=c("green","gray","white"),main=paste(mclabels[i],"methylation annotation"))
    tmp <- list("numberGpart"   = getTargetAnnotationStats(gAnnots.i,percentage=F),
                "percGpart"     = getTargetAnnotationStats(gAnnots.i),
                "numberFeature" = getTargetAnnotationStats(fAnnots.i,percentage=F),
                "percFeature"   = getTargetAnnotationStats(fAnnots.i))
    gmat <- numeric()
    fmat <- numeric()
    for (j in names(annotsobj)) {
      cat(" ",j)
      k <- subsetByOverlaps(annotsobj[[i]], annotsobj[[j]])
      gAnnots.k                <- annotate.WithGenicParts(k, gene.obj)
      overlappingSitesGpart    <- getTargetAnnotationStats(gAnnots.k, percentage=F)
      fAnnots.k                <- annotate.WithFeature.Flank(k, cpg.obj$CpGi, cpg.obj$shores, feature.name="CpGi",flank.name="shores")
      overlappingSitesFeature  <- getTargetAnnotationStats(fAnnots.k, percentage=F)
      gmat <- rbind(gmat, overlappingSitesGpart)
      fmat <- rbind(fmat, overlappingSitesFeature)
    }
    rownames(gmat) <- names(annotsobj)
    rownames(fmat) <- names(annotsobj)
    tmp <- c(tmp, list("numOverlapGpart" = gmat, "numOverlapFeature" = fmat, "numCpGs" = length(annotsobj[[i]])))
    annotStats[[i]] <- tmp
    cat(".\n")
  }
  return(annotStats)

}


plotPercentGenePart <- function(annotStats) {

  acols <- c("red","green","cyan","purple")

  x <- sapply(annotStats, "[[", "percGpart")
  x <- cbind(x[,1:2],NA,x[,3:4],NA,x[,5:6],NA,x[,7:8])
  names.arg <- c(rep(c("A","B",""),3),"A","B")
  platform  <- c("ERRBS", "SSMethylSeq", "CpGiant", "WGBS")

  b <- barplot(x,
    space=0,
    ylab="percent (%)",
    col=acols,
    names.arg=names.arg,
    legend=rownames(x),
    border=rgb(0.5,0.5,0.5,alpha=0.6),
    args.legend=list(x=9,y=8, border=rgb(0.5,0.5,0.5,alpha=0.6), bg=rgb(1,1,1,alpha=0.5),horiz = T,box.lwd=0.1),
    axes=F)
  axis(2,seq(0,100,by=10),las=2)
  mtext(platform, side=1, at=seq(1,length(b),by=3), line=2) # platform labels
  title(main="Proportion of CpGs annotated by gene parts")

}


plotPercentFeature <- function(annotStats) {

  acols <- c("lightgreen","gray","white")

  x <- sapply(annotStats, "[[", "percFeature")
  x <- cbind(x[,1:2],NA,x[,3:4],NA,x[,5:6],NA,x[,7:8])
  names.arg <- c(rep(c("A","B",""),3),"A","B")
  platform  <- c("ERRBS", "SSMethylSeq", "CpGiant", "WGBS")

  b <- barplot(x,
    space=0,
    ylab="percent (%)",
    col=acols,
    names.arg=names.arg,
    legend=rownames(x),
    border=rgb(0.5,0.5,0.5,alpha=0.6),
    args.legend=list(x=7,y=8, border=rgb(0.5,0.5,0.5,alpha=0.6), bg=rgb(1,1,1,alpha=0.5),horiz = T,box.lwd=0.1),
    axes=F)
  axis(2,seq(0,100,by=10),las=2)
  mtext(platform, side=1, at=seq(1,length(b),by=3), line=2) # platform labels
  title(main="Proportion of CpGs annotated by CpG feature")

}


plotAnnotsAbsoluteNumbers <- function(annotStats) {

  x <- rbind(sapply(annotStats, "[[", "numberGpart"),sapply(annotStats, "[[", "numberFeature"))
  xaxis <- pretty(range(x))

  acols <- c("red", "green", "cyan", "purple", "lightgreen", "gray", "white")
  oldmar <- par()$mar
  par(mar=c(5,7,4,2))
  b <- barplot(cbind(x),
    beside=T,
    horiz=T,
    xlim=c(0,tail(pretty(range(x)),n=1)),
    xlab="counts",
    ylab="",
    col=acols,
    border=rgb(0.4,0.4,0.4,alpha=0.7),
    las=1,
    legend=rownames(x),
    args.legend=list(x="right",border=rgb(0.5,0.5,0.5,alpha=0.6), bg="white",box.lwd=0.1),
    axes=F)

  segments(xaxis[-1],0,xaxis[-1],sum(max(b),min(b)),col=rgb(0.5,0.5,0.5,alpha=0.2))
  axis(1,at=xaxis)
  mtext(paste("(",unlist(lapply(annotStats, "[[", "numCpGs")),")",sep=""), side=2, 1, at=b[3,]-0.5, las=2)
  title(main="Number of CpGs by annotation (gene part/CpG feature)")
  title(ylab="(total CpGs)",line=6)

  par(mar=oldmar)
}



vioplot2 <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
    horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
    lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
    at, add = FALSE, wex = 1, drawRect = TRUE) 
{
    require(vioplot)
    datas <- list(x, ...)
    n <- length(datas)
    if (missing(at)) 
        at <- 1:n
    upper <- vector(mode = "numeric", length = n)
    lower <- vector(mode = "numeric", length = n)
    q1 <- vector(mode = "numeric", length = n)
    q3 <- vector(mode = "numeric", length = n)
    med <- vector(mode = "numeric", length = n)
    base <- vector(mode = "list", length = n)
    height <- vector(mode = "list", length = n)
    baserange <- c(Inf, -Inf)
    args <- list(display = "none")
    if (!(is.null(h))) 
        args <- c(args, h = h)
    for (i in 1:n) {
        data <- datas[[i]]
        data.min <- min(data)
        data.max <- max(data)
        q1[i] <- quantile(data, 0.25)
        q3[i] <- quantile(data, 0.75)
        med[i] <- median(data)
        iqd <- q3[i] - q1[i]
        upper[i] <- min(q3[i] + range * iqd, data.max)
        lower[i] <- max(q1[i] - range * iqd, data.min)
        est.xlim <- c(min(lower[i], data.min), max(upper[i], 
            data.max))
        smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
            args))
        hscale <- 0.4/max(smout$estimate) * wex
        base[[i]] <- smout$eval.points
        height[[i]] <- smout$estimate * hscale
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1], t[1])
        baserange[2] <- max(baserange[2], t[2])
    }
    if (!add) {
        xlim <- if (n == 1) 
            at + c(-0.5, 0.5)
        else range(at) + min(diff(at))/2 * c(-1, 1)
        if (is.null(ylim)) {
            ylim <- baserange
        }
    }
    if (is.null(names)) {
        label <- 1:n
    }
    else {
        label <- names
    }
    boxwidth <- 0.05 * wex
    if (!add)
        plot.new()
    if (!horizontal) {
        if (!add) {
            plot.window(xlim = xlim, ylim = ylim)
            axis(2)
            axis(1, at = at, label = label)
        }
        #box()
        for (i in 1:n) {
            polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])),
                c(base[[i]], rev(base[[i]])), col = col, border = border,
                lty = lty, lwd = lwd)
            if (drawRect) {
                lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
                  lty = lty)
                rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2,
                  q3[i], col = rectCol, border = rectCol)
                points(at[i], med[i], pch = pchMed, col = colMed)
            }
        }
    }
    else {
        if (!add) {
            plot.window(xlim = ylim, ylim = xlim)
            axis(1)
            axis(2, at = at, label = label)
        }
        #box()
        for (i in 1:n) {
            polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]],
                rev(at[i] + height[[i]])), col = col, border = border,
                lty = lty, lwd = lwd)
            if (drawRect) {
                lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd,
                  lty = lty)
                rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] +
                  boxwidth/2, col = rectCol, border = rectCol)
                points(med[i], at[i], pch = pchMed, col = colMed)
            }
        }
    }
    box()
    invisible(list(upper = upper, lower = lower, median = med,
        q1 = q1, q3 = q3))
}


plotAnnotPercCoverageViolin <- function(cpgcovmat, totC, what="?what?",h=5) {
  
  a                 <- cpgcovmat
  a[ a == 0 ]       <- NA
  totC[ totC == 0 ] <- NA
  b                 <- a / totC * 100
  omar = par()$mar
  par(mar=c(7,4,4,2))
  plot(0, type="n",
    main=paste("Percent coverage of CpGs in", what),
    xlab="",
    ylab="percent (%)",
    xlim=c(0, ncol(cpgcovmat)+1), ylim=c(0,100), axes=F)
  axis(2, las=2)
  axis(1, at=seq(1, ncol(cpgcovmat)), colnames(cpgcovmat), las=2)
  for (i in 1:ncol(b)) {
    vioplot2(b[!is.na(b[,i]),i],at=i,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  }
  # vioplot2(b[!is.na(b[,1]),1],at=1,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  # vioplot2(b[!is.na(b[,2]),2],at=2,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  # 
  # vioplot2(b[!is.na(b[,3]),3],at=3,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  # vioplot2(b[!is.na(b[,4]),4],at=4,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  # 
  # vioplot2(b[!is.na(b[,5]),5],at=5,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  # vioplot2(b[!is.na(b[,6]),6],at=6,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  # 
  # vioplot2(b[!is.na(b[,7]),7],at=7,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  # vioplot2(b[!is.na(b[,8]),8],at=8,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  par(mar=omar)
}

getAnnotCovMat <- function(mcgrlist, annotgr) {
  cpgcovmat <- matrix(0, nrow=length(annotgr), ncol=length(mcgrlist))
  colnames(cpgcovmat) <- names(mcgrlist)
  for (i in 1:length(mcgrlist)) {
    a <- findOverlaps(annotgr, mcgrlist[[i]])
    b <- table(queryHits(a))
    cpgcovmat[as.numeric(names(b)),i] <- as.numeric(b)
  }
  return (cpgcovmat)
}




plotAnnotCoverage <- function(x, what="?what?") {
  x         <- c(x[1:2],NA,x[3:4],NA,x[5:6],NA,x[7:8])
  names.arg <- c(rep(c("A","B",""),3),"A","B")
  platform  <- c("ERRBS", "SSMethylSeq", "CpGiant", "WGBS")
  cols      <- c("lightblue","lightblue","white","salmon","salmon","white","plum2","plum2","white","lightgreen","lightgreen")
  b <- barplot(x,
    ylim=c(0,100),
    main=paste("Percent of", what, "covered by CpGs"),
    ylab="percent (%)",
    xlab="platform",
    sub="coverage >10x",
    names.arg=names.arg,
    col=cols,
    las=1,
    space=0)
  mtext(platform, side=1, at=seq(1,length(b),by=3), line=2) # platform labels
  text(b,x,round(x,digits=1),pos=1, cex=0.8) # value in bars
}





plotAnnotCoverageAll <- function(annotcoverage) { 
  cols <- c("blue1","blue3","salmon1","salmon3","plum","plum4","green1","green3")
  beeswarm(annotcoverage,
    labels=names(annotcoverage),
    main="Percent feature region CpGs covered",
    xlab="genomic feature",
    ylab="percent (%)",
    ylim=c(0,100),
    pwcol=rep(cols,5),
    pch=21,
    pwbg=rep(cols,5),
    cex=1.5)
  #boxplot(annotcoverage,ylim=c(0,100),boxwex=0.1, border=rgb(0.6,0.6,0.6,alpha=0.5),pars=list(medcol="black",medlwd=1),axes=F,add=T)
  legend("bottomright", names(annotcoverage[[1]]), fill=cols, border=cols, box.col=rgb(0.6,0.6,0.6,alpha=0.5))
}














getStrandParity <- function(mcgrlist) {
  require(GenomicRanges)
  strandTable <- sapply(mcgrlist, function(i) {
    table(as.vector(strand(i)))
  })
  totals <- sapply(mcgrlist, length)
  strandPerc <- t(strandTable) / totals * 100
  return (list(strandTable=strandTable,strandPerc=strandPerc))
}


plotStrandParity <- function(mat) {
  x <- t(mat)
  x <- cbind(x[,1:2],NA,x[,3:4],NA,x[,5:6],NA,x[,7])
  names.arg <- c(rep(c("A","B",""),3),"A")
  platform  <- c("ERRBS", "SSMethylSeq", "CpGiant", "WGBS")
  bcols <- c("hotpink","lightblue")

  b <- barplot(x,
    space=0,
    xlab="platform",
    ylab="percent (%)",
    col=bcols,
    names.arg=names.arg,
    legend=c("forward (+)","reverse (-)"),
    border=rgb(0.5,0.5,0.5,alpha=0.6),
    args.legend=list(x="bottom", inset=0.02, horiz=T, fill=bcols, border=rgb(0.5,0.5,0.5,alpha=0.6), bg=rgb(1,1,1,alpha=0.5),cex=0.9, box.lwd=0.1),
    axes=F)
  axis(2,seq(0,100,by=10),las=2)
<<<<<<< HEAD
  at <- seq(1,length(b),by=3)
  at[length(at)] <- b[length(b)]
  mtext(platform, side=1, at=at, line=2) # platform labels
=======
  mtext(platform, side=1, at=seq(1,length(b),by=3), line=2) # platform labels
>>>>>>> 0717e4c7192348ea949c0343b073fc4ab0e02b88
  title(main="Strand parity")
}

# if(F) {
 # pdf(paste(outdir, "strandparity.pdf", sep="/"))
  # plotStrandParity(strandparity$strandPerc)
 # dev.off()
# }





if(F) {

  
# unique gene.obj categories
# names(gene.obj)
# [1] "exons"     "introns"   "promoters" "TSSes"

x <- gene.obj$exons
y <- unique(x)
gene.obj$exons <- y

x <- gene.obj$introns
y <- unique(x)
gene.obj$introns <- y

x <- gene.obj$promoters
y <- unique(x)
gene.obj$promoters <- y

  
  
  
  
  
  
  
#  overlap of gene parts
exonCoverage <- lapply(mcgrlist, function(i) { # list of exon regions covered by (overlapping) the CpG sites
    subsetByOverlaps(gene.obj$exons, i)
})
percExon <- sapply(exonCoverage, length) / length(gene.obj$exons) * 100
plotAnnotCoverage(percExon, what="Exons")


# prototype
# x <- findOverlaps(exonCoverage[[1]],exonCoverage[[2]])
# c(length(unique(queryHits(x)))/length(exonCoverage[[1]]), length(unique(subjectHits(x)))/length(exonCoverage[[2]]))

exonoverlap <- numeric()
a <- character()
for (i in 1:(length(exonCoverage)-1)) {
  for (j in (i+1):length(exonCoverage)) {
    cat(i, ":", j, "\n")
    x <- findOverlaps(exonCoverage[[i]],exonCoverage[[j]])
    exonoverlap <- rbind(exonoverlap, c(length(unique(queryHits(x)))/length(exonCoverage[[i]]), length(unique(subjectHits(x)))/length(exonCoverage[[j]])))
    a <- c(a, paste(names(exonCoverage)[i], "::", names(exonCoverage)[j],sep=""))
  }
}
rownames(exonoverlap) <- a
save(exonoverlap, file="exonoverlap.rda")


outdir <- "./"
makepdf <- T

omar <- par()$mar # c(5.1, 4.1, 4.1, 2.1)
if (makepdf) { pdf(paste(outdir,"exonCoverageOverlap.pdf",sep="/")) } else { dev.new() }

par(mar=c(11,4,4,2))
plot(0,type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,85), ylim=c(0,100))
segments(0, c(10,30,50,70,90), 85, col=rgb(0.5,0.5,0.5,alpha=0.5))
barplot(t(exonoverlap*100),
  main="percent overlap of exons covered by platform (pairwise)",
  ylab="percent (%)",
  beside=T,
  ylim=c(0,100),
  las=2,
  col=c("lightblue","pink"),
  cex.names=0.9,
  axes=F,
  add=TRUE)
axis(2,at=seq(0,100,by=10),las=2)

if (makepdf) { dev.off() }

par(mar=omar)



# what fraction of region cpgs are covered by
x              <- exonCoverage[[1]]
grplus         <- x[which(strand(x) == "+")]
end(grplus)    <- end(grplus)+1
grminus        <- x[which(strand(x) == "-")]
start(grminus) <- start(grminus)-1
xplus          <- getSeq(Hsapiens, grplus)
xminus         <- reverseComplement(getSeq(Hsapiens, grminus))
grplus$numC    <- vcountPattern("CG", xplus, fixed=T)
grminus$numC   <- vcountPattern("CG", xminus, fixed=T)
end(grplus)    <- end(grplus)-1
start(grminus) <- start(grminus)+1
errbsagr       <- c(grplus,grminus)

x              <- gene.obj$exons
grplus         <- x[which(strand(x) == "+")]
end(grplus)    <- end(grplus)+1
grminus        <- x[which(strand(x) == "-")]
start(grminus) <- start(grminus)-1
xplus          <- getSeq(Hsapiens, grplus)
xminus         <- reverseComplement(getSeq(Hsapiens, grminus))
grplus$numC    <- vcountPattern("CG", xplus, fixed=T)
grminus$numC   <- vcountPattern("CG", xminus, fixed=T)
end(grplus)    <- end(grplus)-1
start(grminus) <- start(grminus)+1
exonsgr        <- c(grplus,grminus)





apply(exoncpgcovmat, 2, function(i) { length(which(i > 0)) }) / length(exonsgr) * 100  # identical to percExon


plotAnnotPercCoverageViolin(exoncpgcovmat, exonsgr$numC , what="Exon")





















x <- findOverlaps(mcgrlist[[1]], errbsagr)

# y <- subjectHits(x)
# yy <- as.numeric(table(y))
# errbsagr$numCcovered <- yy

#errbsagr$numCcovered <- as.numeric(table(subjectHits(x)))
numCcovered <- as.numeric(table(subjectHits(x)))
numC <- errbsagr$numC

diffC <- cbind(numC, numCcovered, numC - numCcovered)
if (length(which(diffC[,3] < 0) != 0)) {
  stop("some regions are covered more than there are CpGs\n")
}
# errbsagr[which(diffC[,3] < 0)]

perC <- numCcovered / numC
hist(perC,breaks=50)


#
 # errbsagr[which((numC - numCcovered) < 0)]
# GRanges object with 623 ranges and 3 metadata columns:
        # seqnames                 ranges strand   |     score         name
           # <Rle>              <IRanges>  <Rle>   | <numeric>  <character>
    # [1]     chr1 [ 21890534,  21890709]      +   |         4 NM_001177520
    # [2]     chr1 [ 21890534,  21890709]      +   |         6    NM_000478
    # [3]     chr1 [ 21890534,  21890709]      +   |         5 NM_001127501
    # [4]     chr1 [ 29616195,  29616224]      +   |        15    NM_005704
    # [5]     chr1 [156896994, 156897062]      +   |         6    NM_144702
    # ...      ...                    ...    ... ...       ...          ...
  # [619]    chr22   [38484755, 38484971]      -   |        10    NM_025045
  # [620]    chr22   [42537544, 42537685]      -   |         6    NR_002570
  # [621]    chr22   [50706249, 50706378]      -   |         2    NM_002751
  # [622]    chr22   [50665428, 50665506]      -   |         6    NM_020461
  # [623]    chr22   [50699596, 50699725]      -   |         2    NM_002969
#
#
# reverseComplement( getSeq(Hsapiens, "chr22", 38484755, 38484971) )
#
#






x <- findOverlaps(mcgrlist[[1]], errbsagr)


numCcovered <- as.numeric(table(subjectHits(x)))
numC <- errbsagr$numC

diffC <- cbind(numC, numCcovered, numC - numCcovered)
if (length(which(diffC[,3] < 0) != 0)) {
  stop("some regions are covered more than there are CpGs\n")
}


perC <- numCcovered / numC
hist(perC,breaks=50)















exonsgr <- getRegionsTotalCpGs(gene.obj$exons)

cpgcovmat <- matrix(0, nrow=length(exonsgr), ncol=length(mcgrlist))
colnames(cpgcovmat) <- names(mcgrlist)
for (i in 1:length(mcgrlist)) {
  a <- findOverlaps(exonsgr, mcgrlist[[i]])
  b <- table(queryHits(a))
  cpgcovmat[as.numeric(names(b)),i] <- as.numeric(b)
}






i=4
a                 <- cpgcovmat[,i]
numC <- exonsgr$numC

diffC <- cbind(numC, a, numC - a)
which(diffC[,3] < 0) 

a[ a == 0 ]       <- NA
numC[ numC == 0 ] <- NA
b                 <- a / numC * 100

hist(b)





h=5
plot(0,type="n",xlim=c(0,9),ylim=c(0,100),axes=F)
axis(2,las=2)
vioplot2(b[!is.na(b[,1]),1],at=1,add=T,h=h,col=NA, rectCol="gray", colMed="red", border="blue")
vioplot2(b[!is.na(b[,2]),2],at=2,add=T,h=h,col=NA, rectCol="gray", colMed="red", border="blue")

vioplot2(b[!is.na(b[,3]),3],at=3,add=T,h=h,col=NA, rectCol="gray", colMed="red", border="blue")
vioplot2(b[!is.na(b[,4]),4],at=4,add=T,h=h,col=NA, rectCol="gray", colMed="red", border="blue")

vioplot2(b[!is.na(b[,5]),5],at=5,add=T,h=h,col=NA, rectCol="gray", colMed="red", border="blue")
vioplot2(b[!is.na(b[,6]),6],at=6,add=T,h=h,col=NA, rectCol="gray", colMed="red", border="blue")

vioplot2(b[!is.na(b[,7]),7],at=7,add=T,h=h,col=NA, rectCol="gray", colMed="red", border="blue")
vioplot2(b[!is.na(b[,8]),8],at=8,add=T,h=h,col=NA, rectCol="gray", colMed="red", border="blue")












}





