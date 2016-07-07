# methylome sequencing region analysis

# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/regionAnalysis.R")

# for more infromation see:
#  /home/thk2008/bin/reformatBismarkBam.pl
#  /home/thk2008/bin/methylseqplatformcomparison/scripts/old/recreateWGBSfragments.R

require(GenomicRanges)
require(vioplot)
require("BSgenome.Hsapiens.UCSC.hg19")
library(beeswarm)

# w1 <- list()
# for (i in 1:length(regionsList)) {
  # w1[[i]] <- width(regionsList[[i]])
# }
# dev.new()
# hist(x,breaks=100)
# boxplot()

makeViolinPlot <- function(regionsList) {
  w2 <- list()
  for (i in 1:length(regionsList)) {
    x <- width(regionsList[[i]])
    w2[[i]] <- x[x<5000]
  }
  #aCols <- c("red","green","magenta","blue","purple")
  aCols <- rep("gray",5)
  h <- 50

  plot(0, xlim=c(0.5,5.5), ylim=c(0,1000),main="Regions Width distributions",xlab="platform",ylab="width (bp)",type="n",axes=F)
  box()
  axis(2)
  axis(1, at=1:5,labels=names(regionsList))
  vioplot(w2[[1]], at=1,add=T,h=h,col=aCols[1])
  vioplot(w2[[2]], at=2,add=T,h=h,col=aCols[2])
  vioplot(w2[[3]], at=3,add=T,h=h,col=aCols[3])
  vioplot(w2[[4]], at=4,add=T,h=h,col=aCols[4])
  vioplot(w2[[5]], at=5,add=T,h=h,col=aCols[5])
}




panel.overlap <- function(i, j, ymax) {
    labelA <- names(i)
    labelB <- names(j)
    i <- i[[1]]
    j <- j[[1]]
    totalcpgA <- length(i)
    totalcpgB <- length(j)

    # this is typically not symmetrical
    x <- findOverlaps(i, j)
    numoverlapA <- length(unique(queryHits(x)))
    numoverlapB <- length(unique(subjectHits(x)))

    x <- rbind(c(numoverlapA,           numoverlapB),
               c(totalcpgA-numoverlapA, totalcpgB-numoverlapB))
    colnames(x) <- c(labelA, labelB)

    yaxis <- pretty(c(0,ymax))

    b <- barplot(x,
                main="Regions overlap",
                xlab="platform",
                ylab="counts",
                ylim=range(yaxis),
                col=c("salmon","lightblue"),
                axes=F)
    axis(2, at=yaxis,las=2,cex.axis=0.7)
    text(b, x[1,],paste(signif(c(numoverlapA / totalcpgA * 100, numoverlapB / totalcpgB * 100),digits=3), "%"), pos=3)

    mat <- rbind(c(totalcpgA,numoverlapA), c(totalcpgB,numoverlapB))
    colnames(mat) <- c("total", "overlap")
    rownames(mat) <- c(labelA, labelB)
    return(mat)
}


panel.hist <- function(regiongr) {
  x <- width(regiongr[[1]])
  x <- x[x<5000]
  h <- hist(x,
    main=names(regiongr),
    xlab="length (bp)",
    breaks=100,
    col="green")
}


makeOverlapPlot1 <- function(regionsList) {
  oldmfrow <- par()$mfrow
  par(mfrow = c(length(regionsList), length(regionsList)))
  for (i in 1:length(regionsList)) {
    for (j in 1:length(regionsList)) {
      if (i == j){ # diag
        panel.hist(regionsList[i])
      }
      if (i < j) { # upper tri
        panel.overlap(regionsList[i], regionsList[j])
      }
      if (i > j) { # lower tri
        plot(0, type="n", axes=F, xlab="", ylab="")
      }
    }
  }
  par(mfrow=oldmfrow)
}


makeOverlapPlot2 <- function(regionsList) {

  dev.new(width=14, height=14)
  oldmfrow <- par()$mfrow
  par(mfrow = c(length(regionsList)+1, length(regionsList)+1))
  for (i in 1:(length(regionsList)+1)) {
    for (j in 1:(length(regionsList)+1)) {
      if ((j == 1) && (i < 6)) {
          plot(0,0,type="n", axes=F, xlab="",ylab="")
          text(0,0,names(regionsList)[i])
      } else if ((i == 6) && (j > 1)) {
          plot(0,0,type="n", axes=F, xlab="",ylab="")
          text(0,0,names(regionsList)[j-1])
      } else {
        if (j == i+1) { # diag
          # plot(0,0,type="n", axes=F, xlab="",ylab="")
          # text(0,0,"hist")
          panel.hist(regionsList[i])
        }
        if (i < j && j != 6) { # upper tri
          # plot(0,0,type="n", axes=F, xlab="",ylab="")
          # text(0,0,"overlap")
          panel.overlap(regionsList[i],regionsList[j])
        }
        if (i == j || i > j) { # lower tri
          plot(0,0,type="n", axes=F, xlab="",ylab="")
          #text(0,0,"empty")
        }
      }
    }
  }
  par(mfrow=oldmfrow)

}


makeOverlapPlot3 <- function(regionsList) {

  oldmfrow <- par()$mfrow
  par(mfrow = c(length(regionsList)+1, length(regionsList)+1))
  for (i in 1:(length(regionsList)+1)) {
    for (j in 1:(length(regionsList)+1)) {
      if ((j == 1) && (i < 6)) {
          plot(0,0,type="n", axes=F, xlab="",ylab="")
          text(0,0,names(regionsList)[i])
      } else if ((i == 6) && (j > 1)) {
          plot(0,0,type="n", axes=F, xlab="",ylab="")
          text(0,0,names(regionsList)[j-1])
      } else {
        if (j == i+1) { # diag
          # plot(0,0,type="n", axes=F, xlab="",ylab="")
          # text(0,0,"hist")
          panel.hist(regionsList[i])
        }
        if (i < j && j != 6) { # upper tri
          # plot(0,0,type="n", axes=F, xlab="",ylab="")
          # text(0,0,"overlap")
          panel.overlap(regionsList[i],regionsList[j])
        }
        if (i == j || i > j) { # lower tri
          plot(0,0,type="n", axes=F, xlab="",ylab="")
          #text(0,0,"empty")
        }
      }
    }
  }
  par(mfrow=oldmfrow)
}







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



getAnnotCovList <- function(mcgrlist, annotgr, pairmat) {

  cpgcovlist <- list()

  for (i in 1:nrow(pairmat)) {
    j <- pairmat[i,]
    cat(names(mcgrlist)[ j[1] ],"-",names(annotgr)[ j[2] ],"\n")
    x <- numeric( length(annotgr[[ j[2] ]]) )
    if ( length(grep("Agilent",names(mcgrlist)[ j[1] ])) != 0 ) {
      # since agilent is ss we are goint to convert any '+' sites to '-' sites
      z             <- mcgrlist[[ j[1] ]]
      zplus         <- z[strand(z) == "+"]
      zminus        <- z[strand(z) == "-"]
      strand(zplus) <- "-"
      end(zplus)    <- end(zplus)+1
      start(zplus)  <- start(zplus)+1
      z             <- unique(c(zplus,zminus))
      a             <- findOverlaps(annotgr[[ j[2] ]] , z)
    } else {
      a <- findOverlaps(annotgr[[ j[2] ]] , mcgrlist[[ j[1] ]])
    }
    b <- table(queryHits(a))
    x[ as.numeric(names(b)) ] <- as.numeric(b)
    cpgcovlist[[ i ]] <- x
  }

  names(cpgcovlist) <- names(mcgrlist)
  return (cpgcovlist)
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
            axis(1, at = at, labels = label)
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
            axis(2, at = at, labels = label)
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
    #box()
    invisible(list(upper = upper, lower = lower, median = med,
        q1 = q1, q3 = q3))
}


plotAnnotPercCoverageViolin <- function(cpgcovlist, regtotalcpggr, pairmat, what="?what?",h=5) {
  omar = par()$mar
  par(mar=c(7,4,4,2))
  plot(0, type="n",
    main=paste("Percent coverage of CpGs in", what),
    xlab="",
    ylab="percent (%)",
    xlim=c(0, length(cpgcovlist)+1), ylim=c(0,100), axes=F)
  axis(2, las=2)
  axis(1, at=seq(1, length(cpgcovlist)), names(cpgcovlist), las=2)


  for (i in 1:nrow(pairmat)) {
    j <- pairmat[i,]
    a                 <- cpgcovlist[[ j[1] ]]
    totC              <- regtotalcpggr[[ j[2] ]]$numC
    a[ a == 0 ]       <- NA
    totC[ totC == 0 ] <- NA
    b                 <- a / totC * 100
    vioplot2(b[!is.na(b)],at=i,add=T,h=h,col=NA, rectCol="blue", colMed="red", border="gray50")
  }
  par(mar=omar)
}


getRegionCoverage <- function(regionsList, mcgrlist) {
  mat <- cbind(c(3,3,2,2,1,1), c(1,2,3,4,5,6))
  perccovered <- apply(mat, 1, function(i) {
    x <- findOverlaps(regionsList[[ i[1] ]], mcgrlist[[ i[2] ]])
    y <- length(unique(queryHits(x)))
    perccovered <- y / length(regionsList[[ i[1] ]]) * 100
  })
  names(perccovered) <- names(mcgrlist)[ mat[,2] ]
  return(perccovered)
}



plotRegionCoverage <- function(perccovered) {
  x <- c(perccovered[1:2], NA, perccovered[3:4], NA, perccovered[5:6])
  names.arg <- c("A","B","","A","B","","A","B")
  platform  <- c("ERRBS","Agilent","NimbleGen")
  bcols     <- c("lightblue","lightblue","white","salmon","salmon","white","plum2","plum2")

  b <- barplot(x,
    main="Percent of target regions covered",
    xlab="platform",
    ylab="percent covered (%)",
    names.arg=names.arg,
    space=0,
    ylim=c(0,100),
    col=bcols)
  mtext(platform, side=1, at=seq(1,length(b),by=3), line=2) # platform labels
  text(b,x,round(x,digits=1),pos=1)
}


plotOnTarget <- function(props) {
  x <- c(props[1:2],NA,props[3:4],NA,props[5:6])
  names.arg <- c(rep(c("A","B",""),2),"A","B")
  platform  <- c("ERRBS","Agilent","NimbleGen")
  cols      <- c("lightblue","lightblue","white","salmon","salmon","white","plum2","plum2")
  b <- barplot(x,
    ylim=c(0,100),
    main="Percent of target CpGs covered", # was "On-target CpG coverage per respective regions"
    ylab="percent (%)",
    xlab="platform",
    names.arg=names.arg,
    col=cols,
    space=0)
  mtext(platform, side=1, at=seq(1,length(b),by=3), line=2) # platform labels
  text(b,x,round(x,digits=1),pos=1, cex=0.9) # cpg counts in bars
}


plotRegionCoverageSwarm <- function(percRegCovered, ontarget) {
  cols <- c("blue1","blue3","salmon1","salmon3","plum","plum4")
  plot(0,
    type="n",
    main="Percent region coverage",
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

  beeswarm(list(percRegCovered,ontarget),
    pwcol=rep(cols,2),
    pch=21,
    pwbg=rep(cols,2),
    cex=1.5,
    las=1,
    add=T)
  legend("topright",
    names(percRegCovered),
    bg="snow",
    fill=cols,
    border=cols,
    box.col=rgb(0.6,0.6,0.6,alpha=0.5))
}










if(F) {


# Added to the  plot:
# par(mar=c(3.1, 3.1, 1.1, 2.1))
# hist(data,xlim=c(-4,4), col = "pink")
# boxplot(data, horizontal=TRUE,  outline=TRUE,  ylim=c(-4,4), frame=F, col = "green1", add = TRUE)


x <- width(regionsList[[1]])
x <- x[x<5000]
# h <- hist(x,
  # main=names(regionsList)[1],
  # xlab="length (bp)",
  # breaks=100,
  # col="green")
# boxplot(x,horizontal=T)



hist(x, xlim=c(0,5000), col = "pink")
boxplot(x, horizontal=TRUE, frame=F, col = "green1", add=TRUE, outline=T, boxwex=10000, pch="+", outcol=rgb(0.5,0.5,0.5,alpha=0.2))



# --------

load("regionsCpGList.rda")

i=1; j=2
x <- findOverlaps(regionsCpGList[[j]], regionsCpGList[[i]])
length(unique(queryHits(x))) / length(regionsCpGList[[j]])   # 0.6591556
length(unique(subjectHits(x))) / length(regionsCpGList[[i]]) # 0.3695883

i=1; j=3
x <- findOverlaps(regionsCpGList[[j]], regionsCpGList[[i]])
length(unique(queryHits(x))) / length(regionsCpGList[[j]])   # 0.2713194
length(unique(subjectHits(x))) / length(regionsCpGList[[i]]) # 0.3438763


i=2; j=3
x <- findOverlaps(regionsCpGList[[j]], regionsCpGList[[i]])
length(unique(queryHits(x))) / length(regionsCpGList[[j]])   # 0.1562911
length(unique(subjectHits(x))) / length(regionsCpGList[[i]]) # 0.3532851


# --------


library(GenomicRanges)
# library(beeswarm)
# library(vioplot)

load("regionsCpGList.rda")

source("/home/thk2008/bin/methylseqplatformcomparison/scripts/regionAnalysis.R")

ymax <- tail(pretty(c(0, max(sapply(regionsCpGList, length)))),n=1)
for (i in 1:(length(regionsCpGList)-1)) {
  for (j in (i+1):(length(regionsCpGList))) {
    dev.new()
    panel.overlap(regionsCpGList[i], regionsCpGList[j], ymax)
  }
}




# --------

load("mcgrlist.rda")
load("regionsCpGList.rda")

# names(mcgrlist)
# [1] "ERRBS_A"     "ERRBS_B"     "Agilent_A"   "Agilent_B"   "NimbleGen_A"
# [6] "NimbleGen_B" "WGBS_A"      "WGBS_B"
#
# names(regionsCpGList)
# [1] "NimbleGen"   "Agilent"     "MspI_84-334"


i=3; j=2
z             <- mcgrlist[[ i ]]
zplus         <- z[strand(z) == "+"]
zminus        <- z[strand(z) == "-"]
# since region coords are in '+' orientation, convert '-' sites to '+'
start(zminus) <- start(zminus) - 1
end(zminus) <- start(zminus)
z <- unique(c(zplus,zminus))
x <- findOverlaps(z, regionsCpGList[[j]])
length(unique(queryHits(x))) / length(z)                     # 0.7718851
length(unique(subjectHits(x))) / length(regionsCpGList[[j]]) # 0.9188417

i=4; j=2
z             <- mcgrlist[[ i ]]
zplus         <- z[strand(z) == "+"]
zminus        <- z[strand(z) == "-"]
# since region coords are in '+' orientation, convert '-' sites to '+'
start(zminus) <- start(zminus) - 1
end(zminus) <- start(zminus)
z <- unique(c(zplus,zminus))
x <- findOverlaps(z, regionsCpGList[[j]])
length(unique(queryHits(x))) / length(z)                     # 0.805891
length(unique(subjectHits(x))) / length(regionsCpGList[[j]]) # 0.9020747



i=5; j=1
z             <- mcgrlist[[ i ]]
x <- findOverlaps(z, regionsCpGList[[j]])
length(unique(queryHits(x))) / length(z)                     # 0.8335514
length(unique(subjectHits(x))) / length(regionsCpGList[[j]]) # 0.8371862

i=5; j=1
z             <- mcgrlist[[ i ]]
zplus         <- z[strand(z) == "+"]
zminus        <- z[strand(z) == "-"]

# zplus         <- z[strand(z) == "+"]
# x <- findOverlaps(zplus, regionsCpGList[[j]])
# length(unique(queryHits(x)))   # 2356361
# length(unique(subjectHits(x))) # 2356361
#
# zplus         <- z[strand(z) == "+"]
# start(zplus) <- start(zplus) + 1
# end(zplus) <- start(zplus)
# x <- findOverlaps(zplus, regionsCpGList[[j]])
# length(unique(queryHits(x)))   # 2356361
# length(unique(subjectHits(x))) # 2356361
#
# zplus         <- z[strand(z) == "+"]
# start(zplus) <- start(zplus) - 1
# end(zplus) <- start(zplus)
# x <- findOverlaps(zplus, regionsCpGList[[j]])
# length(unique(queryHits(x)))   # 155530
# length(unique(subjectHits(x))) # 155530
#
# zminus         <- z[strand(z) == "-"]
# x <- findOverlaps(zminus, regionsCpGList[[j]])
# length(unique(queryHits(x)))   # 2363333
# length(unique(subjectHits(x))) # 2363333
#
# zminus         <- z[strand(z) == "-"]
# start(zminus) <- start(zminus) + 1
# end(zminus) <- start(zminus)
# x <- findOverlaps(zminus, regionsCpGList[[j]])
# length(unique(queryHits(x)))   # 156079
# length(unique(subjectHits(x))) # 156079
#
# zminus         <- z[strand(z) == "-"]
# start(zminus) <- start(zminus) - 1
# end(zminus) <- start(zminus)
# x <- findOverlaps(zminus, regionsCpGList[[j]])
# length(unique(queryHits(x)))   # 2363333
# length(unique(subjectHits(x))) # 2363333


i=5; j=1
z             <- mcgrlist[[ i ]]
x <- findOverlaps(z, regionsCpGList[[j]])
length(unique(queryHits(x))) / length(z)                     # 0.8335514
length(unique(subjectHits(x))) / length(regionsCpGList[[j]]) # 0.8371862

i=6; j=1
z             <- mcgrlist[[ i ]]
x <- findOverlaps(z, regionsCpGList[[j]])
length(unique(queryHits(x))) / length(z)                     # 0.8282598
length(unique(subjectHits(x))) / length(regionsCpGList[[j]]) # 0.8855895


i=1; j=3
z             <- mcgrlist[[ i ]]
x <- findOverlaps(z, regionsCpGList[[j]])
length(unique(queryHits(x))) / length(z)                     # 0.614592
length(unique(subjectHits(x))) / length(regionsCpGList[[j]]) # 0.5702503

i=2; j=3
z             <- mcgrlist[[ i ]]
x <- findOverlaps(z, regionsCpGList[[j]])
length(unique(queryHits(x))) / length(z)                     # 0.6076746
length(unique(subjectHits(x))) / length(regionsCpGList[[j]]) # 0.5695287





} # end if F





