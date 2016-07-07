# correlation

require(GenomicRanges)
require(affy)
require(preprocessCore)


panelCorPearson <- function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method="pearson"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panelCorSpearman <- function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method="spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panelSmooth<-function(x, y, col = par("col"), bg = NA, pch = par("pch"), cex = 1, col.smooth = "darkgreen", span = 2/3, iter = 3, ...) {
  par(new = TRUE) #par(usr = c(usr[1:2], 0, 1.5) )
  smoothScatter(x, y,colramp=colorRampPalette(topo.colors(100)), bg = bg)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
  abline(lm(y[ok]~x[ok]), col="red")
  lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = col.smooth, ...)
  qqq <- qqplot(x,y,pch=".",plot.it=F)
  points(qqq,pch=".")
}

panelHist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h      <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB     <- length(breaks)
  y      <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}



if(F) {
  makeMAplotPairs(mcgrlist)
}

makeMAplotPairs <- function(mcgrlist) {
  
  stop("this doesn't work right. graphs are squished. margins are too large.")
  # i suppose we could make a gigantic device to plot it on

  # x <- (length(mcgrlist)^2 - length(mcgrlist) ) / 2
  
  oldmfrow <- par()$mfrow
  #par(mfrow = c(length(mcgrlist), length(mcgrlist)))
  par(mfrow = c(7,3))
  for (i in 1:(length(mcgrlist)-1)) {
    for (j in (i+1):length(mcgrlist)) {
      makeMAplot(mcgrlist[i],mcgrlist[j])
    }
  }
  
  # for (i in 1:(length(mcgrlist)-1)) {
    # for (j in (i+1):length(mcgrlist)) {
      # labelA <- names(mcgrlist)[i]
      # labelB <- names(mcgrlist)[j]
      # grA = mcgrlist[[i]]
      # grB = mcgrlist[[j]]
      # commonsites <- findOverlaps(grA, grB)
      # x <- grA[ queryHits(commonsites)   ]$numC
      # y <- grB[ subjectHits(commonsites) ]$numC
      # x[x == 0] <- 0.01
      # y[y == 0] <- 0.01
      # zz <- cbind(x,y)
      # ma.plot( rowMeans(log2(zz)),           # A
               # log2(zz[, 1])-log2(zz[, 2]),  # M
               # cex=1,
               # add.loess=F,
               # plot.method="smoothScatter")
      # title(paste(labelA,"::",labelB, "  MAplot",sep=""))
    # }
  # }
  
  par(mfrow=oldmfrow)
}



makeSmoothScatter <- function(grA, grB) {

  labelA <- names(grA)
  labelB <- names(grB)
  grA = grA[[1]]
  grB = grB[[1]]
  commonsites <- findOverlaps(grA, grB)
  x <- grA[ queryHits(commonsites)   ]$numC
  y <- grB[ subjectHits(commonsites) ]$numC

  corvals <- c(cor(x, y, method="pearson"),
               cor(x, y, method="spearman"))
  
  smoothScatter(x, y,
    xlab=paste(labelA," (",length(grA),")",sep=""),
    ylab=paste(labelB," (",length(grB),")",sep="")
  )
  abline(lm(y ~ x), col="red")
  lines(stats::lowess(x, y, f=2/3, iter=3), col="darkgreen")
  points(qqplot(x, y, pch=".", plot.it=F), pch=".")
  
  title(paste("Methylation level correlation ", labelA,"::",labelB, sep=""))
  mtext( paste(
    "common sites: ",     length(commonsites),
    ", cor (pearson): ",  round(corvals[1],digits=2),
    ", cor (spearman): ", round(corvals[2],digits=2),
    sep="")
  )

}


# rmse <- function(m, o) { sqrt( mean( (m - o)^2 ) ) }
rmsd <- function(x, y) { sqrt( mean( (x - y)^2 ) ) }


makeMAplot <- function(grA, grB) {

  labelA <- names(grA)
  labelB <- names(grB)
  grA = grA[[1]]
  grB = grB[[1]]
  commonsites <- findOverlaps(grA, grB)
  x <- grA[ queryHits(commonsites)   ]$numC
  y <- grB[ subjectHits(commonsites) ]$numC
  
  # corvals <- c(cor(x, y, method="pearson"),
               # cor(x, y, method="spearman"))

  x[x == 0] <- 0.01
  y[y == 0] <- 0.01
  zz <- cbind(x,y)
  M  <- log2(zz[, 1]) - log2(zz[, 2])

  ma.plot( rowMeans(log2(zz)), # A average methylation level
           M,                  # M log ratio of methylation level
           cex=1,
           add.loess=F,
           show.statistics=F,
           plot.method="smoothScatter")
  title(paste(labelA,"::",labelB, "  MAplot",sep=""))
  legend(6,14.8,
    c(paste("IQR:" , signif(IQR(M), digits=3)),
      paste("MAD:" , signif(mad(M), digits=3)),
      paste("RMSD:", signif(rmsd(x,y), digits=3))),
    bty="n", cex=1.2, adj=1)
  legend(6,-9.2,
    c(paste("common sites: ", length(commonsites)),
      paste(labelA, "total sites: ", length(grA)),
      paste(labelB, "total sites: ", length(grB))),
    bty="n", cex=1.2, adj=1)
  
  # ?normalize.quantiles
   # This method is based upon the concept of a quantile-quantile plot
   # extended to n dimensions. No special allowances are made for
   # outliers. If you make use of quantile normalization please cite
   # Bolstad et al, Bioinformatics (2003).

  # zzz <- normalize.quantiles(zz)
  # #dev.new()
  # ma.plot( rowMeans(log2(zzz)), log2(zzz[, 1])-log2(zzz[, 2]), cex=1, plot.method="smoothScatter" )
  # title(paste("Post Norm",labelA,labelB))
  # cat("post norm ma done\n")

}


getStatsMat <- function(mcgrlist) {
  statsmat <- numeric()
  rlabels  <- character()
  for (i in 1:(length(mcgrlist)-1)) {
    for (j in (i+1):length(mcgrlist)) {
      rlabels <- c(rlabels, paste(names(mcgrlist)[i], "::", names(mcgrlist)[j], sep=""))
      cat(rlabels, "\n")
      
      grA = mcgrlist[[i]]
      grB = mcgrlist[[j]]
      
      commonsites <- findOverlaps(grA, grB)
      x <- grA[ queryHits(commonsites)   ]$numC
      y <- grB[ subjectHits(commonsites) ]$numC
      
      pcor <- cor(x, y, method="pearson")
      scor <- cor(x, y, method="spearman")
    
      x[x == 0] <- 0.01
      y[y == 0] <- 0.01
      zz <- cbind(x,y)
      M <- log2(zz[, 1]) - log2(zz[, 2])        
        
     statsmat <- rbind(
       statsmat,
       c(
         length(grA),                 # total sites A
         length(grB),                 # total sites B
         length(commonsites),         # number of common sites (overlap)
         pcor,                        # pearson correlation
         scor,                        # spearman correlation
         signif(IQR(M), digits=3),    # IQR
         signif(mad(M), digits=3),    # MAD
         signif(rmsd(x,y), digits=3))) # RMSD
    }
  }
  colnames(statsmat) <- c("totalA","totalB","common","pearson","spearman","iqr","mad","rmsd")
  rownames(statsmat) <- rlabels
  
  return(statsmat)
}


saveCsv <- function(mat, name="outmat", outdir=".") {
  stopifnot(class(mat) == "matrix")
  afile <- paste(outdir, paste(name, "csv", sep="."), sep="/")
  cat(",",file=afile)
  ow <- options("warn")
  options(warn = -1)
  write.table(mat,file=afile,append=T,quote=F,sep=",")
  options(ow) # reset
}


# raplot

# grA <- mcgrlist[1]
# grB <- mcgrlist[2]
# labelA <- names(grA)
# labelB <- names(grB)
# grA = grA[[1]]
# grB = grB[[1]]
# commonsites <- findOverlaps(grA, grB)
# a <- grA[ queryHits(commonsites)   ]
# b <- grB[ subjectHits(commonsites) ]
# a <- round((a$coverage * a$numC) / 100)
# b <- round((b$coverage * b$numC) / 100)
# library(caroline)
# raPlot(a[1:1000000],b[1:1000000])

# cc <- cbind(a,b)
# M <- log2(cc[, 1]) - log2(cc[, 2])
# A <- rowMeans(log2(cc))




plotOverlap <- function(grA, grB, ymax) {
    labelA <- names(grA)
    labelB <- names(grB)
    i <- grA[[1]]
    j <- grB[[1]]
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
            main="CpG overlap",
            xlab="",
            ylab="counts",
            ylim=range(yaxis),
            col=c("salmon","lightblue"),
            axes=F)
    axis(2, at=yaxis, las=2, cex.axis=0.7)    
    text(b, x[1,],paste(signif(c(numoverlapA / totalcpgA * 100, numoverlapB / totalcpgB * 100),digits=3), "%"), pos=3)
    
    mat <- rbind(c(totalcpgA,numoverlapA), c(totalcpgB,numoverlapB))
    colnames(mat) <- c("total", "overlap")
    rownames(mat) <- c(labelA, labelB)
    return(mat)
}




# what are the characteristics of non-overlapping sites?
# where are they covering?





# i="corvsma.Agilent_A.Agilent_B.pdf"
# j=${i%.pdf} && j=${j#corvsma.}

# convert -density 300 -depth 8 -quality 100 ${i}[0] cor.${j}.png
# convert -density 300 -depth 8 -quality 100 ${i}[1] prenormma.${j}.png
# convert -density 300 -depth 8 -quality 100 ${i}[2] postnormma.${j}.png



# cbind(statsmat[,3] / statsmat[,1], statsmat[,3] / statsmat[,2])

# -----------------------------------------------
# Agilent_A::Agilent_B     0.91580655 0.9783680
# -----------------------------------------------
# ERRBS_A::Agilent_A       0.19853694 0.3414184
# ERRBS_A::Agilent_B       0.18869135 0.3466539
# ERRBS_B::Agilent_A       0.20580000 0.3574843
# ERRBS_B::Agilent_B       0.19553398 0.3628543

# Agilent_A::WGBS_A        0.11148417 0.1253167
# Agilent_A::WGBS_B        0.09979163 0.1300546
# Agilent_B::WGBS_A        0.11319698 0.1191056
# Agilent_B::WGBS_B        0.10147218 0.1237884

# -----------------------------------------------

# Agilent_A::NimbleGen_A   0.59534970 0.4053562
# Agilent_A::NimbleGen_B   0.63589838 0.4067019
# Agilent_B::NimbleGen_A   0.61329092 0.3908704
# Agilent_B::NimbleGen_B   0.65387314 0.3914565

# -----------------------------------------------
# NimbleGen_A::NimbleGen_B 0.94257200 0.8853975
# -----------------------------------------------
# NimbleGen_A::WGBS_A      0.11570868 0.1910280
# NimbleGen_A::WGBS_B      0.10311682 0.1973769
# NimbleGen_B::WGBS_A      0.11423817 0.2007791
# NimbleGen_B::WGBS_B      0.10209606 0.2080425

# ERRBS_B::NimbleGen_A     0.29732320 0.3516455
# ERRBS_B::NimbleGen_B     0.32720933 0.3635178
# ERRBS_A::NimbleGen_A     0.28757277 0.3367116
# ERRBS_A::NimbleGen_B     0.31570963 0.3472337

# -----------------------------------------------
# ERRBS_A::ERRBS_B         0.89829762 0.8893123
# -----------------------------------------------
# ERRBS_A::WGBS_A          0.09065643 0.1752428
# ERRBS_A::WGBS_B          0.08115504 0.1818832
# ERRBS_B::WGBS_A          0.08947770 0.1747118
# ERRBS_B::WGBS_B          0.08050650 0.1822527
# -----------------------------------------------
# WGBS_A::WGBS_B           0.44844990 0.5199357





