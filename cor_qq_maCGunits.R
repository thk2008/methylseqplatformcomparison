# correlation

if(FALSE) {
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/cor_qq_maCGunits.R")
}

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




makeSmoothScatter <- function(grA, grB) {

  labelA <- names(grA)
  labelB <- names(grB)
  grA = grA[[1]]
  grB = grB[[1]]
  commonsites <- findOverlaps(grA, grB)
  x <- grA[ queryHits(commonsites)   ]$freqC
  y <- grB[ subjectHits(commonsites) ]$freqC

  corvals <- c(cor(x, y, method="pearson"),
               cor(x, y, method="spearman"))
  
  smoothScatter(x, y,
    xlab=paste(labelA," (",length(grA),")",sep=""),
    ylab=paste(labelB," (",length(grB),")",sep=""),
    transformation = function(x) x^0.70)
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
  x <- grA[ queryHits(commonsites)   ]$freqC
  y <- grB[ subjectHits(commonsites) ]$freqC
  
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
           plot.method="smoothScatter",
           transformation = function(x) x^0.70)
  title(paste(labelA,"::",labelB, " MAplot",sep=""))
  
  # legend("topright",
    # c(paste("IQR:" , signif(IQR(M), digits=3)),
      # paste("MAD:" , signif(mad(M), digits=3)),
      # paste("RMSD:", signif(rmsd(x,y), digits=3))),
    # bty="n", cex=1.2, adj=1)
  legend("topright",
         paste("MAD:" , signif(mad(M), digits=3)),
         bty="n",
         cex=1.5)
    
  # legend("bottomright",
    # c(paste("common sites: ", length(commonsites)),
      # paste(labelA, "total sites: ", length(grA)),
      # paste(labelB, "total sites: ", length(grB))),
    # bty="n", cex=1.2, adj=1)
  legend("bottomright",
         paste("common sites: ", length(commonsites)),
         bty="n",
         cex=1.5)
  
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
      cat(tail(rlabels,n=1), "\n")
      
      grA = mcgrlist[[i]]
      grB = mcgrlist[[j]]
      
      commonsites <- findOverlaps(grA, grB)
      x <- grA[ queryHits(commonsites)   ]$freqC
      y <- grB[ subjectHits(commonsites) ]$freqC
      
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
            main="CpG unit overlap",
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





########
# MAIN #
########


if (!exists("mcgrlist")) {
  cat("loading mcgrlist from mcgrlist_CGunits10x.rda\n")
  load("mcgrlist_CGunits10x.rda")
}

  
if (!exists("statsmat")) {
  rda <- "statsmatCGunits.rda"
  if (!file.exists(rda)){
    cat("computing statsmat\n")
    statsmat <- getStatsMat(mcgrlist)
    cat("saving statsmat as", rda, "\n")
    save(statsmat,file=rda)
    cat("saving statsmat as csv\n")
    saveCsv(statsmat, name="statsmatCGunits")
  } else {
    cat("loading stats mat",rda,"\n")
    load(rda)
  }
}


makepdf <- T
withWGBS <- F

atime  <- format(Sys.time(), "%Y%m%d")
outdir <- paste("cor-qq-ma_CGunits",atime, sep="-")
if (!file.exists(outdir)) {
  dir.create(outdir)
}
cat("output directory is [", outdir, "]\n")

# abtext <- c("A","B","A_opt","B_min","A_opt","B_min","A")

n <- length(mcgrlist)
if (! withWGBS) { n <-  n - 1 }


runMAplots <- function() {
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cat(i,"::",j,"MA-plot\n")
      if (makepdf) { pdf(paste(outdir, paste(names(mcgrlist)[i],names(mcgrlist)[j],"CGunitMAplot.pdf",sep="-"), sep="/")) }
        makeMAplot(mcgrlist[i],mcgrlist[j])
      if (makepdf) { dev.off() }
    }
  }
}



runSmoothScatterPlots <- function() {
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cat(i,"::",j,"smooth scatter cor plot\n")
      if (makepdf) { pdf(paste(outdir, paste(names(mcgrlist)[i],names(mcgrlist)[j],"CGunitSmoothscattercor.pdf",sep="-"), sep="/")) }
        makeSmoothScatter(mcgrlist[i],mcgrlist[j])
      if (makepdf) { dev.off() }
    }
  }
} 
  
runCGunitOverlap <- function() {
  ymax <- tail(pretty(c(0, max(sapply(mcgrlist[1:n], length)))),n=1)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cat(i,"::",j,"cpg overlap barplot\n")
      if (makepdf) { pdf(paste(outdir, paste(names(mcgrlist)[i],names(mcgrlist)[j],"CGunitOverlap.pdf",sep="-"), sep="/")) }
        plotOverlap(mcgrlist[i],mcgrlist[j], ymax)
      if (makepdf) { dev.off() }
    }
  }
}

# convert to png
# for i in *CGunitMAplot.pdf; do j=${i/pdf/png}; convert -density 300 -depth 8 -quality 100 ${i} ${j}; echo "${j} done"; done
# for i in *CGunitOverlap.pdf; do j=${i/pdf/png}; convert -density 300 -depth 8 -quality 100 ${i} ${j}; echo "${j} done"; done
# for i in *CGunitSmoothscattercor.pdf; do j=${i/pdf/png}; convert -density 300 -depth 8 -quality 100 ${i} ${j} &  echo "${j}"; done
