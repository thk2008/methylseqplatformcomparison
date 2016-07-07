

if(FALSE){
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/strandSymmetry.R")
}


require(GenomicRanges)
require(affy)
require(preprocessCore)









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
           plot.method="smoothScatter")
  title(paste(labelA,"::",labelB, "  MAplot",sep=""))
  legend("topright",
    c(paste("IQR:" , signif(IQR(M), digits=3)),
      paste("MAD:" , signif(mad(M), digits=3)),
      paste("RMSD:", signif(rmsd(x,y), digits=3))),
    bty="n", cex=1.2, adj=1)
  legend("bottomright",
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




run <- function() {

  if (!exists("mcgrlist")) {
    load("mcgrlistCpGsites10x.rda")
  }
  print(names(mcgrlist))
  
  makepdf <- T
  
  atime  <- format(Sys.time(), "%Y%m%d")
  outdir <- paste("strandSymmetry",atime, sep="-")
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  cat("output directory is [", outdir, "]\n")
  
  for (i in 1:length(mcgrlist)) {
    sname    <- names(mcgrlist)[i]
    cat(sname, ", ")
    x        <- mcgrlist[[i]]
    totalsites <- length(x)
    xplus  <- x[strand(x) == "+"]
    xminus <- x[strand(x) == "-"]  
    start(xminus)  <- start(xminus) - 1
    end(xminus)    <- end(xminus) - 1
    strand(xminus) <- "+"
    overlaps <- findOverlaps(xplus,xminus)
    x <- xplus[queryHits(overlaps)]
    y <- xminus[subjectHits(overlaps)]

    x <- x$numC / x$coverage * 100
    y <- y$numC / y$coverage * 100
    
    labelA <- "plus"
    labelB <- "minus"
  
   if(F) {
    corvals <- c(cor(x, y, method="pearson"),
                 cor(x, y, method="spearman"))

    cat("smoothscatter , ")
    if (makepdf) { pdf(paste(outdir, paste(sname,"strandSymmetrySmoothscatter.pdf",sep="-"), sep="/")) }
                 
    smoothScatter(x, y,
      xlab=labelA,
      ylab=labelB)
    abline(lm(y ~ x), col="red")
    lines(stats::lowess(x, y, f=2/3, iter=3), col="darkgreen")
    points(qqplot(x, y, pch=".", plot.it=F), pch=".")
    
    title(paste( sname, " strand symmetry methylation level correlation", sep=""))
    mtext( paste(
      "total sites: ",           totalsites,
      ", complementary sites: ", length(overlaps), "(",round(length(overlaps)/length(mcgrlist[[i]])*100,digits=2),"%)",
      ", r(p)= ",      round(corvals[1],digits=2),
      ", r(s)= ",      round(corvals[2],digits=2),
      sep=""),cex=0.9)

    if (makepdf) { dev.off() }
   }
   
    x[x == 0] <- 0.01
    y[y == 0] <- 0.01
    zz <- cbind(x,y)
    M  <- log2(zz[, 1]) - log2(zz[, 2])
    A  <- rowMeans(log2(zz))
    
    cat("ma-plot\n")
    if (makepdf) { pdf(paste(outdir, paste(sname,"strandSymmetryMAplot.pdf",sep="-"), sep="/")) }
    
    ma.plot( A,  # A average methylation level
             M,  # M log ratio of methylation level
             add.loess=F,
             show.statistics=F,
             plot.method="smoothScatter",
             transformation = function(x) x^0.70)
    title(paste( sname, " strand symmetry methylation level concordance", sep=""), cex.main=1)
    # legend("topright", c(paste("IQR:" , signif(IQR(M), digits=3)), paste("MAD:" , signif(mad(M), digits=3)), paste("RMSD:", signif(rmsd(x,y), digits=3))), bty="n", cex=1.2, adj=1)
    legend("topright", paste("MAD:" , signif(mad(M), digits=3)), bty="n", cex=1.5)
    # legend("bottomright", c(paste("complementary sites: ", length(overlaps)), paste("total sites: ", totalsites)), bty="n", cex=1.2, adj=1)      
    legend("bottomright", c(paste("complementary sites: ", length(overlaps), " (",round(length(overlaps) / totalsites * 100, digits=1),"%)", sep="")), bty="n", cex=1.5)      
    
    if (makepdf) { dev.off() }
      
  }

  print("done.")
}

# for (i in 1:(length(mcgrlist)-1)) {
  # for (j in (i+1):length(mcgrlist)) {
    # cat(i,"::",j,"MA-plot\n")
    # if (makepdf) { pdf(paste(outdir, paste(names(mcgrlist)[i],names(mcgrlist)[j],"CGunitMAplot.pdf",sep="-"), sep="/")) }
      # makeMAplot(mcgrlist[i],mcgrlist[j])
    # if (makepdf) { dev.off() }
  # }
# }
# 
# for (i in 1:(length(mcgrlist)-1)) {
  # for (j in (i+1):length(mcgrlist)) {
    # cat(i,"::",j,"smooth scatter cor plot\n")
    # if (makepdf) { pdf(paste(outdir, paste(names(mcgrlist)[i],names(mcgrlist)[j],"CGunitSmoothscattercor.pdf",sep="-"), sep="/")) }
      # makeSmoothScatter(mcgrlist[i],mcgrlist[j])
    # if (makepdf) { dev.off() }
  # }
# }


# what are the clouds in the MAplots?
# not a factor of coverage, but where the methylation on one strand is 0
# and the other strand is 0 < x < 50
# small percentage of all sites ~10% (this one was 8%)

if (F) {

i=1

sname    <- names(mcgrlist)[i]
cat(sname, ", ")
x        <- mcgrlist[[i]]
totalsites <- length(x)
xplus  <- x[strand(x) == "+"]
xminus <- x[strand(x) == "-"]  
start(xminus)  <- start(xminus) - 1
end(xminus)    <- end(xminus) - 1
strand(xminus) <- "+"
overlaps <- findOverlaps(xplus,xminus)
x <- xplus[queryHits(overlaps)]
y <- xminus[subjectHits(overlaps)]

x <- x$numC
y <- y$numC

labelA <- "plus"
labelB <- "minus"

corvals <- c(cor(x, y, method="pearson"),
             cor(x, y, method="spearman"))


x[x == 0] <- 0.01
y[y == 0] <- 0.01
zz <- cbind(x,y)

smoothScatter(x, y, xlab=labelA, ylab=labelB)   

M  <- log2(zz[, 1]) - log2(zz[, 2])
A  <- rowMeans(log2(zz)) 
ma.plot(A, M, cex=1, add.loess=F, show.statistics=F, plot.method="smoothScatter")

####################################################################################
####################################################################################
####################################################################################

x <- xplus[queryHits(overlaps)]
y <- xminus[subjectHits(overlaps)]

x$numC[x$numC == 0] <- 0.01
y$numC[y$numC == 0] <- 0.01

M  <- log2(x$numC) - log2(y$numC)
A  <- rowMeans(log2(cbind(x$numC,y$numC)))
ma.plot(A, M, cex=1, add.loess=F, show.statistics=F, plot.method="smoothScatter")


M2 <- M[which(A > -5 & A < -0.5)]
A2 <- A[which(A > -5 & A < -0.5)]
ma.plot(A2, M2, cex=1, add.loess=F, show.statistics=F, plot.method="smoothScatter")

##################################

x1 <- x[which(A > -5 & A < -0.5)]
y1 <- y[which(A > -5 & A < -0.5)]
z1 <- cbind(x1$numC,y1$numC)

x2 <- x1$numC
y2 <- y1$numC
M3 <- log2(x2) - log2(y2)
A3 <- rowMeans(log2(cbind(x2,y2)))
ma.plot(A3, M3, cex=1, add.loess=F, show.statistics=F, plot.method="smoothScatter")

} # end if(F)

