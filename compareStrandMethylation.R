
if (FALSE) {
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/compareStrandMethylation.R")
}


require(GenomicRanges)
require(affy)
require(preprocessCore)


rmsd <- function(x, y) { sqrt( mean( (x - y)^2 ) ) }

makepdf <- F

load("errbsAmcgr.rda")


# watson = top = plus
# crick = bottom = minus

# mcplus <- errbsA[strand(errbsA) == '+']
# mcminus <- errbsA[strand(errbsA) == '-']

# remove methylation levels between 25-75%
# do we have to make a cutoff somewhere?
# mcgr <- errbsA[-which(errbsA$numC > 25 & errbsA$numC < 75)]

# mcgr <- mcgr[-which(mcgr$numC == 0)]

mcplus <- mcgr[strand(mcgr) == '+']
mcminus <- mcgr[strand(mcgr) == '-']

# convert minus into plus
start(mcminus)  <- start(mcminus)-1
end(mcminus)    <- end(mcminus)-1
strand(mcminus) <- "+"

overlaps <- findOverlaps(mcminus,mcplus)

x <- mcplus[subjectHits(overlaps)]
y <- mcminus[queryHits(overlaps)]
xx <- x$numC
yy <- y$numC

smoothScatter(xx, yy, colramp=colorRampPalette(topo.colors(100)),
xlab="Watson strand methylation level",ylab="Crick strand methylation level")
title("ERRBS_A +/- strand methylation levels")
abline(lm(yy~xx), col="red")
lines(stats::lowess(xx, yy, f = 2/3, iter = 3),col="darkgreen")
qqq <- qqplot(xx, yy, plot.it=F)
points(qqq,pch=".")







dev.new()

xx[xx == 0] <- 0.01
yy[yy == 0] <- 0.01
zz <- cbind(xx,yy)
M <- log2(zz[, 1]) - log2(zz[, 2]) # M log ratio of methylation level
A <- rowMeans(log2(zz))            # A average methylation level
ma.plot( A,
         M,
         cex=1,
         add.loess=F,
         show.statistics=F,
         plot.method="smoothScatter")
title("all of em")
legend("topright",
  c(paste("IQR:" , signif(IQR(M), digits=3)),
    paste("MAD:" , signif(mad(M), digits=3)),
    paste("RMSD:", signif(rmsd(xx,yy), digits=3))),
  bty="n", cex=1.2, adj=1)
  
  
x <- mcplus[subjectHits(overlaps)]
y <- mcminus[queryHits(overlaps)]
xx <- x$numC
yy <- y$numC
zz <- abs(xx - yy)
xxx <- xx[-which(zz == 100)]
yyy <- yy[-which(zz == 100)]
zzz <- zz[-which(zz == 100)]

dev.new()
smoothScatter(xxx, yyy, colramp=colorRampPalette(topo.colors(100)))
title("difference of 100 removed")
abline(lm(yyy~xxx), col="red")
lines(stats::lowess(xxx, yyy, f = 2/3, iter = 3),col="darkgreen")
qqq <- qqplot(xxx, yyy, plot.it=F)
points(qqq,pch=".")

dev.new()

xxx[xxx == 0] <- 0.01
yyy[yyy == 0] <- 0.01
zz <- cbind(xxx,yyy)
M <- log2(zz[, 1]) - log2(zz[, 2]) # M log ratio of methylation level
A <- rowMeans(log2(zz))            # A average methylation level
ma.plot( A,
         M,
         cex=1,
         add.loess=F,
         show.statistics=F,
         plot.method="smoothScatter")
title("difference of 100 removed")
legend("topright",
  c(paste("IQR:" , signif(IQR(M), digits=3)),
    paste("MAD:" , signif(mad(M), digits=3)),
    paste("RMSD:", signif(rmsd(xxx,yyy), digits=3))),
  bty="n", cex=1.2, adj=1)


yyy <- sample(0:100, length(yyy), replace = TRUE)

dev.new()
smoothScatter(xxx, yyy, colramp=colorRampPalette(topo.colors(100)))
title("difference of 100 removed")
abline(lm(yyy~xxx), col="red")
lines(stats::lowess(xxx, yyy, f = 2/3, iter = 3),col="darkgreen")
qqq <- qqplot(xxx, yyy, plot.it=F)
points(qqq,pch=".")

dev.new()

xxx[xxx == 0] <- 0.01
yyy[yyy == 0] <- 0.01
zz <- cbind(xxx,yyy)
M <- log2(zz[, 1]) - log2(zz[, 2]) # M log ratio of methylation level
A <- rowMeans(log2(zz))            # A average methylation level
ma.plot( A,
         M,
         cex=1,
         add.loess=F,
         show.statistics=F,
         plot.method="smoothScatter")
title("difference of 100 removed")
legend("topright",
  c(paste("IQR:" , signif(IQR(M), digits=3)),
    paste("MAD:" , signif(mad(M), digits=3)),
    paste("RMSD:", signif(rmsd(xxx,yyy), digits=3))),
  bty="n", cex=1.2, adj=1)
  
  
  
