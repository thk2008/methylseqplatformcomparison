

if (F) {
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/meanCGcovBox.R")
}

require(GenomicRanges)

plotMeanCpgCovBox <- function(mcgrlist,main="") {
  cgcov <- list()
  for (i in 1:length(mcgrlist)) {
    cgcov[[i]] <- mcgrlist[[i]]$coverage
  }
  meancov <- sapply(cgcov, mean)
  names(cgcov) <- names(mcgrlist)
  bcols <- c("lightblue","lightblue","salmon","salmon","plum2","plum2","lightgreen")
  omar <- par()$mar
  par(mar=c(9.1, 4.1, 4.1, 2.1))
  b <- boxplot(cgcov,
   main=main,
   xlab="",
   ylab="",
   outline=F,
   col=bcols,
   las=2)
  mtext("outliers removed",adj=1)
  points(seq(1:length(meancov)), meancov, pch=3)
  legend("topright",legend="mean",bty="n",pch=3)
  par(mar=omar)
}



run  <- function() {

  makepdf <- T
  
  if (!exists("mcgrlist")) {
    rda <- "mcgrlist.rda"
    if (!file.exists(rda)){
      stop("cannot find mcgrlist.rda [",rda,"]")
    } else {
      cat("loading cg sites [",rda,"]\n")
      load(rda)
    }
  }
  
  
  # if (makepdf) { pdf(paste(outdir, "cgSiteCoverageBox.pdf", sep="/"),pointsize=11) }
  if (makepdf) { pdf("cgSiteCoverageBox.pdf", pointsize=11) }
   cat("plotting box plot of cg sites \n")
   plotMeanCpgCovBox(mcgrlist,main="CpG site coverage >=10x")
  if (makepdf) { dev.off() }
  
  
  
  rda <- "mcgrlist_CGunits10x.rda"
  if (file.exists(rda)) {
    cat("found cg units file [",rda,"]\n")
    load(rda)
  }
  
  if (makepdf) { pdf("cgUnitCoverageBox.pdf", pointsize=11) }
   cat("plotting box plot of cg units \n")
   plotMeanCpgCovBox(mcgrlist,main="CG unit coverage >=10x")
  if (makepdf) { dev.off() }

}
