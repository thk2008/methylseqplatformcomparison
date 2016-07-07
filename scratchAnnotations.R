

library(GenomicRanges)
library(methylKit)
library("BSgenome.Hsapiens.UCSC.hg19")
require(vioplot)

source("/home/thk2008/bin/methylseqplatformcomparison/scripts/annotations.R")

load("mcgrlist.rda")

if(F) {
load("gene.obj.rda")
x <- gene.obj$exons
y <- unique(x)
gene.obj$exons <- y
x <- gene.obj$introns
y <- unique(x)
gene.obj$introns <- y
x <- gene.obj$promoters
y <- unique(x)
gene.obj$promoters <- y
rm(x,y) 


exontotalcpg  <- getRegionsTotalCpGs(gene.obj$exons)
exoncpgcovmat <- getAnnotCovMat(mcgrlist, exontotalcpg)
dev.new()
plotAnnotPercCoverageViolin(exoncpgcovmat, exontotalcpg$numC , what="Exons")


promotertotalcpg  <- getRegionsTotalCpGs(gene.obj$promoters)
promotercpgcovmat <- getAnnotCovMat(mcgrlist, promotertotalcpg)
dev.new()
plotAnnotPercCoverageViolin(promotercpgcovmat, promotertotalcpg$numC , what="Promoters")
  

introntotalcpg  <- getRegionsTotalCpGs(gene.obj$introns)
introncpgcovmat <- getAnnotCovMat(mcgrlist, introntotalcpg)
dev.new()
plotAnnotPercCoverageViolin(introncpgcovmat, introntotalcpg$numC , what="Introns")
}

load("cpg.obj.rda")
x <- cpg.obj$CpGi
y <- unique(x)
cpg.obj$CpGi <- y
x <- cpg.obj$shores
y <- unique(x)
cpg.obj$shores <- y




cpgitotalcpg  <- getRegionsTotalCpGs2(cpg.obj$CpGi, stranded=FALSE)
cpgicpgcovmat <- getAnnotCovMat(mcgrlist, cpgitotalcpg)
dev.new()
plotAnnotPercCoverageViolin(cpgicpgcovmat, cpgitotalcpg$numC, what="CpGi")


shorestotalcpg  <- getRegionsTotalCpGs2(cpg.obj$shores,stranded=FALSE)
shorescpgcovmat <- getAnnotCovMat(mcgrlist, shorestotalcpg)
dev.new()
plotAnnotPercCoverageViolin(shorescpgcovmat, shorestotalcpg$numC , what="Shores")




  #apply(exoncpgcovmat, 2, function(i) { length(which(i > 0)) }) / length(exontotalcpg) * 100  # identical to percExon




