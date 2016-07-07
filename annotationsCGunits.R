
if (FALSE) {
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/annotationsCGunits.R")
}

require(GenomicRanges)
require(methylKit)
require("BSgenome.Hsapiens.UCSC.hg19")
require(vioplot)
require(beeswarm)

getAnnotRegionsTotalCpGs <- function(regiongr) {
# require("BSgenome.Hsapiens.UCSC.hg19")
  gr         <- trim(regiongr)
  # sometimes there is a C on the end of a region that may be a CG
  # which we see in the methylcalls but not in the designed regions
  tmp        <- gr
  start(tmp) <- end(tmp)
  x          <- getSeq(Hsapiens, tmp)
  y          <- vcountPattern("C", x)
  end(tmp[which(y == 1)]) <- end(tmp[which(y == 1)]) + 1
  end(gr)    <- end(tmp)
  regionSeqs <- getSeq(Hsapiens, gr)
  numCGss    <- vcountPattern("CG", regionSeqs)
  return(numCGss)
}

if(F){

cat("reading gene obj from refseq.hg19.bed.txt\n")
gene.obj <- read.transcript.features(paste("data", "refseq.hg19.bed.txt", sep="/"))

x <- gene.obj$exons
y <- unique(x)
gene.obj$exons <- y

x <- gene.obj$introns
y <- unique(x)
gene.obj$introns <- y

x <- gene.obj$promoters
y <- unique(x)
gene.obj$promoters <- y

cat("saving gene obj to gene.obj.rda\n")
save(gene.obj, file="gene.obj.rda")

cat("reading cpg obj from cpgis.hg19.bed.txt\n")
cpg.obj  <- read.feature.flank(paste("data", "cpgis.hg19.bed.txt",sep="/"), feature.flank.name=c("CpGi","shores"))

x <- cpg.obj$CpGi
y <- unique(x)
cpg.obj$CpGi <- y

x <- cpg.obj$shores
y <- unique(x)
cpg.obj$shores <- y

cat("saving cpg obj to cpg.obj.rda\n")
save(cpg.obj, file="cpg.obj.rda")


numCGss           <- list()
cat("computing num CGs in exons\n")
numCGss$exons     <- getAnnotRegionsTotalCpGs(gene.obj$exons)
cat("computing num CGs in introns\n")
numCGss$introns   <- getAnnotRegionsTotalCpGs(gene.obj$introns)
cat("computing num CGs in promoters\n")
numCGss$promoters <- getAnnotRegionsTotalCpGs(gene.obj$promoters)
cat("computing num CGs in cpg islands\n")
numCGss$CpGi      <- getAnnotRegionsTotalCpGs(cpg.obj$CpGi)
cat("computing num CGs in shores\n")
numCGss$shores    <- getAnnotRegionsTotalCpGs(cpg.obj$shores)

cat("saving numCGss to numCGssAnnotationRegions.rda\n")
save(numCGss, file="numCGssAnnotationRegions.rda")


cat("loading mcgrlist from mcgrlist_CGunits10x.rda\n")
load("mcgrlist_CGunits10x.rda")

} # end if(F)



for (i in 1:length(mcgrlist)) {

  cat("computing overlap",names(mcgrlist)[i],":  ")
  mcgrlist[[i]]$exons     <- 0
  mcgrlist[[i]]$introns   <- 0
  mcgrlist[[i]]$promoters <- 0
  mcgrlist[[i]]$cpgi      <- 0
  mcgrlist[[i]]$shores    <- 0

  cat("exons, ")
  x <- findOverlaps(mcgrlist[[i]], gene.obj$exons, ignore.strand=TRUE)
  mcgrlist[[i]][ queryHits(x) ]$exons <- 1

  cat("introns, ")
  x <- findOverlaps(mcgrlist[[i]], gene.obj$introns, ignore.strand=TRUE)
  mcgrlist[[i]][ queryHits(x) ]$introns <- 1

  cat("promoters, ")
  x <- findOverlaps(mcgrlist[[i]], gene.obj$promoters, ignore.strand=TRUE)
  mcgrlist[[i]][ queryHits(x) ]$promoters <- 1

  cat("cpgi, ")
  x <- findOverlaps(mcgrlist[[i]], cpg.obj$CpGi, ignore.strand=TRUE)
  mcgrlist[[i]][ queryHits(x) ]$cpgi <- 1

  cat("shores")
  x <- findOverlaps(mcgrlist[[i]], cpg.obj$shores, ignore.strand=TRUE)
  mcgrlist[[i]][ queryHits(x) ]$shores <- 1
  cat("\n")

  save(mcgrlist, file="mcgrlistCGunitsAnnotations.rda")
}


if(F) {
  annotCountsMat <- sapply(mcgrlist, function(i) { sapply(elementMetadata(i)[6:10], sum) })
  num.annots     <- lapply(mcgrlist, function(i) { apply(as.matrix(values(i)[6:10]), 1, sum) })
  unannotated    <- sapply(num.annots, function(i) { length(which(i == 0)) })
  annotCountsMat <- rbind(annotCountsMat, unannotated)
  totalCG        <- sapply(mcgrlist, length)
  percAnnotMat   <- t(annotCountsMat) / totalCG * 100

  save(annotCountsMat, num.annots, percAnnotMat, file="numberAnnotated.rda")
}




load("mcgrlistCGunitsAnnotations.rda")
load("cpg.obj.rda")
load("gene.obj.rda")
load("numCGssAnnotationRegions.rda")
load("numberAnnotated.rda")
load("genomicAndPercRegionCoverage.rda")

makepdf <- T

atime  <- format(Sys.time(), "%Y%m%d")
outdir <- paste("annotationsCGunits",atime, sep="-")
if (!file.exists(outdir)) {
  dir.create(outdir)
}
cat("output directory is [", outdir, "]\n")


x <- list()
for (i in 1:ncol(percAnnotMat)) {
  x[[i]] <- percAnnotMat[1:6,i] # no WGBS
}
names(x) <- colnames(percAnnotMat)

bscols <- c("blue1","blue3","salmon1","salmon3","plum","plum4")

if (makepdf) { pdf(paste(outdir, "annotationsSampleCGunit.pdf", sep="/")) }


#b <- beeswarm(x,do.plot=F)
#plot(b$x,b$y,type="n")

beeswarm(x,
  labels=names(x),
  main="Sample CG units annotations",
  xlab="genomic feature",
  ylab="percent (%)",
  ylim=c(5,50),
  pwcol=rep(bscols,6),
  pch=21,
  pwbg=rep(bscols,6),
  cex=1.5)

legend("topright", legend=rownames(percAnnotMat)[1:6], fill=bscols, border=bscols, box.col=rgb(0.6,0.6,0.6,alpha=0.5))

if (makepdf) { dev.off() }
if (makepdf) { pdf(paste(outdir, "annotationsSampleCGunit0-100.pdf", sep="/")) }

beeswarm(x,
  labels=names(x),
  main="Sample CG units annotations",
  xlab="genomic feature",
  ylab="percent (%)",
  ylim=c(0,100),
  pwcol=rep(bscols,6),
  pch=21,
  pwbg=rep(bscols,6),
  cex=1.5)
legend("topright", legend=rownames(percAnnotMat)[1:6], fill=bscols, border=bscols, box.col=rgb(0.6,0.6,0.6,alpha=0.5))

if (makepdf) { dev.off() }


# region CG coverage

#load("genomicAndPercRegionCoverage.rda")

genomicCoverage <- list()
for (i in c("exons","introns","promoters")) {
  cat(i,"\n")
  regiongr <- gene.obj[[i]]
  regiongr$numCGss <- numCGss[[i]]
  genomicCoverage[[i]]$totalRegions <- length(regiongr)
  regCov <- numeric()
  percRegCov <- numeric()
  for (j in 1:length(mcgrlist)) {
    x <- findOverlaps(regiongr, mcgrlist[[ j ]])
    y <- length(unique(queryHits(x)))
    regCov <- c(regCov, y)
    percRegCov <- c(percRegCov, y / length(regiongr) * 100)
  }
  names(regCov) <- names(mcgrlist)
  names(percRegCov) <- names(mcgrlist)
  genomicCoverage[[i]]$regionCoverage <- regCov
  genomicCoverage[[i]]$percRegionCoverage <- percRegCov
}

for (i in c("CpGi","shores")) {
  cat(i,"\n")
  regiongr <- cpg.obj[[i]]
  regiongr$numCGss <- numCGss[[i]]
  genomicCoverage[[i]]$totalRegions <- length(regiongr)
  regCov <- numeric()
  percRegCov <- numeric()
  for (j in 1:length(mcgrlist)) {
    x <- findOverlaps(regiongr, mcgrlist[[ j ]])
    y <- length(unique(queryHits(x)))
    regCov <- c(regCov, y)
    percRegCov <- c(percRegCov, y / length(regiongr) * 100)
  }
  names(regCov) <- names(mcgrlist)
  names(percRegCov) <- names(mcgrlist)
  genomicCoverage[[i]]$regionCoverage <- regCov
  genomicCoverage[[i]]$percRegionCoverage <- percRegCov
}


percRegionCoverageMat <- sapply(genomicCoverage, "[[", "percRegionCoverage")


x <- list()
for (i in 1:ncol(percRegionCoverageMat)) {
  x[[i]] <- percRegionCoverageMat[1:6,i] # no WGBS
}
names(x) <- colnames(percRegionCoverageMat)

bscols <- c("blue1","blue3","salmon1","salmon3","plum","plum4")

if (makepdf) { pdf(paste(outdir, "annotationsRegionCoverage.pdf", sep="/")) }

beeswarm(x,
  labels=names(x),
  main="Region coverage",
  xlab="genomic feature",
  ylab="percent (%)",
  ylim=c(0,100),
  pwcol=rep(bscols,5),
  pch=21,
  pwbg=rep(bscols,5),
  cex=1.5)
legend("bottomright", legend=rownames(percRegionCoverageMat)[1:6], fill=bscols, border=bscols, box.col=rgb(0.6,0.6,0.6,alpha=0.5))

if (makepdf) { dev.off() }





regionCGcoverage <- list()
for (i in c("exons","introns","promoters","CpGi","shores")) {
  cat(i," ")
  regCov <- numeric()
  for (j in 1:length(mcgrlist)) {
    cat(j, " ")
    if (i == "CpGi") {
      regCov <- c(regCov, sum(values(mcgrlist[[j]])[["cpgi"]]) / sum(numCGss[[i]]) * 100)
    } else {
      regCov <- c(regCov, sum(values(mcgrlist[[j]])[[i]]) / sum(numCGss[[i]]) * 100)
    }
  }
  names(regCov) <- names(mcgrlist)
  regionCGcoverage[[i]] <- regCov
  cat("\n")
}



percRegionCGcoverageMat <- numeric()
for (i in 1:length(regionCGcoverage)) {
  percRegionCGcoverageMat <- cbind(percRegionCGcoverageMat, regionCGcoverage[[i]])
}
colnames(percRegionCGcoverageMat) <- names(regionCGcoverage)



x <- list()
for (i in 1:ncol(percRegionCGcoverageMat)) {
  x[[i]] <- percRegionCGcoverageMat[1:6,i] # no WGBS
}
names(x) <- colnames(percRegionCGcoverageMat)

bscols <- c("blue1","blue3","salmon1","salmon3","plum","plum4")

if (makepdf) { pdf(paste(outdir, "annotationsRegionCGUCoverage.pdf", sep="/")) }

beeswarm(x,
  labels=names(x),
  main="Region CGU coverage",
  xlab="genomic feature",
  ylab="percent (%)",
  ylim=c(0,100),
  pwcol=rep(bscols,5),
  pch=21,
  pwbg=rep(bscols,5),
  cex=1.5)
legend("topleft", legend=rownames(percRegionCGcoverageMat)[1:6], fill=bscols, border=bscols, box.col=rgb(0.6,0.6,0.6,alpha=0.5))

if (makepdf) { dev.off() }









mat <- rbind(
  c(1,2), c(1,3), c(1,4), c(1,5), c(1,6),
  c(2,3), c(2,4), c(2,5), c(2,6),
  c(3,4), c(3,5), c(3,6),
  c(4,5), c(4,6),
  c(5,6))
  
rownames(mat) <- paste(names(mcgrlist)[mat[,1]], names(mcgrlist)[mat[,2]],sep=":")
cgOverlap <- list()

for (j in c("exons","introns","promoters","cpgi","shores")) {
  cat(j,"\n")
  a <- apply(mat, 1, function(i) {
    cat(i, "\n")
    grA      <- mcgrlist[[ i[1] ]][ which(values(mcgrlist[[ i[1] ]])[[j]] == 1) ]
    grB      <- mcgrlist[[ i[2] ]][ which(values(mcgrlist[[ i[2] ]])[[j]] == 1) ]
    overlaps <- findOverlaps(grA, grB)
    x <- length(unique(queryHits(overlaps)))
    y <- length(unique(subjectHits(overlaps)))
    
    c(length(grA),
      length(grB), 
      x, 
      y,
      x / length(grA) * 100,
      y / length(grB) * 100)
  })
  rownames(a) <- c("totalA","totalB","numOverlapA","numOverlapB","percOverlapA","percOverlapB")
  colnames(a) <- rownames(mat)

  cgOverlap[[j]] <- t(a)
}


labels  <- rownames(cgOverlap$promoters)
labelsA <- sapply(strsplit(labels, ":"), "[[" ,1)
labelsB <- sapply(strsplit(labels, ":"), "[[" ,2)

omar <- par()$mar

par(mar=c(9, 9, 5.1, 4.1))
plot(cgOverlap$promoters[,5], cgOverlap$promoters[,6],main="promoters",xlab="",ylab="",axes=F)
box()
axis(3)
axis(4,las=2)
axis(1, at=cgOverlap$promoters[,5], labels=labelsA, las=2)
axis(2, at=cgOverlap$promoters[,6], labels=labelsB, las=2)
segments(0,cgOverlap$promoters[,6],cgOverlap$promoters[,5], col="gray")
segments(cgOverlap$promoters[,5],0,y1=cgOverlap$promoters[,6], col="gray")

par(mar=omar)















labels  <- rownames(cgOverlap$promoters)
labelsA <- sapply(strsplit(labels, ":"), "[[" ,2)
labelsB <- sapply(strsplit(labels, ":"), "[[" ,1)

omar <- par()$mar

par(mar=c(9, 9, 5.1, 4.1))
plot(cgOverlap$promoters[,6], cgOverlap$promoters[,5],main="promoters",xlab="",ylab="",axes=F)
box()
axis(3)
axis(4,las=2)
axis(1, at=cgOverlap$promoters[,6], labels=labelsA, las=2)
axis(2, at=cgOverlap$promoters[,5], labels=labelsB, las=2)
segments(0,cgOverlap$promoters[,5],cgOverlap$promoters[,6], col="gray")
segments(cgOverlap$promoters[,6],0,y1=cgOverlap$promoters[,5], col="gray")

par(mar=omar)







pcols <- c("black","blue","red","green","magenta")
names(pcols) <- c("exons","introns","promoters","cpgi","shores")

labels  <- rownames(cgOverlap[[1]])
labelsA <- sapply(strsplit(labels, ":"), "[[" ,2)
labelsB <- sapply(strsplit(labels, ":"), "[[" ,1)

omar <- par()$mar
par(mar=c(9, 9, 5.1, 4.1))

pdf("sampleCGannotOverlap.pdf")

plot(0, type="n", main="Sample CGU annotation overlap", xlab="", ylab="", xlim=c(20,100), ylim=c(20,100))
legend("bottomright", legend=names(pcols),pch=1, col=pcols)
box()
for (j in names(pcols)) {
  points(cgOverlap[[j]][,6], cgOverlap[[j]][,5], col=pcols[j])
  # axis(3)
  # axis(4,las=2)
  # axis(1, at=cgOverlap$promoters[,6], labels=labelsA, las=2)
  # axis(2, at=cgOverlap$promoters[,5], labels=labelsB, las=2)
  # segments(0,cgOverlap$promoters[,5],cgOverlap$promoters[,6], col="gray")
  # segments(cgOverlap$promoters[,6],0,y1=cgOverlap$promoters[,5], col="gray")
}
par(mar=omar)


dev.new()



# dendrogram
# plot(hclust(dist(cgOverlap$promoters[,5:6])))





#mean(cgOverlap[[1]][c(1,10,15),5:6])

intra <- sapply(cgOverlap, function(i) { i[c(1,10,15),5:6] })

x <- list()
for (i in 1:ncol(intra)) { x[[i]] <- intra[,i] }
names(x) <- colnames(intra)

beeswarm(x,col=pcols,ylim=c(80,100))



inter <- sapply(cgOverlap, function(i) { i[-c(1,10,15),5:6] })

x <- list()
for (i in 1:ncol(inter)) { x[[i]] <- inter[,i] }
names(x) <- colnames(inter)

library(beeswarm)

beeswarm(x,col=pcols)


#intra <- sapply(cgOverlap, function(i) { i[c(1,10,15),5:6] })



# install.packages("UpSetR")
library(UpSetR)

load("mcgrlistCGunitsAnnotations.rda")

listInput <- list()
#i=1
for (i in 1:(length(mcgrlist)-1)) { # exclude WGBS
  cat(names(mcgrlist)[i], "\n")
  x <- as.matrix(values(mcgrlist[[i]])[6:10])
  y <- rowSums(x)
  x <- which(y != 0)
  x <- mcgrlist[[i]][x]
  y <- paste(seqnames(x), start(x),sep=".")
  listInput[[ names(mcgrlist)[i] ]] <- y
}

x <- fromList(listInput)
save(x, file="allAnnotatedCGU.upsetr.rda")

  if (makepdf) { pdf(paste(outdir, paste("pairwiseAllAnnotationOverlap.pdf",sep=""), sep="/"), width=11, height=8.5) }
   upset(x, nintersects = 21, sets = names(x), number.angles = 45, mainbar.y.label = "Annotation intersection size", sets.x.label = "Total annotated CGUs")
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir, paste("allAnnotationOverlap.pdf",sep=""), sep="/"), width=11, height=8.5) }
   upset(x, nintersects = 63, sets = names(x), number.angles = 45, mainbar.y.label = "Annotation intersection size", sets.x.label = "Total annotated CGUs")
  if (makepdf) { dev.off() }














listInput <- list()

for (j in c("exons","introns","promoters","cpgi","shores")) {
  cat(j, "\n")
  listInput <- list()
  for (i in 1:length(mcgrlist)) {
    cat(i, names(mcgrlist)[i], "\n")
    x <- mcgrlist[[ i ]][ which(values(mcgrlist[[ i ]])[[j]] == 1) ]  
    y <- paste(seqnames(x), start(x),sep=".")
    listInput[[ names(mcgrlist)[i] ]] <- y
  }
  
  if (makepdf) { pdf(paste(outdir, paste("pairwiseAnnotationOverlap",j,".pdf",sep=""), sep="/"),width=11, height=8.5) }
    upset(fromList(listInput), nintersects = 21, sets = names(listInput)[1:6], number.angles = 45, mainbar.y.label = "Annotation intersection size", sets.x.label = "Number annotated CGUs")
  if (makepdf) { dev.off() }
  
  if (makepdf) { pdf(paste(outdir, paste("allAnnotationOverlap",j,".pdf",sep=""), sep="/"),width=11, height=8.5) }
    upset(fromList(listInput), nintersects = 63, sets = names(listInput)[1:6], number.angles = 45, mainbar.y.label = "Annotation intersection size", sets.x.label = "Number annotated CGUs")
  if (makepdf) { dev.off() }
  
  cat("\n")
}





  # upset(data, nsets = 5, nintersects = 40, sets = NULL,
  #   set.metadata = NULL, intersections = NULL, matrix.color = "gray23",
  #   main.bar.color = "gray23", mainbar.y.label = "Intersection Size",
  #   mainbar.y.max = NULL, sets.bar.color = "gray23",
  #   sets.x.label = "Set Size", point.size = 2.2, line.size = 0.7,
  #   name.size = 7, mb.ratio = c(0.7, 0.3), expression = NULL,
  #   att.pos = NULL, att.color = main.bar.color, order.by = c("freq",
  #   "degree"), decreasing = c(T, F), show.numbers = "yes",
  #   number.angles = 0, group.by = "degree", cutoff = NULL, queries = NULL,
  #   query.legend = "none", shade.color = "gray88", shade.alpha = 0.25,
  #   matrix.dot.alpha = 0.5, empty.intersections = NULL, color.pal = 1,
  #   boxplot.summary = NULL, attribute.plots = NULL)
  
  
  # upset(fromList(listInput),sets = names(listInput)[1:6], number.angles = 90, mainbar.y.label = "Annotation intersection", sets.x.label = "Number annotated CGUs")
  # upset(fromList(listInput),sets = names(listInput)[1:6], number.angles = 45, mainbar.y.label = "Annotation intersection", sets.x.label = "Number annotated CGUs")
  # upset(fromList(listInput), nintersects = 21, sets = names(listInput)[1:6], number.angles = 45, mainbar.y.label = "Annotation intersection size", sets.x.label = "Number annotated CGUs")
  # upset(fromList(listInput), nintersects = 63, sets = names(listInput)[1:6], number.angles = 45, mainbar.y.label = "Annotation intersection size", sets.x.label = "Number annotated CGUs")




intramat <-  rbind(c(1,2), c(3,4), c(5,6))
intraPercOverlap <- apply(intramat,1, function(i) {
  c(length(intersect(listInput[[ i[1] ]], listInput[[ i[2] ]])) / length(listInput[[ i[1] ]]),
  length(intersect(listInput[[ i[1] ]], listInput[[ i[2] ]])) / length(listInput[[ i[2] ]])) 
})



intermat <- rbind(c(1,3), c(1,4), c(1,5), c(1,6),
                  c(2,3), c(2,4), c(2,5), c(2,6),
                  c(3,5), c(3,6),
                  c(4,5), c(4,6))
interPercOverlap <- apply(intermat,1, function(i) {
  cat(i,"\n")
  c(length(intersect(listInput[[ i[1] ]], listInput[[ i[2] ]])) / length(listInput[[ i[1] ]]) * 100,
  length(intersect(listInput[[ i[1] ]], listInput[[ i[2] ]])) / length(listInput[[ i[2] ]]) * 100) 
})



