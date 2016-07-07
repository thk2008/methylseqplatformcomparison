
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/scratch.R")

require(GenomicRanges)
# require("BSgenome.Hsapiens.UCSC.hg19")
# require(beeswarm)
# require(vioplot)

# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/regionAnalysis.R")
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/annotations.R")
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/cor_qq_ma.R")
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/getDataUtils.R")
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/overlap.R")
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/platformComparison.R")
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/platformComparisonUtils.R")
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/qc.R")
# source("/home/thk2008/bin/methylseqplatformcomparison/scripts/scratchAnnotations.R")


# load("regionsList.rda")
# names(regionsList)
# [1] "SSMethylSeq" "CpGiant"     "MspI_84-334"


# load("mcgrlist_CGunits10x.rda")
load("mcgrlistCpGsites10x.rda")

# names(mcgrlist)
# [1] "ERRBS_A"           "ERRBS_B"           "SSMethylSeq_A_opt"
# [4] "SSMethylSeq_B_min" "CpGiant_A_opt"     "CpGiant_B_min"
# [7] "WGBS"


discordantPerc <- numeric()

for (i in 1:length(mcgrlist)) {
  #i=1
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
  
  
  tmp <- c(totalsites, length(xplus), length(xminus), length(overlaps))
  
  x1 <- xplus[queryHits(overlaps)]
  x2 <- xminus[subjectHits(overlaps)]
  
  
  d <- cbind( x1$numC / x1$coverage * 100,
              x1$numT / x1$coverage * 100,
              x2$numC / x2$coverage * 100,
              x2$numT / x2$coverage * 100)
  colnames(d) <- c("pC1","pT1","pC2","pT2")
  
  # which are exactly discordant:  c1=0 and t2=100, c1=100 and t2=0
  # these are the same:  c2=0 and t1=100, c2=100 and t1=0
  
  dd <- which( (d[,1] == 100 & d[,4] == 100) | (d[,1] == 0 & d[,4] == 0))
  #head(d[dd,],n=10)
  
  # exactly discordant
  #(length(dd)*2) / totalsites * 100   # 3.84692 %
  
  tmp <- c(tmp, length(dd), (length(dd)*2) / totalsites * 100)
  
  dd1 <- which( (d[,1] >= 90 & d[,4] >= 90) | (d[,1] <= 10 & d[,4] <= 10))
  #head(d[dd1,],n=10)
  
  #(length(dd1)*2) / totalsites * 100  # 6.489375
  
  tmp <- c(tmp, length(dd1), (length(dd1)*2) / totalsites * 100)
  
  
  dd1 <- which( (d[,1] >= 80 & d[,4] >= 80) | (d[,1] <= 20 & d[,4] <= 20))
  #head(d[dd1,],n=10)
  
  # (length(dd1)*2) / totalsites * 100  # 8.826113
  
  tmp <- c(tmp, length(dd1), (length(dd1)*2) / totalsites * 100)
  
  discordantPerc <- rbind(discordantPerc, tmp)
  print(discordantPerc)
  
}
colnames(discordantPerc) <- c("totalSites","totalF","totalR","totalConsecutive",
                              "numberDiscordant0-100","percDiscordant0-100",
                              "numberDiscordant10-90","percDiscordant10-90",
                              "numberDiscordant20-80","percDiscordant20-80")
rownames(discordantPerc) <- names(mcgrlist)
save(discordantPerc,file="discordantPerc.rda")

bcols <- c("red","plum","pink")
omar <- par()$mar

pdf("percentDiscordantSites.pdf")

par(mar=c(9.1,4.1,4.1,2.1))
barplot(t(discordantPerc[,c(6,8,10)]),
  beside=T,
  ylim=c(0,10),
  col=bcols,
  main="fraction discordant sites",
  ylab="percent of total sites",
  las=2)

legend("topright",c("exactly discordant 0 vs 100%","discordant <=10 vs >=90%","discordant <=20 vs >=80%"),
fill=bcols)

par(mar=omar)


dev.off()






    x <- x$numC / x$coverage * 100
    y <- y$numC / y$coverage * 100
    
    labelA <- "plus"
    labelB <- "minus"
  
    corvals <- c(cor(x, y, method="pearson"),
                 cor(x, y, method="spearman"))


   colramp = colorRampPalette(c("white", blues9))

    #smoothScatter(x, y, xlab=labelA, ylab=labelB, nbin=512, colramp=colramp,nr)
    smoothScatter(x, y, xlab=labelA, ylab=labelB, colramp=colorRampPalette(c("white", "blue")),transformation = function(x) x^.9 )

    x[x == 0] <- 0.01
    y[y == 0] <- 0.01
    zz <- cbind(x,y)
    M  <- log2(zz[, 1]) - log2(zz[, 2])
    A  <- rowMeans(log2(zz))


   colramp = colorRampPalette(c("white", blues9))

ma.plot( A, M, cex=1,add.loess=F,show.statistics=F,plot.method="smoothScatter",transformation = function(x) x^.75)






