

# rm(list=ls()); source("/home/thk2008/bin/methylseqplatformcomparison/scripts/newplatformComparison.R")


source("/home/thk2008/bin/methylseqplatformcomparison/scripts/platformComparisonUtils.R")
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/getDataUtils.R")


makepdf <- F


atime  <- format(Sys.time(), "%Y%m%d")
outdir <- paste("results",atime, sep="-")
if (!file.exists(outdir)) {
  dir.create(outdir)
}
cat("output directory is [", outdir, "]\n")

# datasets
#
# ERRBS_A
#  140514_SN250_0692_BC4HF1ACXX	EC-AA-2095	IMR_90 (75ng)
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_140514_SN250_0692_BC4HF1ACXX_EC-AA-2095__uid1313/Project_EC-AA-2095
#   /zenodotus/epicore/scratch/download2/batch022/140514_SN250_0692_BC4HF1ACXX/epigenomics_wcmc/pe-errbs/hg19/IMR_90
# ERRBS_B
#  141003_SN914_0431_BC5GLLACXX	EC-AA-2352	IMR90_A_TruSeq
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_141003_SN914_0431_BC5GLLACXX_EC-AA-2352__uid1886/Project_EC-AA-2352
#   /zenodotus/epicore/scratch/download2/batch022/141003_SN914_0431_BC5GLLACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_A_TruSeq
#
# SSMethylSeq_A optimum
#  140514_SN250_0692_BC4HF1ACXX	EC-AA-2047	IMR90_1 (3ug)
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_140514_SN250_0692_BC4HF1ACXX_EC-AA-2047__uid1315/Project_EC-AA-2047
#   /zenodotus/epicore/scratch/download2/batch022/140514_SN250_0692_BC4HF1ACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_1
# SSMethylSeq_B minimum
#  140801_SN250_0702_AC59NCACXX	EC-AA-2047	IMR90-_2 (1ug)
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_140801_SN250_0702_AC59NCACXX_EC-AA-2047__uid1763/Project_EC-AA-2047
#   /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90-_2
#
# CpGiant_A optimum
#  140801_SN250_0702_AC59NCACXX	EC-AA-2190	IMR90_Roche_1ug
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_140801_SN250_0702_AC59NCACXX_EC-AA-2190__uid1762/Project_EC-AA-2190
#   /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_1ug
# CpGiant_B minimum
#  140801_SN250_0702_AC59NCACXX	EC-AA-2190	IMR90_Roche_0_25ug
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_140801_SN250_0702_AC59NCACXX_EC-AA-2190__uid1762/Project_EC-AA-2190
#   /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_0_25ug
#
# individual
# WGBS_A
#  150115_SN250_0725_BHF3KJADXX  EC-AA-2553  R_WGBS_EpiG_IMR90_100_12_I12
#   /zenodotus/epicore/scratch/sequencing_monitor/store042/demux_182_150115_SN250_0725_BHF3KJADXX_EC-AA-2553__uid2382/Project_EC-AA-2553
#   /zenodotus/epicore/scratch/download2/batch024/150115_SN250_0725_BHF3KJADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12
# WGBS_B
#  150323_SN250_0735_AHGLYKADXX EC-AA-2751 R_WGBS_EpiG_IMR90__100__12__i12
#   /zenodotus/epicore/scratch/sequencing_monitor/store042/demux_182_150323_SN250_0735_AHGLYKADXX_EC-AA-2751__uid2811/Project_EC-AA-2751
#   /zenodotus/epicore/scratch/download2/batch024/150323_SN250_0735_AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_i12
#
# combined
# WGBS_AB
# 150115_SN250_0725_BHF3KJADXX + 150323_SN250_0735_AHGLYKADXX   EC-AA-2553 + EC-AA-2751   R_WGBS_EpiG_IMR90__100__12__i12
# /zenodotus/epicore/scratch/download2/batch025/BHF3KJADXX-AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12


# /zenodotus/epicore/scratch/download2/batch022/140514_SN250_0692_BC4HF1ACXX/epigenomics_wcmc/pe-errbs/hg19/IMR_90
# /zenodotus/epicore/scratch/download2/batch022/141003_SN914_0431_BC5GLLACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_A_TruSeq
# /zenodotus/epicore/scratch/download2/batch022/140514_SN250_0692_BC4HF1ACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_1
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90-_2
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_1ug
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_0_25ug
# individual
# /zenodotus/epicore/scratch/download2/batch024/150115_SN250_0725_BHF3KJADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12
# /zenodotus/epicore/scratch/download2/batch024/150323_SN250_0735_AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_i12
# combined
# /zenodotus/epicore/scratch/download2/batch025/BHF3KJADXX-AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12

# /zenodotus/epicore/scratch/download2/batch022/140514_SN250_0692_BC4HF1ACXX/epigenomics_wcmc/pe-errbs/hg19/IMR_90/methylcall.CpG.IMR_90.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/141003_SN914_0431_BC5GLLACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_A_TruSeq/methylcall.CpG.IMR90_A_TruSeq.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/140514_SN250_0692_BC4HF1ACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_1/methylcall.CpG.IMR90_1.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90-_2/methylcall.CpG.IMR90-_2.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_1ug/methylcall.CpG.IMR90_Roche_1ug.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_0_25ug/methylcall.CpG.IMR90_Roche_0_25ug.mincov0.txt
#
# /zenodotus/epicore/scratch/download2/batch024/150115_SN250_0725_BHF3KJADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12/methylcall.CpG.R_WGBS_EpiG_IMR90_100_12_I12.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch024/150323_SN250_0735_AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_i12/methylcall.CpG.R_WGBS_EpiG_IMR90__100__12__i12.mincov0.txt
#
# /zenodotus/epicore/scratch/download2/batch025/BHF3KJADXX-AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12/methylcall.CpG.R_WGBS_EpiG_IMR90_100_12_I12.mincov0.txt


if (!exists("propsmat")) {
  rda <- "propsmat.rda"
  if (!file.exists(rda)){
    cat("creating properties matrix\n")
    propsmat <- getPropertiesMatrix(summaryfiles)
    save(propsmat,file=rda) # should this also be in results dir?
    savePropsCsv(propsmat,outdir)
  } else {
    cat("loading properties matrix\n")
    load(rda)
  }
}


if (!exists("regionsList")) {
  rda <- "regionsList.rda"
  if (!file.exists(rda)){
    cat("creating regionsList.rda\n")
    source("/home/thk2008/bin/methylseqplatformcomparison/scripts/regionAnalysis.R")
    regionsList <- getRegionsList(regionFiles)
    regionsList[[1]] <- getRegionsTotalCpGs(regionsList[[1]])
    regionsList[[2]] <- getRegionsTotalCpGs(regionsList[[2]])
    regionsList[[3]] <- getRegionsTotalCpGs(regionsList[[3]])
    save(regionsList, file=rda)
  } else {
    cat("loading",rda,"\n")
    load(rda)
  }
}
# names(regionsList)
# [1] "SSMethylSeq"     "CpGiant"   "MspI_84-334"


if (!exists("mcgrlist")) {
  rda <- "mcgrlistCpGsites10x.rda"
  if (!file.exists(rda)){
    cat("creating CpGgrList\n")
    mcgrlist <- getCpGgrList(cpgfiles)
    save(mcgrlist,file=rda)
  } else {
    cat("loading",rda,"\n")
    load(rda)
  }
}
# names(mcgrlist)
# [1] "ERRBS_A"         "ERRBS_B"         "SSMethylSeq_A_opt"   "SSMethylSeq_B_min"
# [5] "CpGiant_A_opt" "CpGiant_B_min" "WGBS"


abtext <- c("A","B","A_opt","B_min","A_opt","B_min","A")


# targetmat
#                 onTarget offTarget
# ERRBS_A          6629668        NA
# ERRBS_B          6696652        NA
# SSMethylSeq_A_opt    2966772    888421
# SSMethylSeq_B_min    2899634    709040
# CpGiant_A_opt  4985497   1042287
# CpGiant_B_min  4713785    948366
# WGBS            12647963        NA

if (!exists("targetmat")) {
  rda <- "targetmat.rda"
  if (!file.exists(rda)){

    cat("creating targetmat\n")
    targetmat <- numeric()
    for (i in 1:length(mcgrlist)) {
      if (i %in% c(1,2,7)){
        targetmat <- rbind(targetmat, c(length(mcgrlist[[i]]),NA))

      } else if (i %in% c(3,4)) { # SSMethylSeq
        x <- findOverlaps(regionsList[[1]], mcgrlist[[i]])
        numoverlapA <- length(unique(queryHits(x)))    # number of regions overlap cpgs
        numoverlapB <- length(unique(subjectHits(x)))  # number of cpgs overlap regions
        totalC <- sum(regionsList[[2]]$numCGss)
        targetmat <- rbind(targetmat, c( numoverlapB, length(mcgrlist[[i]]) - numoverlapB ) )

      } else if (i %in% c(5,6)) { # CpGiant
        x <- findOverlaps(regionsList[[2]], mcgrlist[[i]])
        numoverlapA <- length(unique(queryHits(x)))    # number of regions overlap cpgs
        numoverlapB <- length(unique(subjectHits(x)))  # number of cpgs overlap regions
        totalC <- sum(regionsList[[1]]$numCGds)
        targetmat <- rbind(targetmat, c( numoverlapB, length(mcgrlist[[i]]) - numoverlapB ) )

      } else {
       stop("what is this?")
      }

    }
    colnames(targetmat) <- c("onTarget","offTarget")
    rownames(targetmat) <- names(mcgrlist)

    save(targetmat, file="targetmat.rda")

  } else {
    cat("loading",rda,"\n")
    load(rda)
  }
}


# regionTotalCG
             # SSMethylSeq CpGiant MspI_84-334
# singleStrand     3149681 2812548     3572590
# doubleStrand     6299362 5625096     7145180
regionTotalCG <- sapply(regionsList, function(i) { c(sum(i$numCGss), sum(i$numCGds)) })
rownames(regionTotalCG) <- c("singleStrand","doubleStrand")


# max values
x <- c( regionTotalCG["doubleStrand","MspI_84-334"], 0, regionTotalCG["doubleStrand","MspI_84-334"], 0, # ERRBS
        NA,                # spacer
        regionTotalCG["singleStrand","SSMethylSeq"], 0, regionTotalCG["singleStrand","SSMethylSeq"], 0, # SSMethylSeq
        NA,                # spacer
        regionTotalCG["doubleStrand","CpGiant"], 0, regionTotalCG["doubleStrand","CpGiant"], 0,         # CpGiant
        NA,                # spacer
        28000000)                                                                                       # WGBS

# on target/offtarget values
xx <- as.numeric(t(targetmat))
xx <- c(xx[1:4],
  NA,
  xx[5:8],
  NA,
  xx[9:12],
  NA,
  xx[13])

platform <- c("ERRBS", "SSMethylSeq", "CpGiant", "WGBS")
yaxis    <- pretty(c(0,max(pretty(range(x,na.rm=T)))),n=10)
bcols    <- c(rep("gray90", 5),
              rep("gray50", 10),
              "gray90")

if (makepdf) { pdf(paste(outdir, "numberCpGs.pdf", sep="/"),pointsize=11) }

b <- barplot(x,
 ylim=range(yaxis),
 space=0,
 main="Number of C's in CpG context",
 xlab="platform",
 ylab=expression(paste("number of C's (",10^6,")")),
 col=bcols,
 axes=F)

axis(2, at=yaxis, labels=yaxis/1000000, las=1)

# platform labels
mtext(platform,side=1,at=c(2,7,12,15.5),line=2)
# A, B labels
mtext(abtext,side=1,at=c(0.5,2.5,6,8,11,13),line=1)

bcols <- c(
  rep(c("orange1","orchid1"),2),"white",
  rep(c("lightblue","salmon","lightblue","salmon","white"),2),
  "orange1")

barplot(xx,
 space=0,
 col=bcols,
 add=T,
 axes=F)

legend("topleft",
 legend=c("region's total C's","covered C's in regions","covered C's outside regions","covered C's","predicted total C's"),
 fill=c("gray50","lightblue","salmon","orange1","gray90"))

if (makepdf) { dev.off() } else { dev.new() }


source("/home/thk2008/bin/methylseqplatformcomparison/scripts/qc.R")

if (makepdf) { pdf(paste(outdir, "totalCpGs.pdf", sep="/"),pointsize=11) }
  plotNumCpgs(propsmat[,1], abtext)
if (makepdf) { dev.off() } else { dev.new() }

if (makepdf) { pdf(paste(outdir, "meanCpGcoverage.pdf", sep="/"),pointsize=11) }
  plotMeanCpgCov(propsmat[,2], abtext)
if (makepdf) { dev.off() } else { dev.new() }

if (makepdf) { pdf(paste(outdir, "mappingEfficiency.pdf", sep="/"),pointsize=11) }
  plotMappingEfficiency(propsmat[,3], abtext)
if (makepdf) { dev.off() } else { dev.new() }

if (makepdf) { pdf(paste(outdir, "alignmentStats.pdf", sep="/"),pointsize=11) }
  plotAlignStats(propsmat, abtext)
if (makepdf) { dev.off() } else { dev.new() }

if (makepdf) { pdf(paste(outdir, "alignmentStatsPerc.pdf", sep="/"),pointsize=11) }
  plotAlignStatsPerc(propsmat, abtext)
if (makepdf) { dev.off() } else { dev.new() }


source("/home/thk2008/bin/methylseqplatformcomparison/scripts/cor_qq_maCGunits.R")


plotMADvsDiffOverlap <- function() {

  overlapNmad <- cbind(statsmat[,3] / statsmat[,1] * 100, statsmat[,3] / statsmat[,2] * 100, statsmat[,7])
  overlapNmadNscore <- cbind(overlapNmad, (100 - overlapNmad[,1]) * overlapNmad[,3], (100 - overlapNmad[,2]) * overlapNmad[,3])
  
  par(mar=c(11.1, 4.1, 4.1, 2.1))
  
  plot(overlapNmadNscore[,4],col="red", axes=F, ylim=c(0,80),xlab="",ylab="some score")
  points(overlapNmadNscore[,5],col="blue")
  axis(2)
  axis(1, at=1:nrow(overlapNmadNscore), labels=rownames(overlapNmadNscore), las=2, cex.axis=0.8)
  title(main="combination of correlaion and overlap in to some kind of score")
  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  
  x <- abs(overlapNmadNscore[,1]-overlapNmadNscore[,2])
  y <- overlapNmadNscore[,3]
  z <- cbind(x,y)
  colnames(z) <- c("difference in overlap fraction","MAD value")
  
  # pdf("corVSmad.pdf")
  
  plot(z, col="red",pch=16)
  tmp <- z[which(z[,2] < 0.6),]
  text(tmp, labels=rownames(tmp),srt=90,pos=4,cex=0.7)
  
  tmp <- z[which(z[,2] > 0.6),]
  text(tmp, labels=rownames(tmp),srt=90,pos=2,cex=0.7)
  
  title(main="some measure of the combination of overlap and correlation")
  
  # dev.off()
}




if (F) {
  
  require(GenomicRanges)
  load("regionsList.rda")
  load("mcgrlist_CGunits10x.rda")
  
  rda <- "cgUnitCoverageDistributions.rda"
  if (file.exists(rda)) {
    load("cgUnitCoverageDistributions.rda")
  } else {
    cgUnitCoverageDistributions <- list()
    
    percCoverage <- numeric()
    cat("mspi errbs a\n")
    for (i in 1:length(rl)) {
      x <- findOverlaps(rl[i], mcgrlist[[1]])
      y <- length(unique(subjectHits(x)))
      percCoverage <- c(percCoverage, y / rl[i]$numCGss * 100)
    }
    cgUnitCoverageDistributions[[1]] <- percCoverage
    
    percCoverage <- numeric()
    cat("mspi errbs b\n")
    for (i in 1:length(rl)) {
      x <- findOverlaps(rl[i], mcgrlist[[2]])
      y <- length(unique(subjectHits(x)))
      percCoverage <- c(percCoverage, y / rl[i]$numCGss * 100)
    }
    cgUnitCoverageDistributions[[2]] <- percCoverage
  
    
    rl <- regionsList[[1]] # SSMethylSeq
    percCoverage <- numeric()
    cat("SSMethylSeq a\n")
    for (i in 1:length(rl)) {
      x <- findOverlaps(rl[i], mcgrlist[[3]])
      y <- length(unique(subjectHits(x)))
      percCoverage <- c(percCoverage, y / rl[i]$numCGss * 100)
    }
    cgUnitCoverageDistributions[[3]] <- percCoverage
    
    percCoverage <- numeric()
    cat("SSMethylSeq b\n")
    for (i in 1:length(rl)) {
      x <- findOverlaps(rl[i], mcgrlist[[4]])
      y <- length(unique(subjectHits(x)))
      percCoverage <- c(percCoverage, y / rl[i]$numCGss * 100)
    }
    cgUnitCoverageDistributions[[4]] <- percCoverage
    
    rl <- regionsList[[2]] # CpGiant
    percCoverage <- numeric()
    cat("CpGiant a\n")
    for (i in 1:length(rl)) {
      x <- findOverlaps(rl[i], mcgrlist[[5]])
      y <- length(unique(subjectHits(x)))
      percCoverage <- c(percCoverage, y / rl[i]$numCGss * 100)
    }
    cgUnitCoverageDistributions[[5]] <- percCoverage
    
    percCoverage <- numeric()
    cat("CpGiant b\n")
    for (i in 1:length(rl)) {
      x <- findOverlaps(rl[i], mcgrlist[[6]])
      y <- length(unique(subjectHits(x)))
      percCoverage <- c(percCoverage, y / rl[i]$numCGss * 100)
    }
    cgUnitCoverageDistributions[[6]] <- percCoverage
  
    
    save(cgUnitCoverageDistributions, file="cgUnitCoverageDistributions.rda")
    cat("done\n")
  }
  
}



