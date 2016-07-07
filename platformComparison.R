

# rm(list=ls()); source("/home/thk2008/bin/methylseqplatformcomparison/scripts/platformComparison.R")

# rm(list=ls()); source("/scratch/thk2008/epicore03.bin/methylseqplatformcomparison/scripts/platformComparison.R")

# library(corrplot)
# library(GenomicRanges)
# library(methylKit)
# library(VennDiagram)
# library(Rsamtools)
# library("BSgenome.Hsapiens.UCSC.hg19")
# library(vioplot)
# library(beeswarm)


# pkgs <- c("corrplot","GenomicRanges","methylKit","VennDiagram","BSgenome.Hsapiens.UCSC.hg19","Rsamtools")
# if (!all(pkgs %in% .packages(all=T))) {
  # cat("not all packages are available, ")
  # cat("some analyses may not run.\n")
  # cat("check packages:\n",paste(pkgs,collapse="\n "), "\n")
# }


# if (Sys.info()["nodename"] == "localhost.localdomain") {
  # datadir <- "/scratch/thk2008/methylomecapture/data"
  # source("/scratch/thk2008/epi3bin/methylseqplatformcomparison/scripts/platformComparisonUtils.R")
# }
# if (Sys.info()["nodename"] == "pc142472.med.cornell.edu") {
  # datadir <- "/scratch001/thk2008/methylomecapture/data"
  # source("/home/thk2008/bin/methylseqplatformcomparison/scripts/platformComparisonUtils.R")
# }
# if (Sys.info()["nodename"] == "epicore03.pbtech") {
  datadir <- "/scratch001/thk2008_dat/methylomecapture/data"
  # source("/home/thk2008/bin/methylseqplatformcomparison/scripts/platformComparisonUtils.R")
# }

source("/home/thk2008/bin/methylseqplatformcomparison/scripts/platformComparisonUtils.R")
source("/home/thk2008/bin/methylseqplatformcomparison/scripts/getDataUtils.R")

# makepdf <- switch( menu(c("yes", "no"),title="do you want to make a graphic pdf?") + 1,
 # FALSE,
 # TRUE,
 # FALSE )
makepdf       <- F
doqc          <- T
doregionqc    <- F
docorrelation <- F
doannotations <- F

dooverlap     <- F
domethstats   <- F
doontarget    <- F
dofraglength  <- F


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
# Agilent_A
#  140514_SN250_0692_BC4HF1ACXX	EC-AA-2047	IMR90_1 (3ug)
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_140514_SN250_0692_BC4HF1ACXX_EC-AA-2047__uid1315/Project_EC-AA-2047
#   /zenodotus/epicore/scratch/download2/batch022/140514_SN250_0692_BC4HF1ACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_1
# Agilent_B
#  140801_SN250_0702_AC59NCACXX	EC-AA-2047	IMR90-_2 (1ug)
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_140801_SN250_0702_AC59NCACXX_EC-AA-2047__uid1763/Project_EC-AA-2047
#   /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90-_2
#
# NimbleGen_A
#  140801_SN250_0702_AC59NCACXX	EC-AA-2190	IMR90_Roche_0_25ug
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_140801_SN250_0702_AC59NCACXX_EC-AA-2190__uid1762/Project_EC-AA-2190
#   /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_0_25ug
# NimbleGen_B
#  140801_SN250_0702_AC59NCACXX	EC-AA-2190	IMR90_Roche_1ug
#   /zenodotus/dat01/epicore_scratch/sequencing_monitor/store041/demux_182_140801_SN250_0702_AC59NCACXX_EC-AA-2190__uid1762/Project_EC-AA-2190
#   /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_1ug
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
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_0_25ug
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_1ug
# individual
# /zenodotus/epicore/scratch/download2/batch024/150115_SN250_0725_BHF3KJADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12
# /zenodotus/epicore/scratch/download2/batch024/150323_SN250_0735_AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_i12
# combined
# /zenodotus/epicore/scratch/download2/batch025/BHF3KJADXX-AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12

# /zenodotus/epicore/scratch/download2/batch022/140514_SN250_0692_BC4HF1ACXX/epigenomics_wcmc/pe-errbs/hg19/IMR_90/methylcall.CpG.IMR_90.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/141003_SN914_0431_BC5GLLACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_A_TruSeq/methylcall.CpG.IMR90_A_TruSeq.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/140514_SN250_0692_BC4HF1ACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_1/methylcall.CpG.IMR90_1.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90-_2/methylcall.CpG.IMR90-_2.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_0_25ug/methylcall.CpG.IMR90_Roche_0_25ug.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch022/140801_SN250_0702_AC59NCACXX/epigenomics_wcmc/pe-errbs/hg19/IMR90_Roche_1ug/methylcall.CpG.IMR90_Roche_1ug.mincov0.txt
#
# /zenodotus/epicore/scratch/download2/batch024/150115_SN250_0725_BHF3KJADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12/methylcall.CpG.R_WGBS_EpiG_IMR90_100_12_I12.mincov0.txt
# /zenodotus/epicore/scratch/download2/batch024/150323_SN250_0735_AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_i12/methylcall.CpG.R_WGBS_EpiG_IMR90__100__12__i12.mincov0.txt
#
# /zenodotus/epicore/scratch/download2/batch025/BHF3KJADXX-AHGLYKADXX/epigenomics_wcmc/pewgbs/hg19/R_WGBS_EpiG_IMR90_100_12_I12/methylcall.CpG.R_WGBS_EpiG_IMR90_100_12_I12.mincov0.txt



########
#      #
#  QC  #
#      #
########

if (doqc) {
  cat("qc started... ")
  source("/home/thk2008/bin/methylseqplatformcomparison/scripts/qc.R")

  # the sumamry files should all have identical format
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

  if (makepdf) { pdf(paste(outdir, "numberCpGs.pdf", sep="/")) }
    plotNumCpgs(propsmat[,1])
  if (makepdf) { dev.off() } else { dev.new() }

  if (makepdf) { pdf(paste(outdir, "meanCpGcoverage.pdf", sep="/")) }
    plotMeanCpgCov(propsmat[,2])
  if (makepdf) { dev.off() } else { dev.new() }

  if (makepdf) { pdf(paste(outdir, "mappingEfficiency.pdf", sep="/")) }
    plotMappingEfficiency(propsmat[,3])
  if (makepdf) { dev.off() } else { dev.new() }

  if (makepdf) { pdf(paste(outdir, "alignmentStats.pdf", sep="/")) }
    plotAlignStats(propsmat)
  if (makepdf) { dev.off() } else { dev.new() }

  if (makepdf) { pdf(paste(outdir, "alignmentStatsPerc.pdf", sep="/")) }
    plotAlignStatsPerc(propsmat)
  if (makepdf) { dev.off() } else { dev.new() }

  # if (makepdf) { pdf(paste(outdir, "numberCpGsVsOtherC.pdf", sep="/")) }
    # plotNumCpgvsOtherC(propsmat)
  # if (makepdf) { dev.off() }

  cat(" done.\n")
}
# convert -density 300 -depth 8 -quality 100 x.pdf x.png



if (doregionqc) {
  cat("region qc started... ")
  source("/home/thk2008/bin/methylseqplatformcomparison/scripts/regionAnalysis.R")

  if (!exists("regionsList")) {
    rda <- "regionsList.rda"
    if (!file.exists(rda)){
      cat("creating regionsList.rda\n")
      regionsList <- getRegionsList(regionFiles)
      save(regionsList, file=rda)
    } else {
      cat("loading",rda,"\n")
      load(rda)
    }
  }

  if (!exists("mcgrlist")) {
    rda <- "mcgrlist.rda"
    if (!file.exists(rda)){
      cat("getting CpGgrList\n")
      mcgrlist <- getCpGgrList(cpgfiles)
      save(mcgrlist,file=rda)
    } else {
      cat("loading",rda,"\n")
      load(rda)
    }
  }

  if (!exists("percRegCovered")) {
    rda <- "percRegCovered.rda"
    if (!file.exists(rda)){
      cat("getting percRegCovered\n")
      percRegCovered <- getRegionCoverage(regionsList, mcgrlist)
      save(percRegCovered,file=rda)
    } else {
      cat("loading",rda,"\n")
      load(rda)
    }
  }



  if (makepdf) { pdf(paste(outdir, "regionwidthvioplot.pdf", sep="/")) }
    cat("making violin plot\n")
    makeViolinPlot(regionsList)
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir, "overlappairs.pdf", sep="/"), width=14, height=14) }
    cat("making overlap plot\n")
    makeOverlapPlot1(regionsList)   # can always add label column and row with image editor
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"regionsCovered.pdf",sep="/")) } else { dev.new() }
    plotRegionCoverage(percRegCovered)
  if (makepdf) { dev.off() }

  if (!exists("regions.total.cpg")) {
    rda <- "regions.total.cpg.gr.rda"
    if (!file.exists(rda)){
      cat("computing regions.total.cpg\n")
      regions.total.cpg.gr <- list()
      regions.total.cpg.gr[[1]]   <- getRegionsTotalCpGs(regionsList[[1]], stranded=F)       # NimbleGen
      regions.total.cpg.gr[[2]]   <- getRegionsTotalCpGs(regionsList[[2]], stranded=F, ds=F) # Agilent
      regions.total.cpg.gr[[3]]   <- getRegionsTotalCpGs(regionsList[[3]], stranded=F)       # MspI_84-334
      names(regions.total.cpg.gr) <- names(regionsList)[1:3]
      save(regions.total.cpg.gr, file="regions.total.cpg.gr.rda")
    } else {
      cat("loading",rda,"\n")
      load(rda)
    }
  }

  if (!exists("cpgcovlist")) {
    rda <- "cpgcovlist.rda"
    if (!file.exists(rda)){
      cat("computing cpgcovlist\n")
      pairmat <- matrix(c(1,2,3,4,5,6, 3,3,2,2,1,1), ncol = 2)
      cpgcovlist <- getAnnotCovList(mcgrlist[1:6], regions.total.cpg.gr, pairmat)
      save(cpgcovlist, file=rda)
    } else {
      cat("loading", rda, "\n")
      load(rda)
    }
  }

  if (makepdf) { pdf(paste(outdir,"regionsCpGcoverage.pdf",sep="/")) } else { dev.new() }
    plotAnnotPercCoverageViolin(cpgcovlist, regions.total.cpg.gr, pairmat, what="respective regions")
  if (makepdf) { dev.off() }






  if (!exists("regionsCpGList")) {
    rda <- "regionsCpGList.rda"
    if (!file.exists(rda)){

      # j=1 # nimblegen
      # x <- regionsList[[j]]
      # start(x) <- start(x) - 1
      # end(x)   <- end(x) + 1
      # x <- trim(x)
      # regseqs <- getSeq(Hsapiens, x)
      #
      # y1 <- vmatchPattern("CG", regseqs, fixed=T)
      # yy1 <- start(y1) + start(x) - 1
      # names(yy1) <- seqnames(x)
      # gr1 <- GRanges(seqnames= names(unlist(yy1)), ranges=IRanges(unlist(yy1),end=unlist(yy1)))
      #
      # y2 <-  vmatchPattern("CG", reverseComplement(regseqs), fixed=T)
      # yy2 <- end(x) + 1 - end(y2) + 1
      # names(yy2) <- seqnames(x)
      # gr2 <- GRanges(seqnames= names(unlist(yy2)), ranges=IRanges(unlist(yy2),end=unlist(yy2)))
      #
      # nimblegen.cpggr <- unique(c(gr1,gr2))
      # save(nimblegen.cpggr, file="nimblegen.cpggr.rda")
      #
      # --------
      #
      # j=2 # agilent
      # x <- regionsList[[j]]
      # start(x) <- start(x) - 1
      # end(x)   <- end(x) + 1
      # x <- trim(x)
      # regseqs <- getSeq(Hsapiens, x)
      #
      # y1 <- vmatchPattern("CG", regseqs, fixed=T)
      # yy1 <- start(y1) + start(x) - 1
      # names(yy1) <- seqnames(x)
      # gr1 <- GRanges(seqnames= names(unlist(yy1)), ranges=IRanges(unlist(yy1),end=unlist(yy1)))
      #
      # agilent.cpggr <- unique(gr1)
      # save(agilent.cpggr, file="agilent.cpggr.rda")
      #
      # --------
      #
      # j=3 # mspi
      # x <- regionsList[[j]]
      # start(x) <- start(x) - 1
      # end(x)   <- end(x) + 1
      # x <- trim(x)
      # regseqs <- getSeq(Hsapiens, x)
      #
      # y1 <- vmatchPattern("CG", regseqs, fixed=T)
      # yy1 <- start(y1) + start(x) - 1
      # names(yy1) <- seqnames(x)
      # gr1 <- GRanges(seqnames= names(unlist(yy1)), ranges=IRanges(unlist(yy1),end=unlist(yy1)))
      #
      # y2 <-  vmatchPattern("CG", reverseComplement(regseqs), fixed=T)
      # yy2 <- end(x) + 1 - end(y2) + 1
      # names(yy2) <- seqnames(x)
      # gr2 <- GRanges(seqnames= names(unlist(yy2)), ranges=IRanges(unlist(yy2),end=unlist(yy2)))
      #
      # mspi.cpggr <- unique(c(gr1,gr2))
      # save(mspi.cpggr, file="mspi.cpggr.rda")
      #
      # --------
      #
      # load("nimblegen.cpggr.rda")
      # load("agilent.cpggr.rda")
      # load("mspi.cpggr.rda")
      # regionsCpGList <- list()
      # regionsCpGList[[1]] <- nimblegen.cpggr
      # regionsCpGList[[2]] <- agilent.cpggr
      # regionsCpGList[[3]] <- mspi.cpggr
      # names(regionsCpGList) <- c("NimbleGen","Agilent","MspI_84-334")
      # save(regionsCpGList, file="regionsCpGList.rda")

    } else {
      cat("loading",rda,"\n")
      load(rda)
    }
  }

  ymax <- tail(pretty(c(0, max(sapply(regionsCpGList, length)))),n=1)
  for (i in 1:(length(regionsCpGList)-1)) {
    for (j in (i+1):(length(regionsCpGList))) {
      if (makepdf) { pdf(paste(outdir,paste("regionsCpGoverlap-",names(regionsCpGList)[i],"-",names(regionsCpGList)[j],".pdf",sep=""),sep="/")) } else { dev.new() }
        panel.overlap(regionsCpGList[i], regionsCpGList[j], ymax)
      if (makepdf) { dev.off() }
    }
  }







  # on-/off-target
  ontarget <- sapply(cpgcovlist, sum) / sapply(mcgrlist[1:6], length) * 100
  if (makepdf) { pdf(paste(outdir,"regionsOntargetCpGcoverage.pdf",sep="/")) } else { dev.new() }
    plotOnTarget(ontarget)
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"targetCpGcoverageSwarm.pdf",sep="/")) } else { dev.new() }
    plotRegionCoverageSwarm(percRegCovered, ontarget)
  if (makepdf) { dev.off() }


  cat("region qc done.\n")
}
# convert -density 300 -depth 8 -quality 100 regionwidthvioplot.pdf regionwidthvioplot.png



if (docorrelation) {

  cat("computing correlation\n")
  source("/home/thk2008/bin/methylseqplatformcomparison/scripts/cor_qq_ma.R")

  if (!exists("mcgrlist")) {
    rda <- "mcgrlist.rda"
    if (!file.exists(rda)){
      cat("getting CpGgrList\n")
      mcgrlist <- getCpGgrList(cpgfiles)
      save(mcgrlist,file=rda)
    } else {
      cat("loading",rda,"\n")
      load(rda)
    }
  }

  if (!exists("statsmat")) {
    rda <- "statsmat.rda"
    if (!file.exists(rda)){
      cat("computing statsmat\n")
      statsmat <- getStatsMat(mcgrlist)
      cat("saving statsmat\n")
      save(statsmat,file=rda)
      saveCsv(statsmat, name="statsmat")
    } else {
      cat("loading",rda,"\n")
      load(rda)
    }
  }


  for (i in 1:(length(mcgrlist)-1)) {
    for (j in (i+1):length(mcgrlist)) {
      cat(i,"::",j,"MA-plot\n")
      if (makepdf) { pdf(paste(outdir, paste(names(mcgrlist)[i],names(mcgrlist)[j],"maplot.pdf",sep="-"), sep="/")) }
        makeMAplot(mcgrlist[i],mcgrlist[j])
      if (makepdf) { dev.off() }
    }
  }

  for (i in 1:(length(mcgrlist)-1)) {
    for (j in (i+1):length(mcgrlist)) {
      cat(i,"::",j,"smooth scatter cor plot\n")
      if (makepdf) { pdf(paste(outdir, paste(names(mcgrlist)[i],names(mcgrlist)[j],"smoothscattercor.pdf",sep="-"), sep="/")) }
        makeSmoothScatter(mcgrlist[i],mcgrlist[j])
      if (makepdf) { dev.off() }
    }
  }

  ymax <- tail(pretty(c(0, max(sapply(mcgrlist, length)))),n=1)
  for (i in 1:(length(mcgrlist)-1)) {
    for (j in (i+1):length(mcgrlist)) {
      cat(i,"::",j,"cpg overlap barplot\n")
      if (makepdf) { pdf(paste(outdir, paste(names(mcgrlist)[i],names(mcgrlist)[j],"CpGoverlap.pdf",sep="-"), sep="/")) }
        plotOverlap(mcgrlist[i],mcgrlist[j], ymax)
      if (makepdf) { dev.off() }
    }
  }

# think about including these
# the ecdf may be appropriate
if (F) {
  i=1
  j=2
  a <- mcgrlist[i]
  b <- mcgrlist[j]
  makeMAplot(a,b)
  makeSmoothScatter(a,b)

  commonsites <- findOverlaps(a[[1]], b[[1]])
  x <- a[[1]][ queryHits(commonsites)   ]$numC
  y <- b[[1]][ subjectHits(commonsites) ]$numC

  hist(x,breaks=50)
  dev.new()
  hist(y,breaks=50)
  dev.new()
  plot(density(x,bw=1))
  lines(density(y,bw=1),col="blue")

  plot(ecdf(x))
  plot(ecdf(y),col="blue",add=T)
}



  # if (!exists("overlapsubset")) {
    # overlapsubset <- mcgrlist[[1]]
    # for (i in 2:length(mcgrlist)) {
      # overlapsubset <- subsetByOverlaps(overlapsubset, mcgrlist[[i]])
    # }
  # }
  # if (!exists("mcgrintersect")) {
    # mcgrintersect <- list()
    # for (i in names(mcgrlist)) {
      # mcgrintersect[[i]] <- subsetByOverlaps(mcgrlist[[i]], overlapsubset)
    # }
  # }
  # if (!exists("methMat")) {
    # methMat <- numeric()
    # for (i in mcgrintersect) {
      # methMat <- cbind(methMat,i$numC)
    # }
    # colnames(methMat) <-  names(mcgrlist)
  # }
  # print(cor(methMat))
#
  # pdf (paste(outdir, paste ("intersectionCorPlot", "pdf", sep="."),sep="/"))
  # pairs(methMat,
        # lower.panel=panelSmooth,
        # upper.panel=panelCorSpearman,
        # diag.panel=panelHist,
        # main="CpG intersection methylation correlation (Spearman)")
        # mtext(paste(nrow(methMat), " intersecting CpGs"),line=1)
  # dev.off()
  # cat("made intersectionCorPlot pdf\n")
#
#
  # if (!exists("annotsobj")) {
    # rda <- paste(datadir,"annots.rda",sep="/")
    # if (!file.exists(rda)){
      # annotsobj <- read (as.list(cpgfiles), sample.id=as.list(names(cpgfiles)),assembly="hg19",treatment=rep(1,length(cpgfiles)))
      # save(annotsobj, file=rda)
    # } else {
      # load(rda)
    # }
  # }
  # # destrand: if TRUE, reads covering both strands of a CpG dinucleotide
          # # will be merged, do not set to TRUE if not only interested in
          # # CpGs (default: FALSE). If the methylRawList object contains
          # # regions rather than bases setting destrand to TRUE will have
          # # no effect.
  # meth=unite(annotsobj, destrand=FALSE)
  # pdf (paste(outdir, paste ("methylKitCorrelationPlot", "pdf", sep="."),sep="/"))
   # getCorrelation(meth, plot=T)
  # dev.off()
  # cat("made correlationPlot pdf\n")



  cat("correlation done.\n")

} # end docorrelation



if (doannotations) {

  source("/home/thk2008/bin/methylseqplatformcomparison/scripts/annotations.R")

  if (!exists("mcgrlist")) {
    rda <- "mcgrlist.rda"
    if (!file.exists(rda)){
      cat("getting CpGgrList\n")
      mcgrlist <- getCpGgrList(cpgfiles)
      save(mcgrlist,file=rda)
    } else {
      cat("loading CpGgrList\n")
      load(rda)
    }
  }

  # !!!
  # annotsobj is identical to mcgrlist
  # !!!
  # if (!exists("annotsobj")) {
    # rda <- "annots.rda"
    # if (!file.exists(rda)){
      # cat("creating data object...\n")
      # annotsobj <- getAnnotsObj(cpgfiles)
      # save(annotsobj, file=rda)
    # } else {
      # cat("loading data object [",rda,"] ... \n")
      # load(rda)
    # }
  # }

  if(F) {
    cat("getting annotations\n")
    gene.obj <- read.transcript.features(paste(datadir, "refseq.hg19.bed.txt", sep="/"))
    save(gene.obj, file="gene.obj.rda")
    cpg.obj  <- read.feature.flank(paste(datadir, "cpgis.hg19.bed.txt",sep="/"), feature.flank.name=c("CpGi","shores"))
    save(cpg.obj, file="cpg.obj.rda")
    cat("done getting annotations\n")
  } else {
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

    load("cpg.obj.rda")

    x <- cpg.obj$CpGi
    y <- unique(x)
    cpg.obj$CpGi <- y

    x <- cpg.obj$shores
    y <- unique(x)
    cpg.obj$shores <- y

  }

  if (!exists("annotStats")) {
    rda <- "annotStats.rda"
    if (!file.exists(rda)){
      cat("creating data object...\n")
      annotStats <- getAnnotStats(mcgrlist, gene.obj, cpg.obj)
      save(annotStats, file=rda)
    } else {
      cat("loading data object [",rda,"] ... \n")
      load(rda)
    }
  }

if(F) {
  library(GenomicRanges)
  library(beeswarm)
  x <- lapply(annotStats, "[[", "percGpart")
  y <-sapply(x, c)
  a <- list()
  for (i in 1:nrow(y)) {
    a[[i]] <- y[i,]
  }
  names(a) <- rownames(y)
  beeswarm(a,ylim=c(0,100))
}














  if (makepdf) { pdf(paste(outdir,"annotstatsGeneParts.pdf",sep="/")) } else { dev.new() }
   plotPercentGenePart(annotStats)
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"annotstatsGenomicFeature.pdf",sep="/")) } else { dev.new() }
    plotPercentFeature(annotStats)
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"annotstatsAbsoluteNumbers.pdf",sep="/")) } else { dev.new() }
    plotAnnotsAbsoluteNumbers(annotStats)
  if (makepdf) { dev.off() }

  cat("computing gene part representation\n")
  promoterCoverage <- lapply(mcgrlist, function(i) { # list of promoter regions covered by (overlapping) the CpG sites
    subsetByOverlaps(gene.obj$promoters, i)
  })
  percPromoter <- sapply(promoterCoverage, length) / length(gene.obj$promoters) * 100

  exonCoverage <- lapply(mcgrlist, function(i) { # list of exon regions covered by (overlapping) the CpG sites
    subsetByOverlaps(gene.obj$exons, i)
  })
  percExon <- sapply(exonCoverage, length) / length(gene.obj$exons) * 100

  intronCoverage <- lapply(mcgrlist, function(i) { # list of intron regions covered by (overlapping) the CpG sites
    subsetByOverlaps(gene.obj$introns, i)
  })
  percIntron <- sapply(intronCoverage, length) / length(gene.obj$introns) * 100

  cpgIcoverage <- lapply(mcgrlist, function(i) { # list of CpGi covered by (overlapping) the CpG sites
    subsetByOverlaps(cpg.obj$CpGi, i)
  })
  percCpGi <- sapply(cpgIcoverage, length) / length(cpg.obj$CpGi) * 100

  shoreCoverage <- lapply(mcgrlist, function(i) { # list of CpGi covered by (overlapping) the CpG sites
    subsetByOverlaps(cpg.obj$shores, i)
  })
  percShores <- sapply(shoreCoverage, length) / length(cpg.obj$shores) * 100

  percannotcoverage <- list(percPromoter,percExon,percIntron,percCpGi,percShores)
  names(percannotcoverage) <- c("promoter","exon","intron","CpG islands","shores")
  save(percannotcoverage, file="percannotcoverage.rda")

  if (makepdf) { pdf(paste(outdir,"promoterCoverage.pdf",sep="/")) } else { dev.new() }
   plotAnnotCoverage(percPromoter, what="Promoters")
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"exonCoverage.pdf",sep="/")) } else { dev.new() }
   plotAnnotCoverage(percExon, what="Exons")
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"intronCoverage.pdf",sep="/")) } else { dev.new() }
   plotAnnotCoverage(percIntron, what="Introns")
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"cpgIslandCoverage.pdf",sep="/")) } else { dev.new() }
   plotAnnotCoverage(percCpGi, what="CpG islands")
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"cpgShoresCoverage.pdf",sep="/")) } else { dev.new() }
   plotAnnotCoverage(percShores, what="CpGi shores")
  if (makepdf) { dev.off() }

  if (makepdf) { pdf(paste(outdir,"allPercAnnotCoverage.pdf",sep="/")) } else { dev.new() }
   plotAnnotCoverageAll(percannotcoverage)
  if (makepdf) { dev.off() }




  strandparity <- getStrandParity(mcgrlist)
  if (makepdf)  { pdf(paste(outdir,"strandparity.pdf",sep="/")) } else { dev.new() }
  plotStrandParity(strandparity$strandPerc)
  if (makepdf)  { dev.off() }


# if (makepdf) { pdf(paste(outdir,"annotstatscorplots.pdf",sep="/")) }
#
#
#
  # parts <- rownames(sapply(annotStats, "[[", "numberGpart"))
# for (i in parts) {
  # totpart <- sapply(annotStats, "[[", "numberGpart")[i,]
  # overlappart <- sapply(lapply(annotStats, "[[", "numOverlapGpart"), function(x) { x[,i] })
  # x <- overlappart / totpart * 100
  # if (!makepdf) dev.new()
  # corrplot(x, method = "number",is.corr=F,col="black",tl.pos="d",diag=F,cl.pos="n", mar=c(1,1,2,1))
  # title(main="Percent overlap of genomic part:",cex=0.9,col.main="gray20")
  # mtext(paste(toupper(substring(i, 1, 1)), substring(i, 2), sep="", collapse=" "),line=-1,cex=1,col="darkblue")
# }
#
# parts <- rownames(sapply(annotStats, "[[", "numberFeature"))
# for (i in parts) {
  # totpart <- sapply(annotStats, "[[", "numberFeature")[i,]
  # overlappart <- sapply(lapply(annotStats, "[[", "numOverlapFeature"), function(x) { x[,i] })
  # x <- overlappart / totpart * 100
  # if (!makepdf) dev.new()
  # corrplot(x, method = "number",is.corr=F,col="black",tl.pos="d",diag=F,cl.pos="n", mar=c(1,1,2,1))
  # title(main="Percent overlap of feature:",cex=0.9,col.main="gray20")
  # mtext(paste(toupper(substring(i, 1, 1)), substring(i, 2), sep="", collapse=" "),line=-1,cex=1,col="darkblue")
# }
#
# if (makepdf) dev.off()


#
# convert -density 300 -depth 8 -quality 96 annotstats.pdf[0] annotstats-1.png
# convert -density 300 -depth 8 -quality 96 annotstats.pdf[1] annotstats-2.png
# convert -density 300 -depth 8 -quality 96 annotstats.pdf[2] annotstats-3.png
#
# convert -density 300 -depth 8 -quality 96 annotstatscorplots.pdf[0] annotstatscorplots-1.png
# convert -density 300 -depth 8 -quality 96 annotstatscorplots.pdf[1] annotstatscorplots-2.png
# convert -density 300 -depth 8 -quality 96 annotstatscorplots.pdf[2] annotstatscorplots-3.png
# convert -density 300 -depth 8 -quality 96 annotstatscorplots.pdf[3] annotstatscorplots-4.png
# convert -density 300 -depth 8 -quality 96 annotstatscorplots.pdf[4] annotstatscorplots-5.png
# convert -density 300 -depth 8 -quality 96 annotstatscorplots.pdf[5] annotstatscorplots-6.png
# convert -density 300 -depth 8 -quality 96 annotstatscorplots.pdf[6] annotstatscorplots-7.png
#


} # end doannotations



















#############
#           #
#  Overlap  #
#           #
#############

# cpgfiles <- c(
 # paste(datadir, "cpg.IMR_90.mincov10.txt",         sep="/"),
 # paste(datadir, "cpg.SGUA_2_1.mincov10.txt",       sep="/"),
 # paste(datadir, "cpg.SGUA_2_2.mincov10.txt",       sep="/"),
 # paste(datadir, "cpg.SGUA_2_1_Roche.mincov10.txt", sep="/"),
 # paste(datadir, "cpg.SGUA_2_2_Roche.mincov10.txt", sep="/")
# )
# snames <- sub(".mincov10.txt","",sub("cpg.","",sapply(strsplit(cpgfiles,"/"),tail,n=1)))
# names(cpgfiles) <- snames

if (dooverlap) {
  cat("overlap started...\n")
  if (!exists("mcgrlist")) {
    rda <- paste(datadir,"mcgrlist.rda",sep="/")
    if (!file.exists(rda)){
      mcgrlist <- getCpGgrList(cpgfiles)
      save(mcgrlist,file=rda)
    } else {
      load(rda)
    }
  }

  # get overlap
  if (!exists("overlaps")) {
    rda <- paste(datadir,"overlaps.rda",sep="/")
    if (!file.exists(rda)){
      overlaps <- getOverlapMat(mcgrlist)
      save(overlaps,file=rda)
    } else {
      load(rda)
    }
  }

  for (i in 1:nrow(overlaps$overlapmat)) {
    cat("plotCpGsiteOverlap", sub("::","-",rownames(overlaps$overlapmat)[i]), "\n")
    if (makepdf) { pdf(paste(outdir,paste("overlap",sub("::","-",rownames(overlaps$overlapmat)[i]),".pdf",sep=""),sep="/")) } else { dev.new() }
     plotCpGsiteOverlap(overlaps$overlapmat[i,], rownames(overlaps$overlapmat)[i])
    if (makepdf) { dev.off() }
  }


} # end overlap


if (domethstats) {
  cat("computing methylation and coverage stats... ")
  # get overlap
  if (!exists("annotsobj")) {
    rda <- paste(datadir,"annots.rda",sep="/")
    if (!file.exists(rda)){
      annotsobj <- read (as.list(cpgfiles), sample.id=as.list(names(cpgfiles)),assembly="hg19",treatment=rep(1,length(cpgfiles)))
      save(annotsobj, file=rda)
    } else {
      load(rda)
    }
  }
  for (i in 1:length(annotsobj)) {
    pdf (paste(outdir, paste (names(cpgfiles)[i], "methstats", "pdf", sep="."),sep="/"))
     getMethylationStats(annotsobj[[i]],plot=T,both.strands=F)
     getCoverageStats(annotsobj[[i]],plot=T,both.strands=F)
    dev.off()
  }
  cat("done.\n")
}



if (doontarget) {
  require(GenomicRanges)

  if (!exists("regiongrlist")) {
    rda <- paste(datadir,"regiongrlist.rda",sep="/")
    if (!file.exists(rda)){
      regiongrlist <- list()
      cat("getting region gr...")
      for (i in names(rfiles)) {
        cat(" ", i)
        regiongrlist[[i]] <- getGrRegion(rfiles[i])
      }
      save(regiongrlist,file=rda)
    } else {
      load(rda)
    }
  }
  cat("... got regions\n")


  if (!exists("mcgrlist")) {
    rda <- paste(datadir,"mcgrlist.rda",sep="/")
    if (!file.exists(rda)){
      mcgrlist <- getCpGgrList(cpgfiles)
      save(mcgrlist,file=rda)
    } else {
      load(rda)
    }
  }

  for (i in 1:length(mcgrlist)) {
    for (j in 1:length(regiongrlist)) {
      if (makepdf) { pdf(paste(outdir,paste("binnedregions",names(mcgrlist)[i],names(regiongrlist)[j],"pdf",sep="."),sep="/")) }
      plotBinnedRegionsRecovered(mcgrlist[i], regiongrlist[j])
      if (makepdf) {
        dev.off()
      } else {
        dev.new()
      }
    }
  }

  if (makepdf) { pdf(paste(outdir,"regionsrecovered.pdf",sep="/")) }
   plotRegionsRecovered(mcgrlist, regiongrlist)
  if (makepdf) {
    dev.off()
  } else {
    dev.new()
  }



} # end doontarget


if (dofraglength) {

  for (i in 1:length(bamfiles)) {
    isize <- ""
    if (!exists("isize")) {
      rda <- paste(bamfiles[i],"isize.rda",sep="_")
      if (!file.exists(rda)){
        cat("getting", names(bamfiles)[i], "data")
        isize <- getBamFragmentLengths(bamfiles[i])
        save(isize, file=rda)
        cat("saved ", rda, "\n")
      } else {
        cat("loading ", rda, "\n")
        load(rda)
      }
    }

    isize <- isize[[1]]$isize[ which(isize[[1]]$isize > 0) ]

    if (makepdf) {
      pdf(paste(outdir,paste("fragmentsizedensity-",names(bamfiles)[i],".pdf",sep=""),sep="/"))
    } else {
      dev.new()
    }
    cat(", plotting density")
    plotFragmentSizeDensity(isize, names(bamfiles)[i])

    if (makepdf)  {
      dev.off()
      pdf(paste(outdir,paste("fragmentsizehist-",names(bamfiles)[i],".pdf",sep=""),sep="/"))
    } else {
      dev.new()
    }
    cat(", plotting hist")
    plotFragmentSizeHist(isize, names(bamfiles)[i])

    if (makepdf)  { dev.off() }
    cat(", done\n")
  }
} # end dofraglength



# mc focus
# x <- findOverlaps(mcgrlist[[1]], regiongrlist[[1]])
# length(queryHits(x))           4459273
# length(unique(queryHits(x)))   4459268











if(F) {
# chr1    249250621
# chr2    243199373
# chr3    198022430
# chr4    191154276
# chr5    180915260
# chr6    171115067
# chr7    159138663
# chr8    146364022
# chr9    141213431
# chr10   135534747
# chr11   135006516
# chr12   133851895
# chr13   115169878
# chr14   107349540
# chr15   102531392
# chr16   90354753
# chr17   81195210
# chr18   78077248
# chr19   59128983
# chr20   63025520
# chr21   48129895
# chr22   51304566
# chrM    16571
# chrX    155270560
# chrY    59373566

  # library(Rsamtools)
  # datadir <- "/scratch/thk2008/methylomecapture/data"
#
  # bamfile <- paste(datadir, "IMR_90.sorted.bam", sep="/")
  #
  # param <- ScanBamParam(which=GRanges("chr1", IRanges(1, 150000)),
                        # what=c('strand','pos','seq','qwidth','mpos','mrnm','isize'))
  # tmp <- scanBam(bamfile, param = param, use.name=T)
  #
 #
 # library(GenomicAlignments)
 # param <- ScanBamParam(which=GRanges("chr1", IRanges(1,150000)))
 # tmp <- readGAlignmentsFromBam(bamfile, param=param, use.name=T)
  #



 #
 #
 #
  # [1] "bamfile"
# > bamfile
# [1] "data/IMR_90.sorted.bam"
# > getwd()
# [1] "/scratch001/thk2008_dat/methylomecapture"
#
#
 #
 #
  # bamfile <- "data/IMR_90.sorted.bam"
  #
  # param <- ScanBamParam(which=GRanges("chr1", IRanges(1, 150000)), what=c('isize'))
  # tmp <- scanBam(bamfile, param = param, use.name=T)
#
  # hist(abs(unlist(tmp)))
  # plot(density(abs(unlist(tmp))))
 #
  # param <- ScanBamParam(which=GRanges("chr1", IRanges(1, 249250621)), what=c('isize'))
  # tmp <- scanBam(bamfile, param = param, use.name=T)
#
  # tmp <- unlist(tmp)
  # tmp <- tmp[tmp > 0] #
  # hist(tmp)
  # plot(density(tmp))
 #
  #
  #




# D7ZQJ5M1:692:C4HF1ACXX:8:2311:8170:46284 1:N:0:CTTGT2   179     chr1    11960   255     101M    =       12035   175     CAAACTTTTAAAAAATCACAAAATCTTAATACTATAATCTTCATCTACAAATATCTAACTTCCAACAACTACTAACCTATACCAAAATACAAACTAAACAC   @@@D;2ADFFFFBGDDEE?FEEEDFF<FGDD?BD?DGC@FBFGCFFEFGIFFFFIFF:B?=CDGGDFF@DFF>E>EAE:=C>?D>?=ACCBC>?@BB@AA<   NM:i:34 MD:Z:1AAA1T3A1A1AA5AAA4A2A2A1AA9A2AA1A3A7A5A2AA3A1A3AAA1A3A2A1A3        XM:Z:.zxh.....h.h.hh.....xhh....h..h..x.hh.........x..xh.h...x.......x.....x..xh...x.h...xhh.h...h..x.h...
# D7ZQJ5M1:692:C4HF1ACXX:8:2311:8170:46284 1:N:0:CTTGTA   115     chr1    12035   255     101M    =       11960   -175    CCTATACCAAAATACAAACTAAACACTAAAATAAAATTTTCCTATAAAAAAAAACCATACCTAAAATAAAATAAACCATTATTCATCTTCTAACCCCTATT   ?<?BD?DD8@:DAC4A@EA>C<9AA+AEEEEEEI*?DFID;<B9?:;BDDDIIIID=A3))).;3;AAA>AAA(5;5:555>>AADA5;>;>AA;+:4>>A   NM:i:35 MD:Z:3A1A3AAA1A3A2A1A4AA1A1AA1A7A1AA1A1AA1A4A4A1A1AAA2AAA5A10AA5A2      XM:Z:..x.....hx..........h.....hhh..hhh.h.h....h....h.hh.h.hh.x.......h.hh.h.hx....h.x..h...h.hhx...h.x...

# D7ZQJ5M1:692:C4HF1ACXX:8:1104:3069:47817 1:N:0:CTTGTA   67      chr1    12355   255     101M    =       12487   232     CGGTTGGAGGGAGGGGTTTAGTAGGTTTGGTTTTGGTTTTGGGAGAGTAGGTGGAAGATTAGGTAGGTTATTGTTGTTATAGAATTTAGTGGATTGGTTTA   BB@FFFFFHHHFHIJJGFHFFGIIJJGHFHFGHIFFFHHJCEEHHFHHFFFDFDDEEEDDDCDCDDDDCDDDDDDADDDEEDDDDDDDDCDDDDDDDDDDC   NM:i:24 MD:Z:3T12T1T2T4T3T5TTT8T11T3T3TT2T1T2TT1T4TTT10TT2      XM:Z:Z..x............h.x..x....x...h.....hhx........x...........x...x...hh..z.x..hh.x....hhx..........hh..
# D7ZQJ5M1:692:C4HF1ACXX:8:1104:3069:47817 1:N:0:CTTGT2   131     chr1    12487   255     101M    =       12355   -232    TGGGTGGTAGGTGTAGAGACGGGAGGGGTAGAGTTGTAGGTATAGTTAAGAGGGTTGAAGAAATGGTAGAACGGAGTAGTTGGTGATGTGTGGGTTTATTG  BB@DDDEFHHHHHJJJJJJJJJJJJJJJJJIJJJJJJIIJJGHIJJJJJJBHHIJIJJJJJJJJJIIJJHHHHHHFFFDDDDBDDDDDCCCCDDCBDDDDD   NM:i:17 MD:Z:13T14T4TT1T3T1T2TT7T21T2T14TTT1TT1 XM:Z:.zx.hhh..............x..x....Z................x.......hh..x.h...x.zx....x........Z.....x.............

# D7ZQJ5M1:692:C4HF1ACXX:8:1205:6102:26278 1:N:0:CTTGTA   67      chr1    12355   255     101M    =       12487   232     CGGTTGGAGGGAGGGGTTTAGTAGGTTTGGTTTTGGTTTTGGGAGAGTAGGTGGAAGATTAGGTAGGTTATTGTTGTTATAGAATTTAGTGGATTGGTTTA   BBBFFFFFHGHGHIJIECGHGGIJIIHHEHIHIJFFFFHIFEBHHHHHFFFDCCDEEECDDDDCCDDDCDDDDDDBDDDEEDDDDDDDD@CDDDDCDCBDC   NM:i:24 MD:Z:3T12T1T2T4T3T5TTT8T11T3T3TT2T1T2TT1T4TTT10TT2      XM:Z:Z..x............h.x..x....x...h.....hhx........x...........x...x...hh..z.x..hh.x....hhx..........hh..
# D7ZQJ5M1:692:C4HF1ACXX:8:1205:6102:26278 1:N:0:CTTGT2   131     chr1    12487   255     101M    =       12355   -232    TGGGTGGTAGGTGTAGAGACGGGAGGGGTAGAGTTGTAGGTATAGTTAAGAGGGTTGAAGAAATGGTAGAACGGAGTAGTTGGTGATGTGTGGGTTTATTG  CCCFFFFFHHHHHJJJJJJJJJJIJJJJJJJJJJ*BHIJJJJJJJIIIJJBHIJIIJJJJJJJIJJJJJHHHHHHFFFDDDDDDDDDEDDDDDDCDDDDDD   NM:i:17 MD:Z:13T14T4TT1T3T1T2TT7T21T2T14TTT1TT1 XM:Z:.zx.hhh..............x..x....Z................x.......hh..x.h...x.zx....x........Z.....x.............

# D7ZQJ5M1:692:C4HF1ACXX:8:1213:9905:10333 1:N:0:CTTGTA   67      chr1    12355   255     101M    =       12487   232     CGGTTGGAGGGAGGGGTTTAGTAGGTTTGGTTTTGGTTTTGGGAGAGTAGGTGGAAGATTAGGTAGGTTATTGTTGTTATAGAATTTAGTGGATTGGTTTA   BBBFDDFDFHHDFIIJ?CFGHHIIJIGHGHGFHIDFFGHIFGHGHEEEDFFEFDDEEDCDCCDDDDDC>:ACDDDBDDDEEDDDDDDDDDDDDDDDDDDDC   NM:i:24 MD:Z:3T12T1T2T4T3T5TTT8T11T3T3TT2T1T2TT1T4TTT10TT2      XM:Z:Z..x............h.x..x....x...h.....hhx........x...........x...x...hh..z.x..hh.x....hhx..........hh..
# D7ZQJ5M1:692:C4HF1ACXX:8:1213:9905:10333 1:N:0:CTTGT2   131     chr1    12487   255     101M    =       12355   -232    TGGGTGGTAGGTGTAGAGACGGGAGGGGTAGAGTTGTAGGTATAGTTAAGAGGGTTGAAGAAATGGTAGAACGGAGTAGTTGGTGATGTGTGGGTTTATTG  :?@BDDDFHHHHHJIJJIJJJJJJJJJJJHCGIJJJJIJIIDGIJIIIIJ?CGIIIJJJJJJJIIIJGIHHHHHHFFDDDDDBDDCDCDACDDDDDDDDDD   NM:i:17 MD:Z:13T14T4TT1T3T1T2TT7T21T2T14TTT1TT1 XM:Z:.zx.hhh..............x..x....Z................x.......hh..x.h...x.zx....x........Z.....x.............













}


cat("\ndone.\n")



if(F) {
# regions recovered

library(GenomicRanges)
load("data/mcgrlist.rda")
load("data/regiongrlist.rda")
  regrecovered <- numeric()
  w <- width(regiongrlist[[1]])
  x <- findOverlaps(regiongrlist[[1]], mcgrlist[[1]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[1]]) * 100)
  ssregions <- regiongrlist[[1]][which(w>=50 & w<=500)]
  x <- findOverlaps(ssregions, mcgrlist[[1]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(ssregions) * 100)
  ssregions <- regiongrlist[[1]][which(w>=84 & w<=334)]
  x <- findOverlaps(ssregions, mcgrlist[[1]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(ssregions) * 100)
  x <- findOverlaps(regiongrlist[[1]], mcgrlist[[2]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[1]]) * 100)
  ssregions <- regiongrlist[[1]][which(w>=50 & w<=500)]
  x <- findOverlaps(ssregions, mcgrlist[[2]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(ssregions) * 100)
  ssregions <- regiongrlist[[1]][which(w>=84 & w<=334)]
  x <- findOverlaps(ssregions, mcgrlist[[2]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(ssregions) * 100)
  x <- findOverlaps(regiongrlist[[2]], mcgrlist[[3]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[2]]) * 100)
  x <- findOverlaps(regiongrlist[[2]], mcgrlist[[4]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[2]]) * 100)
  x <- findOverlaps(regiongrlist[[3]], mcgrlist[[5]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[3]]) * 100)
  x <- findOverlaps(regiongrlist[[3]], mcgrlist[[6]])
  regrecovered <- c(regrecovered, length(unique(queryHits(x))) / length(regiongrlist[[3]]) * 100)
  names(regrecovered) <- c(paste(names(mcgrlist)[1],"full", names(regiongrlist)[1]),
                           paste(names(mcgrlist)[1], names(regiongrlist)[1], "50-500bp"),
                           paste(names(mcgrlist)[1], names(regiongrlist)[1], "84-334bp"),
                           paste(names(mcgrlist)[2],"full", names(regiongrlist)[1]),
                           paste(names(mcgrlist)[2], names(regiongrlist)[1], "50-500bp"),
                           paste(names(mcgrlist)[2], names(regiongrlist)[1], "84-334bp"),
                           paste(names(mcgrlist)[3], names(regiongrlist)[2]),
                           paste(names(mcgrlist)[4], names(regiongrlist)[2]),
                           paste(names(mcgrlist)[5], names(regiongrlist)[3]),
                           paste(names(mcgrlist)[6], names(regiongrlist)[3]))







  par(mar=c(12,4,4,2))
  b <- barplot(regrecovered,
  main="Region recovery",
  ylab="percent (%)",
  xlab="",
  ylim=c(0,100),
  las=2)
  #title(xlab="samples and region datasets",line=10)
  text(b,regrecovered,round(regrecovered,digits=1),pos=1)


x <- rbind(regrecovered[c(1,2,3,7,9)],regrecovered[c(4,5,6,8,10)])
colnames(x) <- NULL
cols <- c(rep("lightblue",6),rep("salmon",2),rep("plum",2))

if (makepdf) { pdf(paste(outdir,"regionsrecovered.pdf",sep="/")) }
par(mar=c(6,4,4,2))
b <- barplot(x,
  main="Region recovery",
  ylab="percent recovered (%)",
  beside=T,ylim=c(0,100),
  col=cols)
text(b,x,round(x,digits=1),pos=1)
text(b,0,c("A","B"),pos=3)
text(cex=1, x=colMeans(b)+0.25, y=-2, c("ERRBS full MpsI","ERRBS 50-500bp","ERRBS 84-334bp","Agilent","NimbleGen"), xpd=TRUE, srt=45, adj=c(1,0))

if (makepdf) { dev.off() } else { dev.new() }
par(mar=c(5.1,4.1,4.1,2.1))

}


if(F) {
# on/off target

library(GenomicRanges)


# names(regiongrlist)
# [1] "MspI"      "Agilent" "NimbleGen"

# names(mcgrlist)
# [1] "ERRBS_A"     "ERRBS_B"     "Agilent_A" "Agilent_B" "NimbleGen_A"
# [6] "NimbleGen_B"   "WGBS"


mat <- matrix(NA, ncol=length(mcgrlist),nrow=2,dimnames=list(c("onTarget","offTarget"),names(mcgrlist)))


j=1
load("data/recreateFragments/imr90frags.rda")
regions  <- reduce(imr90frags)
i        <- mcgrlist[[j]]
overlaps <- findOverlaps(i, regions)
mat[1,j] <- length(unique(queryHits(overlaps)))  # sites in regions
mat[2,j] <- length(i[-unique(queryHits(overlaps))]) # sites not in regions
rm(imr90frags)

j=2
load("data/recreateFragments/imr90Afrags.rda")
regions <- reduce(imr90Afrags)
i        <- mcgrlist[[j]]
overlaps <- findOverlaps(i, regions)
mat[1,j] <- length(unique(queryHits(overlaps)))  # sites in regions
mat[2,j] <- length(i[-unique(queryHits(overlaps))]) # sites not in regions
rm(imr90Afrags)



j=1
load("data/recreateFragments/imr90frags.rda")
load("data/recreateFragments/imr90Afrags.rda")
x <- c(imr90frags,imr90Afrags)
rm(imr90frags)
rm(imr90Afrags)
regions  <- reduce(x)
i        <- mcgrlist[[j]]
overlaps <- findOverlaps(i, regions)
mat[1,j] <- length(unique(queryHits(overlaps)))  # sites in regions
mat[2,j] <- length(i[-unique(queryHits(overlaps))]) # sites not in regions
j=2
i        <- mcgrlist[[j]]
overlaps <- findOverlaps(i, regions)
mat[1,j] <- length(unique(queryHits(overlaps)))  # sites in regions
mat[2,j] <- length(i[-unique(queryHits(overlaps))]) # sites not in regions





j=3
k=2
i         <- mcgrlist[[j]]
regions   <- regiongrlist[[k]]
overlaps  <- findOverlaps(i, regions)
mat[1,j]  <- length(unique(queryHits(overlaps)))  # sites in regions
mat[2,j] <- length(i[-unique(queryHits(overlaps))]) # sites not in regions

j=4
k=2
i         <- mcgrlist[[j]]
regions   <- regiongrlist[[k]]
overlaps  <- findOverlaps(i, regions)
mat[1,j]  <- length(unique(queryHits(overlaps)))  # sites in regions
mat[2,j] <- length(i[-unique(queryHits(overlaps))]) # sites not in regions


j=5
k=3
i         <- mcgrlist[[j]]
regions   <- regiongrlist[[k]]
overlaps  <- findOverlaps(i, regions)
mat[1,j]  <- length(unique(queryHits(overlaps)))  # sites in regions
mat[2,j] <- length(i[-unique(queryHits(overlaps))]) # sites not in regions

j=6
k=3
i         <- mcgrlist[[j]]
regions   <- regiongrlist[[k]]
overlaps  <- findOverlaps(i, regions)
mat[1,j]  <- length(unique(queryHits(overlaps)))  # sites in regions
mat[2,j] <- length(i[-unique(queryHits(overlaps))]) # sites not in regions

j=7
k=1
i         <- mcgrlist[[j]]
regions   <- regiongrlist[[k]]
overlaps  <- findOverlaps(i, regions)
mat[1,j]  <- length(unique(queryHits(overlaps)))  # sites in regions
mat[2,j] <- length(i[-unique(queryHits(overlaps))]) # sites not in regions


yaxis <- pretty(range(mat))

pdf("onofftargetcpgs.pdf")

par(mar=c(7,5,4,2))
b <- barplot(mat,
 main="Proportions of on/off target CpG's",
 ylim=range(yaxis),
 ylab=expression(paste("number of CpG's (",10^6,")")),
 col=c("lightblue","salmon"),
 axes=F,las=2)
axis(2,at=yaxis,labels=yaxis/1000000,las=2)


dev.off()





}


if(F) {


mcgr <- mcgrlist
regiongr <- regiongrlist


mat <- matrix(c(1,1, 2,1, 3,2, 4,2, 5,3, 6,3), ncol=2, byrow=T)
y <- numeric()
y1 <- numeric()
y2 <- numeric()
for (i in 1:nrow(mat)){
  print(i)
  x <- findOverlaps(regiongr[[mat[i,2]]], mcgr[[mat[i,1]]])
  xx <- regiongr[[mat[i,2]]][unique(queryHits(x))]
  x <- width(xx)
  recovbins <- c(
                length(x[x>0 & x<50]),
                length(x[x>=50 & x<80]),
                length(x[x>=80 & x<100]),
                length(x[x>=100 & x<200]),
                length(x[x>=200 & x<400]),
                length(x[x>=400 & x<800]),
                length(x[x>=800 & x<1000]),
                length(x[x>=1000 & x<2000]),
                length(x[x>=2000 & x<4000]),
                length(x[x>=4000 & x<5000]),
                length(x[x>5000])
                )
  names(recovbins) <- c("1-49","50-79","80-99","100-199","200-399",
                        "400-799","800-999","1K-2K","2K-4K","4K-5K",
                        ">5K")
  y1 <- rbind(y1, recovbins)

  x <- width(regiongr[[mat[i,2]]])
  allbins <- c(
                length(x[x>0 & x<50]),
                length(x[x>=50 & x<80]),
                length(x[x>=80 & x<100]),
                length(x[x>=100 & x<200]),
                length(x[x>=200 & x<400]),
                length(x[x>=400 & x<800]),
                length(x[x>=800 & x<1000]),
                length(x[x>=1000 & x<2000]),
                length(x[x>=2000 & x<4000]),
                length(x[x>=4000 & x<5000]),
                length(x[x>5000])
                )
  y2 <- rbind(y2, allbins)
  perc <- (recovbins/allbins) * 100
  y <- rbind(y,perc)
}
rownames(y1) <- names(mcgrlist)

cols=c("lightblue","lightblue3","pink","pink2","plum","plum2")
par(mar=c(5,5,4,2))
barplot(y2/1000,beside=T,xlab="region length (bp)",ylab=expression(paste("counts (",10^3, ")")), col="white",las=2)
barplot(y1/1000,beside=T,col=cols,add=T,axes=F)
title(main="Recovery of platform regions")
legend("topright",legend=rownames(y1),fill=cols)
par(mar=c(5.1,4.1,4.1,2.1))


if (makepdf) { pdf(paste(outdir,"platformRegionsRecovered.pdf",sep="/")) }

par(mar=c(5,5,4,2))
barplot(y2/1000,beside=T,xlab="region length (bp)",ylab=expression(paste("counts (",10^3, ")")), col="white",las=2)
barplot(y1/1000,beside=T,col=cols,add=T,axes=F)
title(main="Recovery of platform regions")
legend("topright",legend=rownames(y1),fill=cols)
par(mar=c(5.1,4.1,4.1,2.1))

if (makepdf) { dev.off() }



}


