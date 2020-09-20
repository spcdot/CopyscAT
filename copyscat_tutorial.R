library("devtools")

#STEP 1: Set your working directory to one above the one where you downloaded the package into
#if you installed into a directory named something other than copyscat, replace the directory in install with that
setwd("~")
install("Copyscat")

#load the package
library(CopyscAT)

#first time, set up the references
#load your favourite genome
library(BSgenome.Hsapiens.UCSC.hg38)
#use it to save references - this will create a genome chrom.sizes file, a cytobands file and a CpG file
generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38",tileWidth = 1e6,outputDir = "~")

#initialize the environment 
initialiseEnvironment(genomeFile="/Users/ananikolic/aGBM_scATAC_newcellranger/4218/outs/chrom_sizes.tsv",
                      cytobandFile="~/genomeNew/hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="~/cpg_counts_hg38.bed",
                      binSize=1e6,
                      minFrags=1e4,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)


initializeStandards(chromSizeFile="/Users/ananikolic/aGBM_scATAC_newcellranger/4218/outs/chrom_sizes.tsv",
                    cytobandFile="~/hg38_1e+06_cytoband_densities_granges.tsv",
                    cpgDataFile="~/cpg_counts_hg38.bed")

#for this tutorial we will use the sample data included in scData
#to create your own use the process_fragment_file.py script included in the package and run it on a fragments.tsv.gz file of your choosing

setOutputFile("~/","samp_dataset")

#step 1 normalize the matrix
scData<-scDataSamp
scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,dividingFactor=1,upperFilterQuantile = 0.95)

#collapse into chromosome arm level
#compiling the summary function should theoretically speed things up
library(compiler)
cutAverageC<-cmpfun(cutAverage)
summaryFunction<-cutAverageC

scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 100,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300) 
#apply additional filters
scData_collapse<-filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.01)

#show unscaled chromosome list
graphCNVDistribution(scData_collapse,outputSuffix = "test_violins")

#compute centers
median_iqr <- computeCenters(scData_collapse,summaryFunction=summaryFunction)

candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.25,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 1, subsetSize=500,fakeCellSD = 0.10, uncertaintyCutoff = 0.55,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3)
candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.0) #= 1.5)
final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "clean_cnv",sdCNV = 0.5,filterResults=TRUE,filterRange=0.8)

#identify double minutes
dmRead<-cmpfun(identifyDoubleMinutes)
dm_candidates<-dmRead(scData_k_norm,minCells=100,qualityCutoff2 = 100)
write.table(x=dm_candidates,file="samp_dm.csv",quote=FALSE,row.names = FALSE,sep=",")

#assess putative LOH regions
loh_regions<-getLOHRegions(scData_k_norm,diffThreshold = 0.20,lossCutoff = -0.75,minLength = 2e6,minSeg=2,targetFun=IQR,signalBoost=0.1,lossCutoffCells = 50,lossCutoffReads = 75,quantileLimit=0.2,cpgCutoff=0)
write.table(x=loh_regions[[1]],file=str_c("~/samp_loss.csv"),quote=FALSE,row.names = FALSE,sep=",")

library(dplyr)
loh_regions[[1]] %>% mutate_at(vars(starts_with('chr')),funs(if_else(.==2,-1,.)))