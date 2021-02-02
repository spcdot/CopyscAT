library("devtools")
#####INITIAL ONE-TIME SETUP #####
#STEP 1: replace  "~/CopyscAT" with wherever you git cloned the repo to
install("~/CopyscAT")

#alternate option: devools install-github option
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
install_github("spcdot/copyscat")

#load the package
library(CopyscAT)

#load your favourite genome
library(BSgenome.Hsapiens.UCSC.hg38)
#use it to save references - this will create a genome chrom.sizes file, a cytobands file and a CpG file
generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38",tileWidth = 1e6,outputDir = "~")

##### REGULAR WORKFLOW #####
#initialize the environment 
initialiseEnvironment(genomeFile="~/hg38_chrom_sizes.tsv",
                      cytobandFile="~/hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="~/hg38_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=1e4,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)

#for this tutorial we will use the sample data included in scData
#to create your own use the process_fragment_file.py script included in the package and run it on a fragments.tsv.gz file of your choosing

setOutputFile("~","samp_dataset")
#step 1 normalize the matrix
#USING SAMPLE DATA FROM PACKAGE
#option: if using your own file replace below with the following
#scData<-readInputTable("myInputFile.tsv")
scData<-scDataSamp
scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,dividingFactor=1,upperFilterQuantile = 0.95)

#collapse into chromosome arm level
summaryFunction<-cutAverage
#powval 0.7 - 0.75
scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 100,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300,powVal=0.73) 
#apply additional filters
scData_collapse<-filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)

#show unscaled chromosome list
graphCNVDistribution(scData_collapse,outputSuffix = "test_violinsn2")

#compute centers
median_iqr <- computeCenters(scData_collapse,summaryFunction=summaryFunction)
#identify chromosome-level amplifications
candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.25,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 0.25,bicMinimum = 0.1, subsetSize=600,fakeCellSD = 0.08, uncertaintyCutoff = 0.55,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3,IQRCutoff= 0.2,medianQuantileCutoff = 0.4)
#cleanup step
candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.5) #= 1.5)
#final results and annotation
final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "clean_cnv",sdCNV = 0.5,filterResults=TRUE,filterRange=0.8)

#identify double minutes / amplifications
#note: this is slow, and may take ~5 minutes
#option to compile this code
library(compiler)
dmRead<-cmpfun(identifyDoubleMinutes)

dm_candidates<-dmRead(scData_k_norm,minCells=100,qualityCutoff2 = 100,minThreshold = 4)
write.table(x=dm_candidates,file=str_c(scCNVCaller$locPrefix,"samp_dm.csv"),quote=FALSE,row.names = FALSE,sep=",")
?seq
paste("X",seq(from=1,to=50),"-1",sep="")
#assess putative LOH regions
t1<-rnorm(n=50,mean=300,sd = 0.1*300)
t1<-data.frame(t1)
rownames(t1)<-paste("X",seq(from=1,to=50),"-1",sep="")
as.matrix(t1)
#400
#TODO: slice the regions
sliceList<-function(listClassification,listMeans,listValues)
{
  #reassign items
  #assume there are only two
  #recompute means
  mean1<-mean(listValues[listClassification==1])
  mean2<-mean(listValues[listClassification==2])
  outliers1<-which(listValues[listClassification==1]>mean2)
  outliers2<-which(listValues[listClassification==2]<mean1)
  #reassign outliers
  listClassification[outliers1]<-2
  listClassification[outliers2]<-1
  return(listClassification)
}

loh_regions<-getLOHRegions(scData_k_norm,diffThreshold = 3,lossCutoff = -0.75,minLength = 2e6,minSeg=2,targetFun=IQR,lossCutoffCells = 200,quantileLimit=0.2,cpgCutoff=100,dummyQuantile=0.6,dummyPercentile=0.4,dummySd=0.1)
write.table(x=loh_regions[[1]],file=str_c(scCNVCaller$locPrefix,"samp_loss.csv"),quote=FALSE,row.names = FALSE,sep=",")
