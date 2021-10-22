library("devtools")
#####INITIAL ONE-TIME SETUP #####
#STEP 1: replace  "~/CopyscAT" with wherever you git cloned the repo to
install("~/CopyscAT")

#alternate option: devools install-github option
#some versions of R throw an error because of code warnings with Github install, the below line should fix this issue
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
install_github("spcdot/copyscat")

#load the package
library(CopyscAT)


#load your favourite genome
library(BSgenome.Hsapiens.UCSC.hg38)
#use it to save references - this will create a genome chrom.sizes file, a cytobands file and a CpG file
#NOTE: updated 21-Oct-2021 with a fallback to the REST API
generateReferences(BSgenome.Hsapiens.UCSC.hg38,genomeText = "hg38",tileWidth = 1e6,outputDir = "~")

##### REGULAR WORKFLOW #####
#initialize the environment (note: this needs to be done with every session in which you run Copy-scAT)
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

#SET OUTPUT DEFAULT DIRECTORY AND NAME
setOutputFile("~","samp_dataset")

#PART 1: INITIAL DATA NORMALIZATION
#step 1 normalize the matrix
#USING SAMPLE DATA FROM PACKAGE
#option: if using your own file replace below with the following
#scData<-readInputTable("myInputFile.tsv")
scData<-scDataSamp
scData_k_norm <- normalizeMatrixN(scData,logNorm = FALSE,maxZero=2000,imputeZeros = FALSE,blacklistProp = 0.8,blacklistCutoff=125,dividingFactor=1,upperFilterQuantile = 0.95)
#when using your own data, please make sure you don't have any excess / alt chromosomes

#collapse into chromosome arm level
summaryFunction<-cutAverage
scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=summaryFunction,binExpand = 1,minimumChromValue = 100,logTrans = FALSE,tssEnrich = 1,logBase=2,minCPG=300,powVal=0.73) 
#apply additional filters
scData_collapse<-filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)

#show unscaled chromosome list
graphCNVDistribution(scData_collapse,outputSuffix = "test_violinsn2")

#compute centers
median_iqr <- computeCenters(scData_collapse,summaryFunction=summaryFunction)

#PART 2: ASSESSMENT OF CHROMOSOME-LEVEL CNVs 
#OPTION 1: identify chromosome-level amplifications using all cells to generate 'normal' control
#identify chromosome-level amplifications
candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.25,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 0.25,bicMinimum = 0.1, subsetSize=600,fakeCellSD = 0.08, uncertaintyCutoff = 0.55,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3,IQRCutoff= 0.2,medianQuantileCutoff = 0.4)
#cleanup step
candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.5) #= 1.5)
#final results and annotation
final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "clean_cnv",sdCNV = 0.5,filterResults=TRUE,filterRange=0.8)



#OPTION 2: automatically identify non-neoplastic cells and use these for control
#note: still somewhat experimental - use if the copy number calls from A don't make sense (e.g. baseline appears incorrect)
#works best if tumor cellularity is <90% (very small populations of non-neoplastic cells may get missed)
#run this after filterCells
#this generates a heatmap in your working directory - check it and adjust parameters as needed
#default parameters: estimatedCellularity=0.8,nmfComponents=5,outputHeatmap=TRUE,cutHeight=0.6
#change cutHeight (height at which target dendrogram is cut) and nmfComponents as necessary to improve segeregation of clusters in your sample
nmf_results<-identifyNonNeoplastic(scData_collapse,methodHclust="ward.D")
#?identifyNonNeoplastic
write.table(x=rownames_to_column(data.frame(nmf_results$cellAssigns),var="Barcode"),file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_clusters.csv"),quote=FALSE,row.names = FALSE,sep=",")
print(paste("Normal cluster is: ",nmf_results$clusterNormal))

#ALTERNATE METHOD FOR CNV CALLING (with normal cells as background)
#you can use normals identified by NMF or known/suspected normals in the sample to refine the CNV calls (i.e. as a vector of barcode strings)
#compute central tendencies based on normal cells only
median_iqr <- computeCenters(scData_collapse %>% select(chrom,nmf_results$normalBarcodes),summaryFunction=summaryFunction)
#setting medianQuantileCutoff to -1 and feeding non-neoplastic barcodes in as normalCells can improve accuracy of CNV calls
candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,propDummy=0.25,minMix=0.01,deltaMean = 0.03,deltaBIC2 = 0.25,bicMinimum = 0.1, subsetSize=800,fakeCellSD = 0.09, uncertaintyCutoff = 0.65,summaryFunction=summaryFunction,maxClust = 4,mergeCutoff = 3,IQRCutoff = 0.25,medianQuantileCutoff = -1,normalCells=nmf_results$normalBarcodes) 
candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,medianIQR = candidate_cnvs[[3]],minDiff=1.0) #= 1.5)

#to save this data you can use annotateCNV4 as per usual
final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,outputSuffix = "clean_cnv",sdCNV = 0.6,filterResults=TRUE,filterRange=0.4)

#the other option is to use annotateCNV4B and feed in the normalBarcodes - this will set the "normal" population to 2 -- if the data is noisy it may lead to false positives so use with caution
#you may also use this version if you have a list of normal barcodes generated elsewhere
final_cnv_list<-annotateCNV4B(candidate_cnvs_clean, nmf_results$normalBarcodes, saveOutput=TRUE,outputSuffix = "clean_cnv_b2",sdCNV = 0.6,filterResults=TRUE,filterRange=0.4,minAlteredCellProp = 0.5)

#PART 2B: smoothing CNV calls with clusters
#data smoothing: can provide CNV as list or as an input file (CSV)
smoothedCNVList<-smoothClusters(scDataSampClusters,inputCNVList = final_cnv_list[[3]],percentPositive = 0.4,removeEmpty = FALSE)

#PART 3: identify double minutes / amplifications
#note: this is slow, and may take ~5 minutes
#if very large dataset, may run on subset of the data to estimate the amplifications in distinct clusters
#option to compile this code
library(compiler)
dmRead<-cmpfun(identifyDoubleMinutes)
#
#minThreshold is a time-saving option that doesn't call changepoints on any cell with a maximum Z score less than 4 - you can adjust this to adjust sensitivity of double minute calls (note - lower value = slower)
dm_candidates<-dmRead(scData_k_norm,minCells=100,qualityCutoff2 = 100,minThreshold = 4) 

write.table(x=dm_candidates,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"samp_dm.csv"),quote=FALSE,row.names = FALSE,sep=",")

#PART 4: assess putative LOH regions
#note: this is in beta, interpret results with caution
loh_regions<-getLOHRegions(scData_k_norm,diffThreshold = 3,lossCutoff = -0.75,minLength = 2e6,minSeg=2,targetFun=IQR,lossCutoffCells = 200,quantileLimit=0.2,cpgCutoff=100,dummyQuantile=0.6,dummyPercentile=0.4,dummySd=0.1)
write.table(x=loh_regions[[1]],file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"samp_loss.csv"),quote=FALSE,row.names = FALSE,sep=",")

#PART 5: quantify cycling cells
#this uses the signal from a particular chromosome (we use chromosome X as typically not altered in our samples) to identify 2N versus 4N DNA content within cells
#if there is a known alteration in X in your samples, try using a different chromosome
barcodeCycling<-estimateCellCycleFraction(scData,sampName="sample",cutoff=1000)
#take max as 
write.table(barcodeCycling[order(names(barcodeCycling))]==max(barcodeCycling),file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_cycling_cells.tsv"),sep="\t",quote=FALSE,row.names=TRUE,col.names=FALSE)

