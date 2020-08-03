#CNV caller libraries and functions
#INITIALISE ENVIRONMENT
#' @export
scCNVCaller <- new.env()

#' @export
.onLoad <- function(libname, pkgname)
{
  # ...
  message("Looking for environment")
  print(rlang::env_parent())
  print(rlang::current_env())
 # if (!exists("scCNVCaller"))
#  {
#     message("Creating new environment")
#     assign("scCNVCaller",rlang::env(rlang::global_env()),rlang::env_parent())
#  }
  # ...
}
####
#' initialiseEnvironment 
#'
#' This function prepares and preloads the genome data required for CNV analysis.
#' @param genomeFile Chromosome size file for genome of interest
#' @param cytobandFile Processed cytoband file
#' @param cpgFile Processed CpG island file
#' @param binSize Size of chromosome bins (defaults to 1e6)
#' @param minFrags Minimum # of fragments per cell (default 1e4)
#' @param cellSuffix Suffix appended to cells
#' @param lowerTrim Lower quantile trim for trimmed mean (default 0.5)
#' @param upperTrim Upper quantile trim for trimmed mean (default 0.8)
#' @keywords configuration
#' @export
initialiseEnvironment <- function(genomeFile,cytobandFile,cpgFile,binSize=1e6,minFrags=1e4,cellSuffix=c("-1"),lowerTrim=0.5,upperTrim=0.8)
{
  scCNVCaller$genomeFile = genomeFile
  #for 1e6 bins
  scCNVCaller$cytobandFile = cytobandFile
  scCNVCaller$cpgFile = cpgFile
  scCNVCaller$binSize=binSize
  scCNVCaller$minFrags=minFrags
  scCNVCaller$cellSuffix=cellSuffix
  scCNVCaller$lowerTrim=lowerTrim
  scCNVCaller$upperTrim=upperTrim
  
  #SUMMARY STAT VARIABLES
  scCNVCaller$meanReadsPerCell<-0
  scCNVCaller$startingCellCount<-0
  scCNVCaller$cellsPassingFilter<-0
  scCNVCaller$blacklistCount<-0
  scCNVCaller$finalFilterCells<-0
}
#' setOutputFile 
#'
#' This function prepares and preloads the genome data required for CNV analysis.
#' @param outputDir Output directory
#' @param outputFileSuffix Suffix to add to output files
#' @keywords configuration
#' @export
setOutputFile <- function(outputDir,outputFileSuffix)
{
  scCNVCaller$locPrefix=outputDir
  scCNVCaller$outPrefix<-outputFileSuffix
}
##INITIALIZE STANDARD FILES
#' initializeStandards 
#'
#' This function prepares and preloads the genome data required for CNV analysis.
#' @param chromSizeFile - a tsv chrom sizes file
#' @param cpgDataFile - a binned cpg density file (see documentation for how to prepare
#' @param cytobandFile - a binned cytoband file for your genome of interest (see documentation for how to prepare)
#' @keywords initialisation
#' @export
#' @examples
#' initializeStandards("hg38.chrom.sizes","hg38.cpg.1e6.bed","hg38.cytoband.1e6.bed"
initializeStandards <- function(chromSizeFile,cpgDataFile,cytobandFile)
{
  #initialize chromosome data
  scCNVCaller$chrom_sizes <- read.delim(chromSizeFile,stringsAsFactors=FALSE,header=FALSE)
  colnames(scCNVCaller$chrom_sizes) <- c("chrom","length")
  
  scCNVCaller$cpg_data<-read.table(cpgDataFile,stringsAsFactors = FALSE,header=FALSE)
  colnames(scCNVCaller$cpg_data) <- c("chrom","start","end","cpg_density")
  
  #initialize cytoband data
  scCNVCaller$cytoband_data<-read.table(cytobandFile,stringsAsFactors = FALSE,header=FALSE)
  for (n in 1:nrow(scCNVCaller$cytoband_data))
  {
    new_weight=0
    setArm=""
    split_arm<-str_split(scCNVCaller$cytoband_data$V4[n],fixed(","),simplify=TRUE)
    
    #the numbers are heterochromatin DUH
    setArm=split_arm[[length(split_arm)]]
    #print(str_sub(setArm,start=1,end=1))
    #print(cytoband_data$V4[n])
    scCNVCaller$cytoband_data$V4[n]<-str_sub(setArm,start=1,end=1)
    if (str_detect(scCNVCaller$cytoband_data$V5[n],"acen"))
    {
      scCNVCaller$cytoband_data$V4[n]<-"cen"
    }
  }
}
#END OF INITIALIZATION BLOCK
#' readInputTable
#'
#' This function prepares and preloads the genome data required for CNV analysis.
#' @param inputFile - preprocessed counts/density file
#' @param sep Separator character (default: tab)
#' @keywords initialisation
#' @keywords loading files
#' @export
#' @examples
readInputTable<-function(inputFile,sep="\t")
{
  scData<-fread(inputFile,header=TRUE,stringsAsFactors = FALSE,sep="\t",data.table = FALSE)
  #TODO: Add masking
  #scData<-scData %>% mutate(Cell_id = paste(Cell_id,scCNVCaller$cellSuffix,sep="")) 
  
  rownames(scData)<-scData$Cell_id
  scData<-scData %>% select(-Cell_id)
  return(scData)
}

#SAVE SUMMARY STATISTICS
#' saveSummaryStats
#'
#' saves summary statistics about loaded files
#' @param outputSuffix suffix to append to output file names. defaults to _stats 
#' @keywords summary
#' @export
saveSummaryStats <- function(outputSuffix="_stats")
{
  summaryStats<-list("meanReadsPerCell" = scCNVCaller$meanReadsPerCell,"startingCellCount" = scCNVCaller$startingCellCount, "cellsPassingFilter" = scCNVCaller$cellsPassingFilter, "blacklistCount" = scCNVCaller$blacklistCount, "finalFilterCells" = scCNVCaller$finalFilterCells)
  write.table(summaryStats,file=str_c(locPrefix,outPrefix,outputSuffix,".tsv"),quote=FALSE)
  return(summaryStats)
}

#' normalizeMatrix
#'
#' This function normalises the input matrix file. Deprecated. Use normalizeMatrixN instead
#' @keywords cats
#' @export
normalizeMatrix <- function(inputMatrix,logNorm=TRUE,maxZero=500,imputeZeros=FALSE,blacklistProp=0.8, priorCount=1,blacklistCutoff=0,dividingFactor=1e6)
{
  scData_t<-t(inputMatrix)
  readsPerCell<- as.data.frame(scData_t) %>% summarise_if(is.numeric,funs(sum)) %>% select(ends_with(scCNVCaller$cellSuffix)) %>% gather(Cell,Fragments,ends_with(scCNVCaller$cellSuffix))
  print(str_c("Total number of starting cells: ",round(length(readsPerCell$Fragments))," Average reads per cell: ",mean(readsPerCell$Fragments)))
  scCNVCaller$meanReadsPerCell<<-round(mean(readsPerCell$Fragments))
  #get low confidence regions
  nCells<-ncol(scData_t)-1
  scCNVCaller$startingCellCount<<-nCells
  #let's say 60% of reads are zeros is blacklist material
  zero_cols<- inputMatrix %>% summarise_if(is.numeric,funs(sum(.<blacklistCutoff))) %>% gather(Chrom,Value) 
  print(zero_cols)
  print(str_c(round(blacklistProp*nCells),"BAH"))
  #zeros greater than X % of all cells
  blacklistRegions<-zero_cols %>% filter(Value>=round(blacklistProp*nCells))
  print(str_c("Blacklist regions: ", length(blacklistRegions$Chrom)))
  print(head(blacklistRegions$Chrom))
  scCNVCaller$blacklistCount<<-length(blacklistRegions$Chrom)
  tmp1<-as.data.frame(scData_t,stringsAsFactors=FALSE)
  #remove excess zeros (poor quality cells)
  #`%notin%` <- Negate(`%in%`)
  zeros<-rownames_to_column(tmp1,var="Chrom") %>% filter(!(Chrom %in% blacklistRegions$Chrom)) %>% summarise_if(is.numeric,funs(sum(.==0))) %>% gather(Cell,Value) %>% mutate(cellPass=(Value<maxZero))
  scCNVCaller$cellsPassingFilter<<-length(which(zeros$cellPass==TRUE))
  print(str_c(cellsPassingFilter, " cells passing filter"))
  #message(blacklistRegions$Chrom)
  #hist(zeros$Value,breaks=100)
  #zeros$Cell[zeros$cellPass==TRUE]
  tmp1 <- tmp1 %>% select(zeros$Cell[zeros$cellPass==TRUE])
  raw_medians<-apply(tmp1,1,median)
  raw_medians
  #tmp1
  tmp2<-cbind.data.frame(tmp1,raw_medians,stringsAsFactors=FALSE)
  #impute zeros in cells passing filter
  if (imputeZeros==TRUE)
  {
     tmp3 <- tmp2 %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),~ if_else(. == 0,raw_medians,.)) %>% select(-raw_medians)
  }
  else
  {
    tmp3 <- tmp2
  }
  #now normalize quantiles to account for differences in coverage
  #scData_n<- normalize.quantiles(scData_t)
  scData_n<-cpm(tmp3,log =logNorm,prior.count = priorCount)
  
  rownames(scData_n)<-rownames(tmp2)
  colnames(scData_n)<-colnames(tmp3)
  scData_n
  head(scData_n)
  #var_int<-apply(scData_n[1:nrow(scData_n),], 1, var, na.rm=TRUE)
  #scData_nv<-cbind(scData_n,nvariance=var_int)
  scData_nc_split <- rownames_to_column(as.data.frame(scData_n),var = "Loc") %>% mutate(blacklist=(Loc %in% blacklistRegions$Chrom)) %>% separate(col=Loc,into=c("chrom","pos"))
  scData_k<- scData_nc_split %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./ dividingFactor))
  print(scCNVCaller$cellSuffix)
  return(scData_k)
}
#' normalizeMatrixN
#'
#' This function normalizes the input matrix. Do this first.
#' @param inputMatrix - processed fragment table from your input
#' @param logNorm - whether to log normalize input file. default: FALSE 
#' @param maxZero - maximum # of zeros per cell after removal of uninformative regions
#' @param imputeZeros - default to FALSE; will impute cells with medians - don't use it
#' @param blacklistProp - proportion of cells that have zeros in an interval to consider that interval blacklisted/uninformative - lower = more stringency
#' @param priorCount - use if log normalizing data - priorCount to add to fields before normalisation
#' @param blacklistCutoff - minimum signal in a bin to count it as "empty" for the blacklist cutoff
#' @param dividingFactor - deprecated
#' @export
normalizeMatrixN <- function(inputMatrix,logNorm=FALSE,maxZero=500,imputeZeros=FALSE,blacklistProp=0.8, priorCount=1,blacklistCutoff=0,dividingFactor=1e6)
{
  sc_t<-data.table(t(inputMatrix))
  #sc_t
  cellReads<-transpose(sc_t[,lapply(.SD,sum)],keep.names="Cell")
  readsCells=cellReads[,mean(V1)]
  nCells<-nrow(cellReads)
  print(str_c("Total number of starting cells: ",nCells," Average reads per cell: ",mean(readsCells)))
  scCNVCaller$meanReadsPerCell<-readsCells
  blacklistPropCutoff=blacklistProp*nCells
  #good
  #find bad columns
  sc_lines<-data.table(scData)
 # head(sc_lines)
  #blacklistCutoff = 500
  sc_lines<-sc_lines[,lapply(.SD,function(x) x<blacklistCutoff)][,lapply(.SD,sum)]
  sc_pos<-transpose(sc_lines,keep.names = "Pos")
  blacklistRegions<-sc_pos[which(sc_pos[,V1>=blacklistPropCutoff]),]$Pos
  print(length(blacklistRegions))
  print(nrow(sc_pos))
  if (length(blacklistRegions)==nrow(sc_pos))
  {
    print("Error: no regions meet cutoff criteria")
    return(NULL)
  }
  sc_t2<-sc_t[which(sc_pos[,V1<blacklistPropCutoff]),lapply(.SD,function(x) as.numeric(x<blacklistCutoff))][,lapply(.SD,sum)]
  print(sc_t[which(sc_pos[,V1<blacklistPropCutoff]),lapply(.SD,function(x) as.numeric(x<blacklistCutoff))][1:10,1:10])
  #print(sc_t2[,1:100])
  
  sc_zeros<-transpose(sc_t2,keep.names="Cells")
  #print(sc_zeros)
  #zero_cutoff=zero_cutoff
  zero_list<-sc_zeros[which(sc_zeros[,V1<maxZero])]$Cells
  zero_list
  #print(zero_list)
  #zeros<-rownames_to_column(tmp1,var="Chrom") %>% filter(!(Chrom %in% blacklistRegions$Chrom)) %>% summarise_if(is.numeric,funs(sum(.==0))) %>% gather(Cell,Value) %>% mutate(cellPass=(Value<maxZero))
  scCNVCaller$cellsPassingFilter<-length(zero_list)
  print(str_c(scCNVCaller$cellsPassingFilter, " cells passing filter"))
  #get low confidence regions
  #message(blacklistRegions$Chrom)
  #hist(zeros$Value,breaks=100)
  #zeros$Cell[zeros$cellPass==TRUE]
  tmp1 <- as.data.frame(sc_t,stringsAsFactors=FALSE) %>% select(zero_list)
  raw_medians<-t(transpose(sc_t[,..zero_list])[,lapply(.SD,median)])
  #tmp1
  tmp1<-cbind.data.frame(tmp1,raw_medians,stringsAsFactors=FALSE)
  #impute zeros in cells passing filter
  if (imputeZeros==TRUE)
  {
    tmp3 <- tmp1 %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),~ if_else(. == 0,raw_medians,.)) %>% select(-raw_medians)
  }
  else
  {
    tmp3 <- tmp1  # %>% select(-raw_medians)
  }
  #now normalize quantiles to account for differences in coverage
  #scData_n<- normalize.quantiles(scData_t)
  scData_n<-cpm(tmp3,log = logNorm,prior.count = priorCount)
  
  rownames(scData_n)<-colnames(inputMatrix)
  colnames(scData_n)<-colnames(tmp3)
  #scData_n
  head(scData_n)
  #var_int<-apply(scData_n[1:nrow(scData_n),], 1, var, na.rm=TRUE)
  #scData_nv<-cbind(scData_n,nvariance=var_int)
  scData_nc_split <- rownames_to_column(as.data.frame(scData_n),var = "Loc") %>% mutate(blacklist=(Loc %in% blacklistRegions)) %>% separate(col=Loc,into=c("chrom","pos"))
  scData_k<- scData_nc_split %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./ dividingFactor))
  return(scData_k)
}

#CALCULATE DMs
#' getDoubleMinutes
#'
#' This function allows you to express your love of cats.
#' @keywords double minutes
#' @export
getDoubleMinutes = function(inputMatrix,targetCell,doPlot=FALSE,penalty_type="SIC",doLosses=FALSE,peakCutoff=3,lossCutoff=-1)
{
  medianDensity<-median(inputMatrix[,targetCell])
  if (medianDensity>0){
    tmp<-inputMatrix[,targetCell]
    #tmp[tmp>4e-6]<-medianDensity
    #consider using IQR Here
    tmp<-as.numeric(scale(tmp,center=median(tmp),scale=IQR(tmp)),stringsAsFactors=FALSE)
    #adjusted minseglength 
    cm<-cpt.meanvar(data=tmp,test.stat="Normal", penalty=penalty_type,method = "PELT")
    #,minseglen = 3)  
    #TODO: need to document double minutes in some way (e.g. presence of each per cell) - this works
    #output cptlist
    if (doPlot==TRUE)
    {
      print(plot(cm,xlab="Chromosome bin",ylim=c(0,50),ylab="Chromosome Z-score"))
    }
    cptlist<-t(rbind(cm@param.est$mean,cm@cpts))
    colnames(cptlist)<-c("Mean","Point")
    cptlist<-as_tibble(cptlist) %>% mutate(Diff = Point - lag(Point))
    #set up list of double minutes; cutoff 0.002
    d_minutes<-cptlist %>% filter(Mean>peakCutoff)
    d_minutes
    d_losses<-cptlist %>% filter(Mean<lossCutoff)
    #message(nrow(d_minutes))
    if (!doLosses)
    {
      if (nrow(d_minutes)>0)
      {
        #message(d_minutes)
        #message("test")
        #for each row - to extend for more than one spike
        d_min_coords<-vector(mode="character",length=nrow(d_minutes))
        
        for (k in 1:nrow(d_minutes))
        {
          if (is.na(d_minutes$Diff[k]))
          {
            #assign if first value to a start of 1
            d_minutes$Diff[k]=d_minutes$Point[k]-1
          }
          #just do column file here
          #message(d_minutes$Point[k])
          #message("b")
          #message(d_minutes$Point[k] - d_minutes$Diff[k])
          #message(d_minutes$Point[k])
          coords_temp<-inputMatrix[(d_minutes$Point[k] - d_minutes$Diff[k]):d_minutes$Point[k],]$loc1
          #message(coords_temp)
          #message("c")
          tmp_coords<-sprintf(fmt="%s.%s",coords_temp[1],coords_temp[length(coords_temp)-1])
          d_min_coords[k]=tmp_coords
        }
        return(d_min_coords)
      }
      else {
        return(NULL)
      }
    }
    else {
      #look for losses
      if (nrow(d_losses)>0)
      {
        #message(d_minutes)
        #message("test")
        #for each row - to extend for more than one spike
        d_loss_coords<-vector(mode="character",length=nrow(d_losses))
        
        for (k in 1:nrow(d_losses))
        {
          if (is.na(d_losses$Diff[k]))
          {
            #assign if first value to a start of 1
            d_losses$Diff[k]=d_losses$Point[k]-1
          }
          #just do column file here
          #message(d_minutes$Point[k])
          #message("b")
          #message(d_minutes$Point[k] - d_minutes$Diff[k])
          #message(d_minutes$Point[k])
          coords_temp<-inputMatrix[(d_losses$Point[k] - d_losses$Diff[k]):d_losses$Point[k],]$loc1
          #message(coords_temp)
          #message("c")
          tmp_coords<-sprintf(fmt="%s.%s",coords_temp[1],coords_temp[length(coords_temp)-1])
          d_loss_coords[k]=tmp_coords
        }
        #return list of double minutes for target cell
        #d_min_coords<-inputMatrix[(d_minutes$Point - d_minutes$Diff):d_minutes$Point,]$loc1 
        #d_min_coords
        
        return(d_loss_coords)
      } 
    }
  }
  else {
    return(NULL)
  }
}
#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @examples
#' cat_function()
testOverlap <- function(interval.1,interval.2){
  isOverlap=FALSE
  leftEnd="0"
  rightEnd="0"
  interval.1.sep<-sapply(str_split(interval.1,fixed(".")),str_split,fixed("_"))
  interval.2.sep<-sapply(str_split(interval.2,fixed(".")),str_split,fixed("_"))
  
  if (as.numeric(interval.1.sep[[1]][2])<=as.numeric(interval.2.sep[[1]][2]))
  {
    leftEnd=interval.1.sep[[1]][2]
    #next check other end; if bigger then overlap is true
    if (as.numeric(interval.1.sep[[2]][2])>=as.numeric(interval.2.sep[[1]][2]))
    {
      isOverlap=TRUE
      rightEnd=interval.1.sep[[2]][2]
      if (as.numeric(interval.2.sep[[2]][2])>=as.numeric(interval.1.sep[[2]][2]))
      {
        rightEnd=interval.2.sep[[2]][2]
      }
    }
  } else {
    #right interval starts first  
    leftEnd=interval.2.sep[[1]][2]
    #next check other end; if bigger then overlap is true
    if (as.numeric(interval.2.sep[[2]][2])>=as.numeric(interval.1.sep[[1]][2]))
    {
      isOverlap=TRUE
      rightEnd=interval.2.sep[[2]][2]
      if (as.numeric(interval.2.sep[[2]][2])<=as.numeric(interval.1.sep[[2]][2]))
      {
        rightEnd=interval.1.sep[[2]][2]
      }
      
    }
  }
  return(c(leftEnd,rightEnd,isOverlap))
}
#head(scData_k)

#PURPOSE: identify double minutes on all chromosomes using changepoint analysis
#' identifyDoubleMinutes
#'
#' Identify double minutes on all chromosomes using changepoint analysis
#' @keywords double minutes
#' @keywords amplification 
#' @keywords ecDNA
#' @export
identifyDoubleMinutes<-function(inputMatrix,minCells=20,qualityCutoff2=40,logTrans=FALSE,cpgTransform=FALSE,peakCutoff=3,lossCutoff=-1,doLosses=FALSE,doPlots=TRUE,imageNumber=250)
{
  dm_per_cell<-data.frame(cellName=colnames(inputMatrix  %>% select(ends_with(scCNVCaller$cellSuffix))),stringsAsFactors=FALSE)
#  print(head(dm_per_cell))
  dm_per_cell[,seq(from=2,to=5)]<-FALSE
  
  #initialize repeat index
  firstRepeatIndex=2
  #MINIMUM # of cells with alteration - 50 is defualt
  qualityCutoff=minCells
  #test - can remove centromeres
  #cytoband_data$V4
  scData_k_cpg<-inputMatrix %>% mutate(cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(arm=scCNVCaller$cytoband_data$V4) %>% filter(arm!="cen", blacklist==0) %>% select(-arm,-blacklist) #%>%  #cpg+
  if (cpgTransform==TRUE)
  {
    if (logTrans==FALSE)
    {
      scData_k_cpg <- scData_k_cpg %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(cpg+1)))
    }
    else
    {
      scData_k_cpg <- scData_k_cpg %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(log(cpg,base=2)+1)))
    }
  }
  for (m in 1:(nrow(scCNVCaller$chrom_sizes)-3))
  {
  
    #loop through chromosomes
    currentChrom = scCNVCaller$chrom_sizes$chrom[m]
    message(str_c("processing ",currentChrom))
   # print(nrow(scData_k_cpg))
    #scData_1h <- scData_k_cpg %>% filter(chrom==currentChrom ) %>% mutate(loc1=paste(chrom,pos,sep="_")) %>% select(-chrom,-pos)
    scData_1h <- scData_k_cpg %>% filter(chrom==currentChrom ) %>% mutate(loc1=paste(chrom,pos,sep="_")) %>% select(-chrom,-pos,-cpg)
  #  print(scData_1h)
    #seleect peak cutoff
    #unlistValues <- unlist((scData_1h %>% select(-loc1)))
    #peakCutoff = quantile(scale(unlistValues,center=median(unlistValues),scale=IQR(unlistValues)),.99)
    
    #scData_1h[,2:ncol(scData_1h)]<-scale(scData_1h[,2:ncol(scData_1h)],center=TRUE,scale=TRUE)
    dm_list=data.frame(dm="null",stringsAsFactors=FALSE)
    if (m>1)
    {
      rm(dm_per_cell_clean)
      rm(recurrentEvents)
    }
    #initialize temporary matrix
    dm_per_cell_temp<-as.data.frame(matrix(nrow = nrow(dm_per_cell),ncol=5))
    colnames(dm_per_cell_temp)<-c("cellName","V2","V3","V4","V5")
    dm_per_cell_temp[,1]<-dm_per_cell$cellName
    #blank out temporary file
    #need to be smarter with this bit here but I can't be bothered RN
    dm_per_cell_temp[,seq(from=2,to=5)]<-FALSE
    
    #START TESTING HERE
    #temp_dm<-getDoubleMinutes(scData_1h,17,doPlot=TRUE)
    #dm_index<-which(dm_list$dm==temp_dm_boundaries)
    #dm_index
    plottableCell=FALSE
    plotEachTime=imageNumber
    for (cellNum in 1:(nrow(dm_per_cell)-1))
    {
    #  print(cellNum)
      # medianDensity<-median(scData_1h[,5])
      # if (medianDensity>0){
      #   tmp<-scData_1h[,5]
      #   #replace telomere values
      #   tmp[tmp>4e-6]<-medianDensity
      #   plot(tmp)
      #   #consider using IQR Here
      #   tmp<-as.numeric(scale(tmp,center=median(tmp),scale=IQR(tmp)),stringsAsFactors=FALSE)
      if ((cellNum %% plotEachTime)==0)
      {
        plottableCell=TRUE
      #  print(cellNum)
      #  print(plottableCell)
      }
      if (doPlots & plottableCell)
      {
        #print(plottableCell)
        pdf(str_c("~/cell",cellNum,"_",currentChrom,".pdf"),width=6,height=4)
        temp_dm<-getDoubleMinutes(scData_1h,cellNum,doLosses=doLosses,peakCutoff = peakCutoff,lossCutoff = lossCutoff,doPlot=TRUE)
        dev.off()
        plottableCell=FALSE
      }
      else
      {
        temp_dm<-getDoubleMinutes(scData_1h,cellNum,doLosses=doLosses,peakCutoff = peakCutoff,lossCutoff = lossCutoff,doPlot=FALSE)
      }
      
      if (!is.null(temp_dm))
      {
        #iterate through each DM identified
        for (j in 1:length(temp_dm))
        {
          #LENGTH CHECK on temp_dm
          #get boundaries - this will be our identifier for this location
          #OBSOLETE: temp_dm_boundaries<-sprintf(fmt="%s.%s",temp_dm[j,1],temp_dm[j,ncol(temp_dm)])
          #this works to ID if this vector already exists
          #,arr.ind=TRUE
          dm_index<-which(dm_list==temp_dm[j])
          #message(dm_index)
          #message("identify index")
          if (length(dm_index)==0)
          {
            #message("adding to list")
            dm_list<-rbind.data.frame(dm_list,data.frame(dm=temp_dm[j],stringsAsFactors = FALSE),stringsAsFactors=FALSE)
            dm_per_cell_temp[cellNum,nrow(dm_list)]=TRUE
            #TODO: get way to name this after the interval (name column name)
            # dm_per_cell<-cbind.data.fraame(dm_per_cell,data.frame(temp_dm_boundaries=rep(FALSE,times=nrow(dm_per_cell))))
          } else {
            
            # dm_index2<-which(dm_list==temp_dm_boundaries,arr.ind=TRUE)
            #because cellName is a column, and first value of dm is always blank, so don't add one
            dm_per_cell_temp[cellNum,dm_index]=TRUE
          }
        }
      }
    }
    names(dm_per_cell_temp)[2:nrow(dm_list)]<-dm_list$dm[2:nrow(dm_list)]
    #dm_per_cell
    #dm_per_cell %>% rename_at(2:(nrow(dm_list)-1),~ dm_list$dm)
    #test for overlaps - start index 2
    #start one smaller
    
    dm_list_test2=dm_list %>% mutate(newIndex=1)
    dm_list_test2
    chrom_name=currentChrom
    #now loop through all intervals - if detected greater than one
    if (nrow(dm_list_test2)>1)
    {
      for (j1 in 2:(nrow(dm_list_test2)))
      {
        if(dm_list_test2$newIndex[j1]==1)
        {
          dm_list_test2$newIndex[j1]=j1
        }
        #(j1+1)
        for (j2 in 2:(nrow(dm_list_test2)))
        {
          if (dm_list_test2$newIndex[j2]!=j1)
          {
            #message(j2)
            interval.1<-dm_list_test2$dm[j1]
            interval.2<-dm_list_test2$dm[j2]
            overlapResult<-testOverlap(interval.1,interval.2)
            #message(overlapResult[3])
            if (overlapResult[3]=="TRUE")
            {
              #update table
              dm_list_test2$dm[j1]=str_c(chrom_name,"_",overlapResult[1],".",chrom_name,"_",overlapResult[2])
              #message(dm_list_test2$dm[j2])
              dm_list_test2$dm[j2]=str_c(chrom_name,"_",overlapResult[1],".",chrom_name,"_",overlapResult[2])
              #message(dm_list_test2$dm[j2])
              dm_list_test2$newIndex[j2]=dm_list_test2$newIndex[j1]
            }
          }
        }
      }
      #NOTE: still some minor kinks in the overlap function; it works but the intervals aren't coming up perfect
      dm_list_test2
      #TODO: fudge indices to empty blanks
      tmp_fact<-factor(dm_list_test2$newIndex)
      levels(tmp_fact)<-seq(1,length(levels(tmp_fact)))
      dm_list_test2$newIndex<-as.numeric(levels(tmp_fact)[tmp_fact])
      dm_list_test2
      #option: can clean up dataframe here and replace NA with FALSE, and then proceed instead of making a giant one up front
      dm_per_cell_temp[is.na(dm_per_cell_temp)]<-FALSE
      dm_per_cell_temp
      #OK that worked, now to amalgamate cell-wide data
      for (j3 in 2:nrow(dm_list_test2))
      {
        #combine and rename columns appropriately
        dm_per_cell_temp[,dm_list_test2$newIndex[j3]]<-(dm_per_cell_temp[,j3] | dm_per_cell_temp[,dm_list_test2$newIndex[j3]])
      }
      dm_per_cell_temp
      #now rename colnames
      for (j3 in 2:nrow(dm_list_test2))
      {
        colnames(dm_per_cell_temp)[dm_list_test2$newIndex[j3]]<-dm_list_test2$dm[j3]
      }
      dm_per_cell_temp<-dm_per_cell_temp[,1:(max(dm_list_test2$newIndex))]
      
      #now summarise and count each; remove any with 5 or fewer events
      
      recurrentEvents <- dm_per_cell_temp  %>% summarize_if(is.logical,sum) %>% select_if(function(x) any(x>qualityCutoff))
      #select only these columns that meet cutoff
      #check that recurrent events have been identified
      if (dim(recurrentEvents)[2]!=0)
      {
        message("copying recurrent events")
        dm_per_cell_clean <- dm_per_cell_temp %>% select(cellName,colnames(recurrentEvents))
        #copy over columns and column names and increment column counter - ncol-2
        dm_per_cell[,firstRepeatIndex:(firstRepeatIndex+ncol(dm_per_cell_clean)-2)]=dm_per_cell_clean[,2:(ncol(dm_per_cell_clean))]
        #check indexes in following line
        colnames(dm_per_cell)[firstRepeatIndex:(firstRepeatIndex+ncol(dm_per_cell_clean)-2)]<-colnames(dm_per_cell_clean[2:(ncol(dm_per_cell_clean))])
        firstRepeatIndex=firstRepeatIndex + ncol(dm_per_cell_clean)-1
      }
    }
    #beautiful
  }
  #now clean this
  dm_per_cell_clean<-dm_per_cell[,1:(firstRepeatIndex-1)]
  dm_per_cell_clean %>% summarize_if(is.logical,sum)
  
  #can cutoff here too (def 50-75)
  secondThreshold=qualityCutoff2
  keepAlterations <- dm_per_cell_clean %>% summarize_if(is.logical,sum) %>% gather(Chrom,Value) %>% filter(Value>secondThreshold)
  dm_per_cell_clean <- dm_per_cell_clean %>% select(cellName, keepAlterations$Chrom)
  return(dm_per_cell_clean)
}
#LOG TRANSFORM



#dynamic cutoff isn't bad

#CAN DO SECOND ITERATION THROUGH THESE REGIONS and sum signal in these regions then assign (for things that don't quite meet cutoff)
#' callDMGaussian
#'
#' Calls DMs identified using Gaussian decomposition of the regions of interest. Useful but not as sensitive
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @examples
callDMGaussian<-function(inputMatrix,dm_results)
{
  #TODO: 
  t1<-as.data.frame(x=colnames(dm_results)[2:ncol(dm_results)])
  colnames(t1)="Alts"
  t1 <- t1 %>% separate(Alts,c(".c"),into=c("loc1","loc2")) %>% separate(loc1,c("_"),into=c("chrom","start")) %>% mutate(loc2=str_replace(loc2,"hr[0-9]+_",""))
  
  #back to this 
  tmp_scaled<-inputMatrix 
  tmp_scaled[,3:ncol(tmp_scaled)]=scale(tmp_scaled[,3:ncol(tmp_scaled)],center=TRUE,scale=TRUE)
  dm_per_cell_vals<-data.frame(cellName=colnames(inputMatrix  %>% select(-chrom,-pos)),stringsAsFactors=FALSE)
  alteration_list=c()
  for (i in 1:nrow(t1))
  {
    #message(i)
    #process 
    #message(t1$loc2[i])
    posList=seq(from=as.numeric(t1$start[i]),to=as.numeric(t1$loc2[i]),by = 1e6)
    posList<-sprintf(fmt="%d",posList)
    posList
    t2<-tmp_scaled %>% filter(chrom==t1$chrom[i],pos %in% posList) %>% summarise_if(is.numeric,sum) # %>% select(-cpg)
    t2d<-t2 %>% gather(Cell,Value)
    fit1<-Mclust(t2d$Value,G=1:3,modelNames=c("V"),initialization = list(noise=TRUE))
    if (fit1$G>2)
    {
      #collapse clusters
      fit1$classification[fit1$data<fit1$parameters$mean[2] & fit1$classification==3]<-2
      fit1$classification[fit1$data>fit1$parameters$mean[3] & fit1$classification==2]<-3
      
      fit1$classification[fit1$classification==2]<-1
      fit1$classification[fit1$classification==3]<-2
      alteration_list<-cbind(alteration_list,colnames(dm_results[i+1]))
      message(fit1$parameters$mean)
      dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification)
      #message(fit1$parameters)
    }
    if (fit1$G==2)
    {
      #check for and collapse cluster assignments
      #cells in 2 with val < 1 and vice versa
      fit1$classification[fit1$data<fit1$parameters$mean[1] & fit1$classification==2]<-1
      fit1$classification[fit1$data>fit1$parameters$mean[2] & fit1$classification==1]<-2
      alteration_list<-cbind(alteration_list,colnames(dm_results[i+1]))
      message(fit1$parameters$mean)
      dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification)
    }
    #plot(fit1)
    #hist(t2d$Value,breaks=50)
  }
  
  colnames(dm_per_cell_vals)[2:ncol(dm_per_cell_vals)]<-alteration_list
  dm_per_cell_vals[,2:ncol(dm_per_cell_vals)]=apply(dm_per_cell_vals[2:ncol(dm_per_cell_vals)],2,factor,labels = c("Unclassified","Normal","Amplified"))
  return(dm_per_cell_vals)
  
}

#write.table(dm_per_cell_vals,file="~/4349_amp_v_bin_nolog.csv",quote=FALSE,row.names = FALSE,sep=",")



#PART 2: CNV CALLING
#PARAMETERS

#options: cutAverage, mean, median, etc

#functions
#' cutAverage
#'
#' Averages regions of vector between prespecified quantiles
#' @param inputVector - vector of numbers
#' @keywords average`
#' @export
#' @examples
cutAverage <- function(inputVector)
{
  #remove top and bottom 5%
  #other option is to use denstiy
  #tmp_val_dens <- density(inputVector)
  #return(tmp_val_dens$x[which.max(tmp_val_dens$y)])
  
  #OPTION 5, use average of quantiles; 20 and 80
  quants_tmp <- quantile(inputVector, c(scCNVCaller$lowerTrim,scCNVCaller$upperTrim),na.rm=TRUE,type=8)
  if (!is.na(quants_tmp[1]))
  {
  #or find median of data within
  # message(str_c(quants_tmp[1],quants_tmp[2],sep=","))
  boolVal<-inputVector>=quants_tmp[1] & inputVector<=quants_tmp[2]
  #message(boolVal)
  #print(quants_tmp)
  trimInput <- inputVector[boolVal]
  #message(trimInput,inputVector)
  if (quants_tmp[2]!=0)
  {
    return(mean(trimInput))
  }
  else
  {
    return(mean(inputVector[boolVal]))
  }
  }
  else
  {
    return(NA)
  }
}
#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @examples
#' cat_function()
scaleMatrixBins <- function(inputMatrix, binSizeRatio,inColumn)
{
  newMatrix<-inputMatrix %>% mutate(pos_b=floor(as.numeric(!!rlang::sym(inColumn))/(scCNVCaller$binSize*binSizeRatio))*(scCNVCaller$binSize*binSizeRatio))
  newMatrix<-newMatrix %>% select(-inColumn) %>% group_by(pos_b,chrom) %>% summarise_if(is.numeric,list(sum)) %>% arrange(chrom,pos_b)
  
  return(newMatrix)
}
#move this elsewhere
isMale=FALSE

#old version of code

#' collapseChrom
#' deprecated
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords deprecated
#' @examples
collapseChrom<-function(inputMatrix,minimumSegments=40,summaryFunction=cutAverage,logTrans=TRUE,binExpand=1,minimumChromValue=2)
{
  #remove centromeric bins
  #these are ordered in the same fashion so should be accurate
  #divide signal by cpg density + 1 (To avoid dividing by zero)
  #remove blacklist regions
  scData_k_norm <- inputMatrix %>% mutate(chromArm=scCNVCaller$cytoband_data$V4,cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(chrom=str_c(chrom,chromArm))
  if (binExpand>1)
  {
    #collapse bins
    scData_k_norm <- scaleMatrixBins(scData_k_norm,binSizeRatio = binExpand,"pos")
    #rename column
    scData_k_norm <- scData_k_norm %>% rename(pos = pos_b)
  }
  #gender determination
  xy_signal <-scData_k_norm %>% filter(chrom %in% c("chrXq","chrYq")) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
 rowwise() %>% mutate(totalcpg=summaryFunction(cpg)) %>% select(-cpg)
  diffCutoff=5e-7
  if (logTrans) {
    xy_signal <- xy_signal %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(totalcpg,base=2)))) 
  }
  else {
    xy_signal <- xy_signal %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+totalcpg))) 
    diffCutoff=0.25
  }
  XYmedians <- xy_signal %>% gather(Cell,Density,ends_with(scCNVCaller$cellSuffix)) %>% group_by(chrom) %>% summarise(med=median(Density))
  diffXY<-XYmedians$med[1]-XYmedians$med[2]
  if (abs(diffXY)<diffCutoff)
  {
    isMale=TRUE
  }
  scData_k_norm <- scData_k_norm %>% mutate(isCentromere=str_detect(chrom,"cen")) %>% filter(isCentromere==FALSE,blacklist==0) %>% select(-blacklist)
  
  message(nrow(scData_k_norm))
  #now append 
  #summaryFunction
  total_cpg<-scData_k_norm %>% select(chrom,cpg) %>%
    group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction))
  #print(tail(total_cpg))
  
  #recompute medians after normalization
  median_chromosome_density<-scData_k_norm %>%
    group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) 
  #can we normalize signal to X chromosome signal?
  
  median_chrom_signal<-median_chromosome_density %>% gather(Cell,Value,2:ncol(median_chromosome_density)) %>% select(-Cell) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction))
  print(median_chrom_signal)
  median_chrom_signal<-median_chrom_signal %>% filter(!is.na(Value))
  print(median_chrom_signal)
  IQR_chrom_signal<-median_chromosome_density %>% filter(chrom %in% median_chrom_signal$chrom) %>% gather(Cell,Value,2:ncol(median_chromosome_density)) %>% select(-Cell) %>% group_by(chrom) %>% summarise_if(is.numeric,list(IQR))
  print(IQR_chrom_signal)
  
  #TODO: cpg of 200 and medianNorm cutoffs only work for 1e6 bins; need to change for others
  if (logTrans==TRUE)
  {
    #5.1e-6 old 
    #5.1e-06, totalcpg 130 for 200k
    scData_prechrom<-scData_k_norm %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    mutate(totalcpg=total_cpg$cpg) %>% select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(totalcpg,base=2)))) %>% filter(chrom %in% median_chrom_signal$chrom) %>% 
    mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>minimumChromValue,totalcpg>(200 * binExpand)) 
    #+log(totalcpg,base=2)
    #was 200
  }
  else
  {
    scData_prechrom<-scData_k_norm %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    mutate(totalcpg=total_cpg$cpg) %>% select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+totalcpg))) %>% filter(chrom %in% median_chrom_signal$chrom) %>%
    mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>0,totalcpg>(200* binExpand)) 
  }
  
  #return the appropriate matrices
  return(scData_prechrom)
}
#' collapseChrom3N
#'
#' This function summarises signal across chromosomes and normalises using the CpG content
#' @param inputMatrix: normalised matrix
#' @param summaryFunction: function to use to summarise signal
#' @param logTrans: whether matrix is log-transformed
#' @param binExpand: whether to merge adjacent bins, if so how many (merges if > 1)
#' @param minimumChromValue: cutoff to keep a chromosome 
#' @keywords CNV
#' @export
collapseChrom3N<-function(inputMatrix,minimumSegments=40,summaryFunction=cutAverage,logTrans=FALSE,binExpand=1,minimumChromValue=2)
{
  #remove centromeric bins
  #these are ordered in the same fashion so should be accurate
  #divide signal by cpg density + 1 (To avoid dividing by zero)
  #remove blacklist regions
  # %>% filter(cpg!=0)
  scData_k_norm <- inputMatrix %>% mutate(chromArm=scCNVCaller$cytoband_data$V4,cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(chrom=str_c(chrom,chromArm)) %>% filter(cpg!=0)
  if (binExpand>1)
  {
    #collapse bins
    scData_k_norm <- scaleMatrixBins(scData_k_norm,binSizeRatio = binExpand,"pos")
    #rename column
    scData_k_norm <- scData_k_norm %>% rename(pos = pos_b)
  }
  # %>% filter(cpg!=0)
  #gender determination
  sckn<-data.table(scData_k_norm)
  total_cpg<-sckn[,summaryFunction(cpg),by="chrom"]
  sckn[,c("blacklist","raw_medians","chromArm","cpg","pos"):=NULL]
  tail(colnames(sckn))
  setkey(sckn,chrom)
  sckn<-sckn[c("chrXq","chrYq"),lapply(.SD,quantile,probs=0.8),by="chrom"]
  nrow(sckn)
  total_cpg[chrom %in% c("chrXq","chrYq")]
  sckn[,cpg:=total_cpg[chrom %in% c("chrXq","chrYq")]$V1]
  #sckn[,print(.SD)]
  xy_signal<-transpose(sckn[,lapply(.SD,"/",1+cpg),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
  
  #t(sckn)[2:ncol(sckn),]
  #xy_signal <- as_tibble(apply(t(sckn)[2:ncol(sckn),],2,as.numeric))
  #xy_signal
  #xy_signal <-scData_k_norm %>% filter(chrom %in% c("chrXq","chrYq")) %>% group_by(chrom) %>% summarise_if(is.numeric,list(quantile),probs=0.8) 
  sexCutoff=5e-7
  logTrans=FALSE
  if (logTrans)
  {
    xy_signal<-transpose(sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
    
    #xy_signal <- xy_signal %>% rowwise() %>% mutate(totalcpg=summaryFunction(cpg)) %>% select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(totalcpg,base=2)))) 
  } else
  {
    xy_signal<-transpose(sckn[,lapply(.SD,"/",1+cpg),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
    
    #xy_signal <- xy_signal %>% rowwise() %>% mutate(totalcpg=summaryFunction(cpg)) %>% select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+totalcpg))) 
    sexCutoff=0.25
  }
  print(xy_signal)
  diffXY<-xy_signal$chrXq-xy_signal$chrYq
  print(abs(diffXY))
  
  if (abs(diffXY)<sexCutoff)
  {
    isMale=TRUE
  }
  #top part works awesome
  
  #now to convert rest
  sckn<-data.table(scData_k_norm)
  print(nrow(sckn))
  sckn <- sckn[chromArm!="cen"][blacklist==FALSE]
  print(nrow(sckn))
  #total_cpg <- sckn[,summaryFunction(cpg),by="chrom"]
  #total_cpg
  #setnames(total_cpg,"V1","totalcpg")
  #remove raw_medians,blacklist,chromArm,cpg
  sckn<-sckn[,c("blacklist","raw_medians","chromArm","pos"):=NULL]
  #mm I see the issue - the cutAverage upper quantile busts / overcorrects the cpg - use a median instead
  cpg_density<-sckn[,cpg]
  median_chromosome_density<-sckn[1:(nrow(sckn)),lapply(.SD,summaryFunction),by="chrom"]
  med_temp<-transpose(median_chromosome_density[,c("cpg"):=NULL],make.names="chrom")
  #head(med_temp)
  transpose(med_temp[,lapply(.SD,summaryFunction)],keep.names="chrom")
  median_chrom_signal<-transpose(med_temp[,lapply(.SD,summaryFunction)],keep.names="chrom")
  print(median_chrom_signal)
  
  
  #median_chromosome_density[,1:5]
  #]
  #median_chromosome_density<-scData_k_norm %>%
  #  group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) 
  message(nrow(scData_k_norm))
  #now append 
  #summaryFunction
  
  
  
  #TODO: cpg of 200 and medianNorm cutoffs only work for 1e6 bins; need to change for others
  if (logTrans==TRUE)
  {
    #5.1e-6 old 
    #5.1e-06, totalcpg 130 for 200k
    #THIS ONE
    #transpose(sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL]
    sckn<-sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL][,lapply(.SD,sum),by="chrom"]
    #scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(cpg,base=2)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    #  select(-cpg) %>% filter(chrom %in% median_chrom_signal$chrom) %>% 
    #  mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>minimumChromValue) 
    #+log(totalcpg,base=2)
    #was 200
  }
  else
  {
    #THIS ONE
    sckn<-sckn[,lapply(.SD,"/",log(1+cpg)),by="chrom"][,cpg:=NULL][,lapply(.SD,sum),by="chrom"]
    #scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+sqrt(cpg)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    #  select(-cpg) %>% filter(chrom %in% median_chrom_signal$chrom) %>%
    #  mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>minimumChromValue) 
  }
  
  sckn<-sckn[,medianNorm:=median_chrom_signal$V1]
  sckn<-sckn[which(sckn[,medianNorm>minimumChromValue])]
  #return the appropriate matrices
  return(as.data.frame(sckn))
}
#' collapseChromN
#'
#' This function collapses bins into chromosomes.
#' @param inputMatrix normalised input matrix
#' @param summaryFunction function to use to collapse chromoseome 
#' @param logTrans is the data log-transformed - TRUE or FALSE (default: FALSE)
#' @param binExpand - are we collapsing adjacent bins? >1 = true; integer
#' @param minimumChromValue - minimum cutoff to keep a chromosome arm as informative
#' @keywords CNV
#' @export
collapseChromN<-function(inputMatrix,summaryFunction=cutAverage,logTrans=FALSE,binExpand=1,minimumChromValue=2)
{
  #remove centromeric bins
  #these are ordered in the same fashion so should be accurate
  #divide signal by cpg density + 1 (To avoid dividing by zero)
  #remove blacklist regions
  # %>% filter(cpg!=0)
  scData_k_norm <- inputMatrix %>% mutate(chromArm=scCNVCaller$cytoband_data$V4,cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(chrom=str_c(chrom,chromArm)) %>% filter(cpg!=0)
  if (binExpand>1)
  {
    #collapse bins
    scData_k_norm <- scaleMatrixBins(scData_k_norm,binSizeRatio = binExpand,"pos")
    #rename column
    scData_k_norm <- scData_k_norm %>% rename(pos = pos_b)
  }
  # %>% filter(cpg!=0)
  #gender determination
  sckn<-data.table(scData_k_norm)
  total_cpg<-sckn[,summaryFunction(cpg),by="chrom"]
  sckn[,c("blacklist","raw_medians","chromArm","cpg","pos"):=NULL]
  tail(colnames(sckn))
  setkey(sckn,chrom)
  sckn<-sckn[c("chrXq","chrYq"),lapply(.SD,quantile,probs=0.8),by="chrom"]
  nrow(sckn)
  total_cpg[chrom %in% c("chrXq","chrYq")]
  sckn[,cpg:=total_cpg[chrom %in% c("chrXq","chrYq")]$V1]
  #sckn[,print(.SD)]
  xy_signal<-transpose(sckn[,lapply(.SD,"/",1+cpg),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
  
  #t(sckn)[2:ncol(sckn),]
  #xy_signal <- as_tibble(apply(t(sckn)[2:ncol(sckn),],2,as.numeric))
  #xy_signal
  #xy_signal <-scData_k_norm %>% filter(chrom %in% c("chrXq","chrYq")) %>% group_by(chrom) %>% summarise_if(is.numeric,list(quantile),probs=0.8) 
  sexCutoff=5e-7
  logTrans=FALSE
  if (logTrans)
  {
    xy_signal<-transpose(sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
    
    #xy_signal <- xy_signal %>% rowwise() %>% mutate(totalcpg=summaryFunction(cpg)) %>% select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(totalcpg,base=2)))) 
  } else
  {
    xy_signal<-transpose(sckn[,lapply(.SD,"/",1+cpg),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
    
    #xy_signal <- xy_signal %>% rowwise() %>% mutate(totalcpg=summaryFunction(cpg)) %>% select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+totalcpg))) 
    sexCutoff=0.25
  }
  print(xy_signal)
  diffXY<-xy_signal$chrXq-xy_signal$chrYq
  print(abs(diffXY))
  
  if (abs(diffXY)<sexCutoff)
  {
    isMale=TRUE
  }
  #top part works awesome
  
  #now to convert rest
  sckn<-data.table(scData_k_norm)
  print(nrow(sckn))
  sckn <- sckn[chromArm!="cen"][blacklist==FALSE]
  print(nrow(sckn))
  #total_cpg <- sckn[,summaryFunction(cpg),by="chrom"]
  #total_cpg
  #setnames(total_cpg,"V1","totalcpg")
  #remove raw_medians,blacklist,chromArm,cpg
  sckn<-sckn[,c("blacklist","raw_medians","chromArm","pos"):=NULL]
  median_chromosome_density<-sckn[1:(nrow(sckn)),lapply(.SD,summaryFunction),by="chrom"]
  med_temp<-transpose(median_chromosome_density[,c("cpg"):=NULL],make.names="chrom")
  #head(med_temp)
  transpose(med_temp[,lapply(.SD,summaryFunction)],keep.names="chrom")
  median_chrom_signal<-transpose(med_temp[,lapply(.SD,summaryFunction)],keep.names="chrom")
  print(median_chrom_signal)
  
  
  #median_chromosome_density[,1:5]
  #]
  #median_chromosome_density<-scData_k_norm %>%
  #  group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) 
  message(nrow(scData_k_norm))
  #now append 
  #summaryFunction
  
  
  total_cpg<-sckn[,summaryFunction(cpg),by="chrom"]
  
  #TODO: cpg of 200 and medianNorm cutoffs only work for 1e6 bins; need to change for others
  if (logTrans==TRUE)
  {
    #5.1e-6 old 
    #5.1e-06, totalcpg 130 for 200k
    #THIS ONE
    #transpose(sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL]
    sckn<-sckn[,lapply(.SD,sum),by="chrom"][,cpg:=total_cpg$V1][,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"]
    #sckn<-sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL][,lapply(.SD,sum),by="chrom"]
    #scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(cpg,base=2)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    #  select(-cpg) %>% filter(chrom %in% median_chrom_signal$chrom) %>% 
    #  mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>minimumChromValue) 
    #+log(totalcpg,base=2)
    #was 200
  }
  else
  {
    #THIS ONE
    #[,lapply(.SD,"/",1+cpg),by="chrom"][,cpg:=NULL]
  
    sckn<-sckn[,lapply(.SD,sum),by="chrom"][,cpg:=total_cpg$V1][,lapply(.SD,"/",1+cpg),by="chrom"]
    #scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+sqrt(cpg)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    #  select(-cpg) %>% filter(chrom %in% median_chrom_signal$chrom) %>%
    #  mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>minimumChromValue) 
  }
  
  sckn<-sckn[,medianNorm:=median_chrom_signal$V1]
  sckn<-sckn[which(sckn[,medianNorm>minimumChromValue])]
  #return the appropriate matrices
  return(as.data.frame(sckn))
}
#' collapseChrom3
#'
#' Deprecated
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @examples
#' cat_function()
collapseChrom3<-function(inputMatrix,minimumSegments=40,summaryFunction=cutAverage,logTrans=TRUE,binExpand=1,minimumChromValue=2)
{
  #remove centromeric bins
  #these are ordered in the same fashion so should be accurate
  #divide signal by cpg density + 1 (To avoid dividing by zero)
  #remove blacklist regions
  # %>% filter(cpg!=0)
  scData_k_norm <- inputMatrix %>% mutate(chromArm=scCNVCaller$cytoband_data$V4,cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(chrom=str_c(chrom,chromArm)) %>% filter(cpg!=0)
  if (binExpand>1)
  {
    #collapse bins
    scData_k_norm <- scaleMatrixBins(scData_k_norm,binSizeRatio = binExpand,"pos")
    #rename column
    scData_k_norm <- scData_k_norm %>% rename(pos = pos_b)
  }
  # %>% filter(cpg!=0)
  #gender determination
  xy_signal <-scData_k_norm %>% filter(chrom %in% c("chrXq","chrYq")) %>% group_by(chrom) %>% summarise_if(is.numeric,list(quantile),probs=0.8) 
  sexCutoff=5e-7
  
  if (logTrans)
  {
      xy_signal <- xy_signal %>% rowwise() %>% mutate(totalcpg=summaryFunction(cpg)) %>% select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(totalcpg,base=2)))) 
  }
  else
  {
    xy_signal <- xy_signal %>% rowwise() %>% mutate(totalcpg=summaryFunction(cpg)) %>% select(-cpg) %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+totalcpg))) 
    sexCutoff=0.25
  }
  print(xy_signal)
  XYmedians <- xy_signal %>% gather(Cell,Density,ends_with(scCNVCaller$cellSuffix)) %>% group_by(chrom) %>% summarise(med=quantile(Density,probs=0.8))
  diffXY<-XYmedians$med[1]-XYmedians$med[2]
  print(abs(diffXY))
  print(XYmedians$med)
  if (abs(diffXY)<sexCutoff)
  {
    isMale=TRUE
  }
  scData_k_norm <- scData_k_norm %>% mutate(isCentromere=str_detect(chrom,"cen")) %>% filter(isCentromere==FALSE,blacklist==FALSE) %>% select(-blacklist)
  
  message(nrow(scData_k_norm))
  #now append 
  #summaryFunction
  total_cpg<-scData_k_norm %>% select(chrom,cpg) %>%
    group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction))
  #print(tail(total_cpg))
  
  #recompute medians after normalization
  median_chromosome_density<-scData_k_norm %>%
    group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) 
  #can we normalize signal to X chromosome signal?
  
  median_chrom_signal<-median_chromosome_density %>% gather(Cell,Value,2:ncol(median_chromosome_density)) %>% select(-Cell) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction))
  print(median_chrom_signal)
  median_chrom_signal<-median_chrom_signal %>% filter(!is.na(Value))
  print(median_chrom_signal)
  
  #TODO: cpg of 200 and medianNorm cutoffs only work for 1e6 bins; need to change for others
  if (logTrans==TRUE)
  {
    #5.1e-6 old 
    #5.1e-06, totalcpg 130 for 200k
    scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(cpg,base=2)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
      select(-cpg) %>% filter(chrom %in% median_chrom_signal$chrom) %>% 
      mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>minimumChromValue) 
    #+log(totalcpg,base=2)
    #was 200
  }
  else
  {
    scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+sqrt(cpg)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
      select(-cpg) %>% filter(chrom %in% median_chrom_signal$chrom) %>%
      mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>minimumChromValue) 
  }
  
  #return the appropriate matrices
  return(scData_prechrom)
}
#STEP 1:  bin signal per chromosome
#' filterCells
#'
#' This function filters the collapsed chromosome signal density matrix based on predefined parameters (minimum non-zero chromosome segments)
#' @param inputMatrix - collapsed, normalised matrix
#' @param minimumSegments  the minimum segments required to keep a cell
#' @param minDensity  minimum segment density to consider useful (defaults to zero)
#' @keywords filter
#' @keywords CNV
#' @export
#' @examples
#' filterCells(inputMatrix,minimumSegments=41)
filterCells <- function(inputMatrix,minimumSegments=40,minDensity=0)
{
  tmp<-inputMatrix
  cellQuality<-tmp %>% 
    gather(Cell,Density,2:ncol(tmp)) %>% 
    spread(chrom,Density) 
  print(head(cellQuality))
  cellQuality <- cellQuality %>% mutate(count=rowSums(.[2:ncol(cellQuality)]>minDensity))
  cellQuality$count
  barcodesPassing <- (cellQuality %>% filter(count>=minimumSegments))$Cell
  print(str_c(length(barcodesPassing)," barcodes passed quality check"))
  
  #filter scData_prechrom for quality checks
  tmp <- tmp %>% select(c("chrom",barcodesPassing))
  
  #average total signal per cell
  aveSignal <- tmp %>% summarize_at(vars(ends_with(scCNVCaller$cellSuffix)),list(sum)) %>% gather(Cell,Average)
  se_aveSignal <- sd(aveSignal$Average)/sqrt(length(aveSignal$Average))
  sd_aveSignal <- sd(aveSignal$Average)
  mean_aveSignal <- mean(aveSignal$Average)
  sd_aveSignal
  mean_aveSignal - 2*sd_aveSignal
  #select cells
  cellsPassingB <- aveSignal %>% filter(between(Average,mean_aveSignal - 2*sd_aveSignal,mean_aveSignal + 2*sd_aveSignal)) %>% select(Cell)
  print(str_c(length(barcodesPassing)," barcodes passed quality check 2"))
  scData_chrom_filtered <- tmp %>% select(c("chrom",cellsPassingB$Cell))
  scCNVCaller$finalFilterCells<-length(barcodesPassing)
  return(scData_chrom_filtered)
}

#recompute range values after normalization and filtering
#' computeCenters
#'
#' Recompute range values after normalization and filtering. Returns a list containing values with summary values and IQRs
#' @param inputMatrix - normalised, filtered matrix
#' @param summaryFunction - summary function to use (defaults to cutAverage)
#' @keywords median
#' @keywords CNV
#' @export
#' @examples
#' listCenters <- computeCenters(inputMatrix,cutAverage)
computeCenters <- function(inputMatrix,summaryFunction=cutAverage)
{
  DT1<-data.table(Chrom=inputMatrix$chrom,Cells=inputMatrix %>% select(ends_with(scCNVCaller$cellSuffix)))
  DT1b<-transpose(DT1,make.names = "Chrom")
  #DT1b<-DT1b[1:(nrow(DT1b)),lapply(.SD,summaryFunction)]
  #DT1b2<-transpose(DT1b,keep.names ="Chrom")
  #median_chrom_signal_filter<-DT1b2[,.(Median=median(V1))]
  #IQR_chrom_signal_filter<-DT1b2[,.(IQR=IQR(V1))]
  DT1b2<-DT1b[1:(nrow(DT1b)),lapply(.SD,summaryFunction)]
  DT1b_IQR<-DT1b[1:(nrow(DT1b)),lapply(.SD,IQR)]
  #rownames_to_column(data.frame(Value=t(DT1b_IQR)),"chrom")

  return(list(
    rownames_to_column(data.frame(Value=t(DT1b2)),"chrom"),
    rownames_to_column(data.frame(Value=t(DT1b_IQR)),"chrom"))
    )
}
#' scaleMatrix
#'
#' Scales the normalised, filtered input matrix based on provided median and IQR data
#' @param inputMatrix normalised, filtered matrix
#' @param median_iqr_list  list containing summary lists with median/trimmed mean and IQR (computed by computeCenters)
#' @param spread   return the list in Tidyverse "spread" format - for plotting
#' @keywords CNV
#' @keywords scale
#' @export
#' @examples
#' scaleMatrix(inputMatrix,medianIQR)
scaleMatrix <- function(inputMatrix,median_iqr_list,spread=FALSE)
{
  #FINAL FILTER: remove Y chromosome
  tmp<-inputMatrix
  median_list<-median_iqr_list[[1]]
  iqr_list<-median_iqr_list[[2]]
  print(median_list)
  print(iqr_list)
  #rep(median(median_list$Value),times=length(median_list$Value))
  tmp[2:ncol(tmp)]=t(scale(t(tmp[2:ncol(tmp)]),center=median_list$Value,scale=iqr_list$Value))
  tmp <- tmp %>% filter(chrom!="chrYp") 
  if (spread==TRUE)
  {
    scData_chrom_spread <- tmp %>% gather(Cell,Density,2:ncol(tmp))
    return(scData_chrom_spread)
  }
  else
  {
    return(tmp)
  }
}


#now release results

#list of chromosomes and cluster means

#make fake cells?
#' makeFakeCells
#'
#' simulate a pseudodiploid control
#' @param x (mean of values)
#' @param num (number of cells)
#' @param sd_fact (scaling factor for cell SD)
#' @export
makeFakeCells<-function(x,num=300,sd_fact=0.1){
  a<-rnorm(n=num,mean=x,sd=sd_fact*x)
  return(a)
}
#deltaMean used to collapse clusters initially
#deltaBICs 
#do clustering
#uncertainty cutoff of less than 1 cuts off more and more cells
#' identifyCNVClusters
#'
#' How the sausage is made - uses Gaussian decomposition to analyse a filtered, normalised input matrix and call putative CNVs
#' @param x (mean of values)
#' @param num (number of cells)
#' @param sd_fact (scaling factor for cell SD)
#' @export
identifyCNVClusters <- function(inputMatrix, median_iqr, useDummyCells = FALSE,propDummy=0.1, maxClust=4,deltaMean=0.10, minDiff=0.25,deltaBIC2=50, bicMinimum=5,minMix=0.3,subsetSize=500,fakeCellSD=0.2,uncertaintyCutoff=1, summaryFunction=cutAverage)
{
  scData_chrom<-inputMatrix
  median_chrom_signal<-median_iqr[[1]]
  median_chrom_signal_filter <- median_chrom_signal %>% filter(chrom %in% scData_chrom$chrom) %>% filter(chrom!="chrX",chrom!="chrY")
  IQR_chrom_signal<-median_iqr[[2]]
  IQR_chrom_signal_filter <- IQR_chrom_signal %>% filter(chrom %in% median_chrom_signal_filter$chrom)
  #
  median_chrom_signal_filter$Value<-median(median_chrom_signal_filter$Value)
  #IQR_chrom_signal_filter$Value<-quantile(IQR_chrom_signal_filter$Value,0.25)
  #IQR_chrom_signal_filter$Value<-1/binSize
  #generate simulated cells
  fakeCellNum=floor(ncol(scData_chrom)*propDummy)
  fakeCells<-sapply(median_chrom_signal_filter$Value,makeFakeCells,fakeCellNum,fakeCellSD)
  fakeCellsSorted<-cbind.data.frame(chrom=median_chrom_signal_filter$chrom,t(fakeCells),stringsAsFactors=FALSE)
  colnames(fakeCellsSorted)[2:ncol(fakeCellsSorted)]=str_c("X",colnames(fakeCellsSorted)[2:ncol(fakeCellsSorted)])
  
  #generate scaffold with simulated cells
  if (useDummyCells)
  {
    scData_chrom2 <- left_join(inputMatrix,fakeCellsSorted,by="chrom")
    scData_chrom2 <- scData_chrom2 %>% filter(chrom!="chrYp")
  }
  else
  {
    scData_chrom2 <- scData_chrom
  }
  #summaryFunction
  
  med_iqr2<-computeCenters(inputMatrix = scData_chrom2,summaryFun = summaryFunction)
  

  median_chrom_signal_filter <- med_iqr2[[1]] %>% filter(Value>0)
  IQR_chrom_signal_filter <- med_iqr2[[2]] %>% filter(chrom %in% median_chrom_signal_filter$chrom)

  #clean up again
  median_chrom_signal_filter$Value<-median(median_chrom_signal_filter$Value)
  IQR_chrom_signal_filter$Value<-quantile(IQR_chrom_signal_filter$Value,0.25)
  print(median_chrom_signal_filter)
  print(IQR_chrom_signal_filter)
  #add simulated cells if needed
  
  if (useDummyCells) {
    # scData_chrom2[2:ncol(scData_chrom2)]=t(scale(t(scData_chrom2[2:ncol(scData_chrom2)]),center=rep(median(median_chrom_signal_filter$Value),times=length(median_chrom_signal_filter$Value)),scale=rep(median(IQR_chrom_signal_filter$Value),times=length(median_chrom_signal_filter$Value))))
    #  scData_chrom2 <- scData_chrom2 %>% filter(chrom!="chrYp") 
    #  scData_chrom_spread <- scData_chrom2 %>% gather(Cell,Density,2:ncol(scData_chrom2))
    scData_chrom_spread <- scaleMatrix(scData_chrom2,median_iqr_list = list(median_chrom_signal_filter,IQR_chrom_signal_filter),spread = TRUE)
  }  else {
    # scData_chrom <- scData_chrom %>% filter(chrom!="chrYp") 
    # scData_chrom_spread <- scData_chrom %>% gather(Cell,Density,2:ncol(scData_chrom))
    scData_chrom_spread <- scaleMatrix(scData_chrom2,median_iqr_list = list(median_chrom_signal_filter,IQR_chrom_signal_filter),spread = TRUE)
  }
  #NOW DO GAUSSIAN CLUSTERING
  print(ggplot(scData_chrom_spread,aes(chrom,Density))+geom_violin(scale="width") + theme(axis.text.x=element_text(angle=-90,hjust = 0,vjust=0.5))) 
  
  #initialize variables for clustering
  cell_ids=unique(scData_chrom_spread$Cell)
  chrom_clusters<-data.frame(Chrom=median_chrom_signal_filter$chrom,V1=0,V2=0,V3=0,V4=0,V5=0,V6=0,stringsAsFactors=FALSE)
  chrom_clusters
  cell_assignments<-data.frame(Cells=cell_ids,stringsAsFactors = FALSE)
  
  #initialization
  set.seed(0)
  defSubset=sample(1:length(unique(scData_chrom_spread$Cell)),subsetSize)
  defSubset
  #can cut down if the mixing proportion minimum is less than X% (e.g. 20)
  for (m5 in 1:(length(median_chrom_signal_filter$chrom)))
  {
    set.seed(0)
    
    #hc doesn't wo
    dens1<-scData_chrom_spread %>% filter(chrom==median_chrom_signal_filter$chrom[m5]) %>% select(Density)
    initParams=list(noise=TRUE,subset=defSubset) #hcPairs=randomPairs(dens1$Density,modelName = "V")) #,noise=sample(1:length(dens1$Density),10)) 
    message(median_chrom_signal_filter$chrom[m5])
    fit = Mclust(dens1$Density,G=1:(maxClust-1),model=c("V"), prior = priorControl(), initialization = initParams) 
    #retype if greater than two components
    if (fit$G>2)
    {
      if ((fit$BIC[fit$G]-fit$BIC[fit$G-1])<deltaBIC2)
      {
        fit = Mclust(dens1$Density,G=1:2,prior = priorControl(),model="V",initialization = initParams)
      }
    }
    if (fit$G==2)
    {
      #if delta BIC less than X do it
      deltaBIC = fit$BIC[2]-fit$BIC[1]
      #10-15 BIC minimum
      #or smaller pop less than predefined percentage
      arrangedMix <- fit$parameters$pro[order(fit$parameters$pro)]
      propCutoff=minMix
      minProp = 0
      if (0 %in% fit$classification)
      {
        minProp = arrangedMix[2]
      }
      else
      {
        minProp = arrangedMix[3]
      }
      if ((deltaBIC < bicMinimum) | (minProp < propCutoff))
      {
        message(median_chrom_signal_filter$chrom[m5])
        
        #redo as single gaussian
        fit = Mclust(dens1$Density,G=1,prior = priorControl(),model="V",initialization = initParams)
      }
      else
      {
        #retype if difference between means is less than cutoff
        diffMean <- (fit$parameters$mean[2]-fit$parameters$mean[1])
        
        varChange <- 0
        if (fit$modelName=="V")
        {
          arrangeVar <- fit$parameters$variance$scale[order(fit$parameters$variance$scale)]
          varChange <- (arrangeVar[2]/arrangeVar[1])
          
        }
        #if the larger slops greater than smaller then we need to fix it
        if ((diffMean < deltaMean)) # | ((varChange > 4) & (diffMean < deltaMean * 4)))
        {
          message(median_chrom_signal_filter$chrom[m5])
          
          fit = Mclust(dens1$Density,G=1,model="V",prior = priorControl(),initialization = initParams)
        }
        else
        {
          if ((varChange > 3) & (diffMean<0.5))
          {
            message(str_c("cutting dist", median_chrom_signal_filter$chrom[m5]))
            print(fit$parameters$mean)
            #cut the distribution
            largerSpread <- order(fit$parameters$variance$scale)[2]
            if (largerSpread==2)
            {
              #transfer small 2 to cluster 1
              fit$classification[fit$classification==2 & fit$data<fit$parameters$mean[1]]<-1
            }
            else
            {
              #transfer larger 1 to cluster 2
              fit$classification[fit$classification==1 & fit$data>fit$parameters$mean[2]]<-2
            }
          }
        }
      }
    }
    #test uncertainty
    #print(fit$uncertainty)
    fit$classification[fit$uncertainty>uncertaintyCutoff]<-0
    print(summary(fit))
    #can do adaptations here to check for overfitting
    chrom_clusters[chrom_clusters$Chrom==median_chrom_signal_filter$chrom[m5],2:(2+length(fit$parameters$mean)-1)]<-fit$parameters$mean
    cell_assignments <- cbind.data.frame(cell_assignments,fit$classification,stringsAsFactors=FALSE)
    
    names(cell_assignments)[m5+1]=median_chrom_signal_filter$chrom[m5]
    #TODO: compute LRT statistic
  }
  names(cell_assignments)[2:(length(median_chrom_signal_filter$chrom)+1)]<-median_chrom_signal_filter$chrom
  #try to assign groups based on imputed clusters
  chrom_clusters
  return(list(cell_assignments,chrom_clusters,med_iqr2))
}

#OPTION: make fake cells to use as a normalizing factor (if low non-tumour content)
# #PARAMETERS for clustering
# useDummyCells=FALSE
# #maximum G value
# maxClust=4
# 
# #minimum difference between cluster means in fit
# deltaMean=0.10
# #as above but for merging clusters
# minDiff=0.25
# 
# #minimums for BIC of G>2 and G=2 clusters
# deltaBIC2 = 50
# bicMinimum = 5


#cluster the CNVs
#' clusterCNVs
#'
#' Cleans the putative CNV calls
#' @keywords CNV
#' @export
clusterCNV<-function(initialResultList,medianIQR,maxClust=4,minDiff=0.25){
  #return a list of CNV
  #MERGE AND COLLAPSE CLUSTERS
  chrom_clusters<-initialResultList[[2]]
  cell_assignments<-initialResultList[[1]]
  distinctClusts=1
  chrom_clusters_final<-chrom_clusters
  
  for (m5 in 1:length(medianIQR[[1]]$chrom)){
    clust_list=chrom_clusters %>% filter(Chrom==medianIQR[[1]]$chrom[m5])

    zeroClusters<-length(which(clust_list[2:maxClust]==0))
    numClusts<-maxClust-zeroClusters-1
    #now order vector
    activeClusters<-clust_list[2:(2+numClusts-1)]
    #reorder the keys and values
    #reorder keys
    print(clust_list[2:(2+numClusts-1)])
    clust_list[2:(2+numClusts-1)]<-activeClusters[order(activeClusters)]
    print(clust_list[2:(2+numClusts-1)])
    #copy to final list
    chrom_clusters_final[m5,2:(2+numClusts-1)]<-activeClusters[order(activeClusters)]
    #
    #now reorder the values by the vector using the factor levels trick
    cell_assignments[,m5+1]<-factor(cell_assignments[,m5+1])
    #message(length(levels(cell_assignments[,m5+1])))
    newOrder<-order(activeClusters)
    hasZero<-length(which(cell_assignments[,m5+1]==0))>0
    if (hasZero) {
      newOrder<-c(0,newOrder)
    }
    newOrder
    #TODO: what about zeros
    levels(cell_assignments[,m5+1])<-newOrder
    #tmp_fact<-factor(dm_list_test2$newIndex)
    #
    cell_assignments[,m5+1]<-as.numeric(levels(cell_assignments[,m5+1])[cell_assignments[,m5+1]])
    lastValue=clust_list[2]
    clust_first=1
    distinctClusts=1
    clustLength=1
    ###
    for (i1 in 3:(1+numClusts)){
  #    message(str_c(clust_list[i1],i1,sep=","))
      if (clust_list[i1]!=0)
      {
        #assigned cluster
        clustDiff=abs(clust_list[i1]-lastValue)
   #     message(str_c(clustDiff,distinctClusts,m5,"boo",sep=","))
        if (clustDiff<minDiff)  {
          #collapse clusters
          #
          #merge i1 to prevIndex
          #lastValue=(lastValue + clust_list[i1])/2
          clustLength = clustLength +1
          lastValue=(lastValue*(clustLength-1) + clust_list[i1])/clustLength
          #can rechagne thi
          for (cMembers in (i1-clustLength+1):i1)
          {
            clust_list[cMembers]<-lastValue
          }
          chrom_clusters_final[m5,i1]<-chrom_clusters_final[m5,clust_first+1]
   #       message(str_c(chrom_clusters_final[m5,i1],lastValue,i1,sep=","))
          cell_assignments[cell_assignments[,m5+1]==(i1-1),m5+1]=clust_first
          #clustMembers = clustMembers + 1
        }
        else {
          #save the previous cluster
     #     message("new")
          distinctClusts = distinctClusts+1
          clust_first=i1-1
          lastValue=clust_list[clust_first+1]
        }
      }
    }
    #refactor the list - or could use dense_rank - todo later
    minCellAssign=1
    if (hasZero){
      minCellAssign=0
    }
    #message(length(le,vels(cell_assignments[,m5+1])))
    cell_assignments[,m5+1]<-factor(cell_assignments[,m5+1])
   # message(levels(cell_assignments[,m5+1]))
    levels(cell_assignments[,m5+1])<-seq(minCellAssign,distinctClusts)
    #tmp_fact<-factor(dm_list_test2$newIndex)
    
    uniq_clusts<-unique(as.numeric(chrom_clusters_final[m5,2:ncol(chrom_clusters_final)]))
    chrom_clusters_final[m5,2:ncol(chrom_clusters_final)]<-0
     chrom_clusters_final[m5,2:(1+length(uniq_clusts))]<-uniq_clusts
    #
   # message(levels(cell_assignments[,m5+1]))
    cell_assignments[,m5+1]<-as.numeric(levels(cell_assignments[,m5+1])[cell_assignments[,m5+1]])
    #TODO: merge cluster values and medians into a list which we can merge (e.g. chr0_0,0.5))
  }
  chrom_clusters_final
  return(list(cell_assignments,chrom_clusters_final))
}
#'
#' annotateCNV3
#' 
#' Annotates the filtered CNV calls and saves the results, and estimates absolute copy numbers
#' @param cnvResults   output of filtered CNVs (set of lists)
#' @param saveOutput    save output to a file on disk
#' @param maxClust2    maximum # of clusters
#' @param outputSuffix   output suffix to save files
#' @keywords CNV
#' @keywords output
#' @export
#' @examples
#' annotateCNV3(cleanCNV,saveOutput=TRUE,outputSuffix="_nolog")
annotateCNV3 <- function(cnvResults,saveOutput=TRUE,maxClust2=4,outputSuffix="_1")
{
  cell_assignments<-cnvResults[[1]]
  cell_assignments<-column_to_rownames(cell_assignments,var = "Cells")
  chrom_clusters_final<-cnvResults[[2]]
  colnames(chrom_clusters_final)=c("Chrom",seq(1:(maxClust2-2)),"0")
  
  #identify possible CNVs - can use bins to test
  possCNV <- cell_assignments %>% summarise_if(is.numeric,list(max)) %>% gather(Chrom,maxClust) %>% filter(maxClust>1)
  possCNV <- possCNV %>% mutate(Type="Unknown")
  print(cell_assignments %>% summarise_if(is.numeric,list(max)))
  #%>% gather(Chrom,maxClust))
  print(possCNV)
  #output CNV list by cluster
  #cell_assignments %>% filter()
  #WITHOUT BLANKS
  #consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  
  #WITH BLANKS
  consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  
  if(saveOutput==TRUE)
  {
    write.table(x=consensus_CNV_clusters,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  
  colnames(chrom_clusters_final)<-c("Chrom","1","2","0","T1","T2","T3")
  possCNV <- possCNV %>% mutate(normCluster=0)
  XChromCount<-2
  if (isMale)
  {
    XChromCount<-1
  }
  for (suspCNV in possCNV$Chrom)
  {
    #message(possCNV)
    clusterList<-unique(chrom_clusters_final %>% filter(Chrom==suspCNV) %>% select(-Chrom))
    print(clusterList)
    abn<-min(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]]))
    #first item is normal one
    print(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]]))
    print(abn)
    print(which(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]])==abn))
    abn_index<-which(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]])==abn)
    #print(order(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]])))
    #compare to median autosomal value
    #if value is less than autosomal, suspect a loss, otherwise suspect a gain to be abnormal
    possCNV
    print(as.numeric(clusterList[abn_index]))
    possCNV$normCluster[possCNV$Chrom==suspCNV] =abn_index
    print(possCNV)
  }
  #merge cell assignments and annotate
  v1<-rownames_to_column(cell_assignments) %>% filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% filter(Chrom %in% possCNV$Chrom)
  v2<-left_join(v1,possCNV,by="Chrom") %>% mutate(clust_text=if_else(clust>normCluster,3,if_else(clust<normCluster,1,2))) %>% mutate(clust_text=if_else(clust==0,2,clust_text))
  v2 <- v2 %>% select(rowname,Chrom,clust_text) %>% spread(Chrom,clust_text) %>% arrange(rowname)
  v2
  if(saveOutput==TRUE)
  {
    write.table(x=v2,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv_anno_binary.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  colnames(possCNV)[1]<-"chrom"
  
  final_cnv_binarized<-inner_join(consensus_CNV_clusters %>% gather(chrom,clust,2:ncol(consensus_CNV_clusters)),possCNV,by="chrom") %>% mutate(clust=if_else(clust==normCluster,0,1)) %>% select(rowname,chrom,clust) %>% spread(chrom,clust)
  return(list(v2,possCNV,final_cnv_binarized))
}
#' graphCNVDistribution
#'
#' Makes graph of normalised, collapsed fragment data
#' @param inputMatrix (the normalised matrix)
#' @param outputSuffix   suffix to append to pdf output file
#' @keywords output
#' @export
#' @examples
#' graphCNVDistribution(inputMatrix,outputSuffix="unfiltered)
graphCNVDistribution <- function(inputMatrix,outputSuffix="_all")
{
  pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_graph.pdf"),width=8,height=6)
  print(ggplot(inputMatrix %>% gather(Cell,Density,2:(ncol(inputMatrix))),aes(chrom,Density))+geom_violin(scale="width") + theme(axis.text.x=element_text(angle=-90)))
  dev.off()
}
#' graphSmoothCNV
#'
#' @keywords output smoothing
#' @export
graphSmoothCNV <- function(clusterOutputs, originalGraph,perplex=20,numNN=0){
  #TODO: cleanup based on geo cluster assignments

  #try plotting without noise cells
  group_clusts<-clusterOutputs[[1]] %>% filter(str_detect(Cells,"X",negate=TRUE))
  set.seed(12345)
  tsn1<-Rtsne(originalGraph,checkDuplicates=FALSE,perplexity=perplex)
  #TODO: convert from factor
  #will do this later
  #https://jmonlong.github.io/Hippocamplus/2018/02/13/tsne-and-clustering/
  est_k<-sqrt(nrow(clusterOutputs[[1]]))
  if (numNN!=0)
  {
    est_k<-numNN
  }
  k = round(est_k)
  if (k %% 2 ==0) {
    k = k + 1
  }
  #in this case Euclidean is appropriate as we are using the tSNE results
  knn.norm = get.knn(as.matrix(tsn1$Y), k = k)
  knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index), 
                                   k), to = as.vector(knn.norm$nn.index), weight = 1/(1 + 
                                                                                        as.vector(knn.norm$nn.dist)))
  nw.norm = graph_from_data_frame(knn.norm, directed = FALSE)
  nw.norm = simplify(nw.norm)
  lc.norm = cluster_louvain(nw.norm)
  
  #candidate_cnvs_clean[[1]]$chr1q
  print(plot(tsn1$Y, col=as.factor(membership(lc.norm)),asp=1))
  new_table<-group_clusts %>% mutate(louvain=as.factor(membership(lc.norm)))
  
  louvain_means<-new_table %>% group_by(louvain) %>% summarise_at(vars(starts_with("chr")),mean)

  louvain_means[2:ncol(louvain_means)]<-louvain_means[2:ncol(louvain_means)]
  louvain_means[2:ncol(louvain_means)]<-apply(louvain_means[2:ncol(louvain_means)],2,round)
  louvain_means

  discretized_cells<-inner_join(new_table %>% select(Cells,louvain),louvain_means)
  #write only relevant output
 # key_chroms<-discretized_cells %>% summarise_if(is.numeric,IQR) %>% gather(chrom,maxval) %>% filter(IQR>1) %>% select(chrom)
 # print(key_chroms)
  #discretized_cells<-discretized_cells %>% select(Cells,louvain,key_chroms$chrom)
  print(discretized_cells)
  
  #NOICE
  write.table(x=discretized_cells  %>% filter(str_detect(Cells,"X",negate=TRUE)),file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_cnv_louvain_nox_50knn.csv"),quote=FALSE,row.names = FALSE,sep=",")
  #save the discretized CNV file
  return(list(louvain_means,discretized_cells,tsn1))
}
#' graphSmoothDM
#'
#' Louvain cluster smooth the DMs
#' @keywords smoothing
#' @export
graphSmoothDM <- function(clusterOutputs,origData,dm_outputs,numNN=0,perplex=20){
  #TODO: cleanup based on geo cluster assignments
  #try plotting without noise cells
  dm_outputs_t<-dm_outputs
  if (!is.numeric(dm_outputs_t[2]))
  {
     dm_outputs_t[2:ncol(dm_outputs)]<-apply(dm_outputs_t[2:ncol(dm_outputs)],2,as.numeric)
  }
  print(dm_outputs_t)
  group_clusts<-clusterOutputs[[1]] %>% filter(str_detect(Cells,"X",negate=TRUE))
  #group_clusts<-inner_join(group_clusts,dm_outputs_t,by="Cells")
 # group_clusts[is.na(group_clusts)]<-0
  print(group_clusts)
  tsn1<-Rtsne(origData,check_duplicates=FALSE,perplexity=perplex)
  #TODO: convert from factor
  #will do this later
  #https://jmonlong.github.io/Hippocamplus/2018/02/13/tsne-and-clustering/
  
  est_k<-sqrt(nrow(clusterOutputs[[1]]))
  if (numNN!=0)
  {
    est_k<-numNN
  }
  print(est_k)
  k = round(est_k)
  if (k %% 2 == 0) {
    k = k + 1
  }
  #in this case Euclidean is appropriate as we are using the tSNE results
  knn.norm = get.knn(as.matrix(tsn1$Y), k = k)
  knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index), 
                                   k), to = as.vector(knn.norm$nn.index), weight = 1/(1 + 
                                                                                        as.vector(knn.norm$nn.dist)))
  nw.norm = graph_from_data_frame(knn.norm, directed = FALSE)
  nw.norm = simplify(nw.norm)
  lc.norm = cluster_louvain(nw.norm)
  
  #candidate_cnvs_clean[[1]]$chr1q
  print(plot(tsn1$Y, col=as.factor(membership(lc.norm)),asp=1))
  new_table<-group_clusts %>% mutate(louvain=as.factor(membership(lc.norm)))
  #print(nrow(dm_outputs))
  new_table<-inner_join(new_table,dm_outputs,by="Cells")
  print(nrow(new_table))
  louvain_means<-new_table %>% group_by(louvain) %>% summarise_at(vars(starts_with("chr")),mean)
  wiggleFactor<-0
  louvain_means[2:ncol(louvain_means)]<-louvain_means[2:ncol(louvain_means)]+wiggleFactor
  print(louvain_means)
  louvain_means[2:ncol(louvain_means)]<-apply(louvain_means[2:ncol(louvain_means)],2,round)
  print(louvain_means)
  
  discretized_cells<-inner_join(new_table %>% select(Cells,louvain),louvain_means)
  #write only relevant output
  key_chroms<-discretized_cells %>% summarise_if(is.numeric,max) %>% gather(chrom,maxval) %>% select(chrom)
  #print(key_chroms)
  discretized_cells<-discretized_cells %>% select(Cells,louvain,key_chroms$chrom)
  print(discretized_cells)
  #NOICE
  write.table(x=discretized_cells,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_dm_louvain.csv"),quote=FALSE,row.names = FALSE,sep=",")
  #save the discretized CNV file
  return(list(louvain_means,discretized_cells))
}
#LOSS CNV
#' getLOHRegions
#'
#' Calls putative regions of LOH
#' @param inputMatrix - normalised input fragment count matrix
#' @param lossCutoff - minimum Z score for losses (default -0.25)
#' @param uncertaintyCutLoss - maximum uncertainty to accept a cell assignment
#' @param diffThreshold - difference threshold between bins
#' @param minLength - minimum length of a LOH region
#' @param minSeg - minimum number of bins for changepoint analysis
#' @param lossCutoffCells - minimum number of cells to accept a LOH (default: 100)
#' @param targetFun - target function to identify LOH regions
#' @param signalBoost - value to add prior to logtransforming signal
#' @param lossCutoffReads
#' @keywords LOH
#' @keywords CNV
#' @export
getLOHRegions <- function(inputMatrixIn,lossCutoff=(-0.25), uncertaintyCutLoss=0.5, diffThreshold=0.7, minLength=3e6, minSeg=3, lossCutoffCells=100,targetFun=IQR,signalBoost=1e-6,lossCutoffReads=100)
{
  #c("E","V")
  inputMatrix<-inputMatrixIn %>% mutate(cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(arm=scCNVCaller$cytoband_data$V4) %>% filter(arm!="cen", blacklist==0, cpg>0) %>% select(-arm,-blacklist) #%>%  #cpg+
  inputMatrix<-inputMatrix %>% filter(chrom!="chrY") 
  #check if raw_medians
  if ("raw_medians" %in% colnames(inputMatrix))
  {
    inputMatrix <- inputMatrix %>% select(-raw_medians)
  }
  #mutate
  #view(inputMatrix$cpg)
  #inputMatrix<-inputMatrix %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(log(cpg,base=2))))
  #print(nrow(inputMatrix))
  chromList<-unique(inputMatrix$chrom)
  #print(chromList)
  dm_per_cell_vals<-data.frame(cellName=colnames(inputMatrix  %>% select(-chrom,-pos,-cpg)),stringsAsFactors=FALSE)
  alteration_list=c()
  last_coords=c(0,0)
  alteration_delta=c()
  #-0.25
  #0.00007
  
  #1e-4,5e-5
  #diffThreshold<-0
  #chromList=c("chr9")
  #chromList=c("chr5")
  #may not need IQR, may work with mean/median instead
  for (targetChrom in chromList)
  {
    #message(targetChrom)
    IQRv<-apply(inputMatrix  %>% filter(chrom==targetChrom) %>% select(ends_with(scCNVCaller$cellSuffix)),1,quantile,0.3,na.rm=TRUE)
    IQRs<-scale(IQRv,center=TRUE,scale=TRUE)
    expectedSignal<-inputMatrix %>% filter(chrom==targetChrom) %>% select(-cpg,-pos,-chrom) %>% summarise_if(is.numeric,median) %>% gather(barcode,value) %>% summarise_if(is.numeric,mean)
    
    cm<-cpt.meanvar(data=as.vector(IQRs),test.stat="Normal", penalty="AIC",method = "PELT",minseglen = minSeg)
    print(plot(cm,xlab="Chromosome bin",ylim=c(-5,5),ylab="Z-score"))
    cptlist<-t(rbind(cm@param.est$mean,cm@cpts))
    colnames(cptlist)<-c("Mean","Point")
    cptlist<-as_tibble(cptlist) %>% mutate(Diff = Point - lag(Point))
    print(cptlist)
    if (is.na(cptlist$Diff[1]))
    {
    #  message(cptlist$Diff[1])
      cptlist$Diff[1]=cptlist$Point[1]-1
    #  message(cptlist$Diff[1])
    }
   # print(cptlist)
    #set up list of double minutes; cutoff 0.002
    #1.05e-4
    coord_list<-vector(mode = "list")
    d_loss<-cptlist %>% filter(Mean<(lossCutoff))
    d_lossb<-d_loss
   # print(d_loss)
    if (nrow(d_loss)>0)
    {
    tmp_coords<-inputMatrix  %>% filter(chrom==targetChrom) %>% select(chrom,pos)
    d_lossb<-d_loss %>% mutate(startCoord=as.numeric(tmp_coords$pos[Point-Diff]),endCoord=as.numeric(tmp_coords$pos[Point]),startChrom=tmp_coords$chrom[Point])
   # print(d_lossb)
    #d_lossb<-d_lossb %>% filter((endCoord-startCoord)>=minLength) 
    tmpb_merged<-d_lossb %>% mutate(touching=(startCoord==lag(endCoord)))
    tmpb_merged$touching[1]<-FALSE #otherwise NA
    tmpb_merged <- tmpb_merged %>% mutate(startCoord=if_else(touching==TRUE,lag(startCoord),startCoord),Diff=if_else(touching==TRUE,Diff+lag(Diff),Diff))
    #print(length(which(tmpb_merged$touching==TRUE)))
  #  print(tmpb_merged)
   # print(which(tmpb_merged$touching==TRUE)-1)
    if (length(which(tmpb_merged$touching==TRUE)) > 0)
    {
      tmp_merged_final<-tmpb_merged %>% filter(!(row_number() %in% (which(tmpb_merged$touching==TRUE)-1))) %>% select(-touching)
      d_lossb<-tmp_merged_final
    }
    d_lossb <- d_lossb %>% arrange(startChrom,as.numeric(startCoord))
    #print(d_lossb)
    }
     #
      #coord_list<-append(coord_list,list(c(last_coords[1],last_coords[2])))
      #change this to loop
    if (nrow(d_lossb)>=1)
    {
      for (i in 1:nrow(d_lossb))
      {
        last_coords[1]<-d_lossb$startCoord[i]
        last_coords[2]<-d_lossb$endCoord[i]
        alterationName<-str_c(targetChrom,sprintf(fmt="%d",last_coords[1]),sprintf(fmt="%d",last_coords[2]),sep="_")
        # message(seg_length)
          
        posList=seq(from=last_coords[1],to=last_coords[2],by = 1e6)
        posList<-sprintf(fmt="%d",posList)
        #posList
        #or max
        #print(head(expectedSignal,n=1))
        t2<-inputMatrix %>% filter(chrom==targetChrom,pos %in% posList) %>% select(-cpg,-pos,-chrom) %>% summarise_if(is.numeric,max) # %>% select(-cpg)
        t2d<-t2 %>% gather(Cell,Value)
       # print(t2d)
       # print(inputMatrix %>% filter(chrom==targetChrom,pos %in% posList) %>% select(-cpg))
       # print(targetChrom)
       # print(posList)
        fit1<-Mclust(-1*log2(t2d$Value+signalBoost),G=2:3,modelNames=c("V"),initialization = list(noise=TRUE))
        if (length(fit1)==0)
        {
           message("broken")
           next()
        }
        message(fit1$G)
          #plot(fit1)
          #if (!is.na(fit1))
          #{
          #print(plot(fit1))
          if (fit1$G>2)
          {
            #collapse clusters
           print(fit1$parameters)
            fit1$classification[fit1$classification==2]<-1
            fit1$classification[fit1$classification==3]<-2
            clust1_mean<-mean(fit1$parameters$mean[1:2])
            clust2_mean<-fit1$parameters$mean[3]
            print(clust2_mean)
            signalDiff=2^(-1*clust2_mean)-expectedSignal$value-signalBoost
                 delta_mean = clust2_mean-clust1_mean
                 #print(clust2_mean-clust1_mean)
            #      delta_mean
            #      message(delta_mean)
                 if (abs(delta_mean)>diffThreshold && (signalDiff<(-1*lossCutoffReads)))
                 {
                    alteration_list<-cbind(alteration_list,alterationName)
                    alteration_delta<-cbind(alteration_delta,signalDiff)
                    fit1$classification[fit1$uncertainty>uncertaintyCutLoss]<-0
            #message(fit1$parameters$mean)
                    dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification)
                  }
            #      #message(fit1$parameters)
          }
          if (fit1$G==2)
          {
            #check for and collapse cluster assignments
            #cells in 2 with val < 1 and vice versa
            fit1$classification[fit1$data<fit1$parameters$mean[1] & fit1$classification==2]<-1
            fit1$classification[fit1$data>fit1$parameters$mean[2] & fit1$classification==1]<-2
            delta_mean = fit1$parameters$mean[2]-fit1$parameters$mean[1]
            #message(str_c(delta_mean,"MEAN"))
            print(delta_mean)
            signalDiff=2^(-1*fit1$parameters$mean[2])-expectedSignal$value-signalBoost
            if (abs(delta_mean)>diffThreshold && (signalDiff<(-1*lossCutoffReads)))
            {
              alteration_list<-cbind(alteration_list,alterationName)
              alteration_delta<-cbind(alteration_delta,signalDiff)
              #add delta mean
              #message(fit1$parameters$mean)
              #uncertainty less than 0.2
              fit1$classification[fit1$uncertainty>uncertaintyCutLoss]<-0
              dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification)
            }
         # }
        #  }
        }
      }
    }
    #cleanup
   
  }
  #final cleanup
  if (length(alteration_list)>0)
  {
  colnames(dm_per_cell_vals)[2:ncol(dm_per_cell_vals)]<-alteration_list
  
  #cutoff for minimum # of cells
  #count # of cells in each cluster
  #rowSums(.[2:ncol(cellQuality)]>0)
  loss_cluster_counts<-dm_per_cell_vals %>% gather(Alteration,Clust,2:ncol(dm_per_cell_vals)) %>% group_by(Alteration) %>% count(Alteration,Clust)
  print(loss_cluster_counts)
  min_loss<-loss_cluster_counts %>% spread(Clust,n) %>% select(-'0') %>% mutate(min=min(`1`,`2`))
  print(min_loss)
  min_loss[is.na(min_loss)]<-0
  #min_loss$Alteration[min_loss$min>lossCutoffCells]
  #cut appropriately
  dm_per_cell_vals <- dm_per_cell_vals %>% select(cellName,min_loss$Alteration[min_loss$min>lossCutoffCells])
  return(list(dm_per_cell_vals,alteration_list,alteration_delta))
  }
  return(list())
}
#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
annotateLosses<-function(lossFile)
{
  coords<-t(as.data.frame(sapply(lossFile[[2]],str_split,"[._]"),stringsAsFactors = FALSE))
  coords
  colnames(coords)<-c("chrom","start","end")
  coords_final<-as_tibble(coords) %>% select(chrom,start,end) %>% mutate(genes="")
  #annotate
  ensembl <- useMart("ensembl")
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  filts=c('chromosome_name','start','end')
  for (i4 in 1:nrow(coords_final))
  {
    chrname=str_remove(coords_final$chrom[i4],"chr")
    vals=list(chromosome_name=chrname,start=coords_final$start[i4],end=coords_final$end[i4])
    genes_found<-getBM(attributes=c('hgnc_symbol','transcript_length'), 
                       filters = filts, 
                       values = vals, 
                       mart = ensembl)
    coords_final$genes[i4]<-str_c(unique((genes_found %>% filter(transcript_length>400) %>% arrange(desc(transcript_length)))$hgnc_symbol),collapse=",")
  }
  write.table(x=coords_final,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_loss_annotations.csv"),quote=FALSE,row.names = FALSE,sep=",")
  return(coords_final)
  
}
#' getLOHRegions2
#'
#' LOH analysis - alternate approach
#' @keywords LOH
#' @export
getLOHRegions2 <- function(inputMatrixIn,lossCutoff=(-0.25), uncertaintyCutLoss=0.5, diffThreshold=0.7, minLength=3e6, minSeg=3, lossCutoffCells=100,targetFun=IQR,signalBoost=1e-6)
{
  #c("E","V")
  inputMatrix<-inputMatrixIn %>% mutate(cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(arm=scCNVCaller$cytoband_data$V4) %>% filter(arm!="cen", blacklist==0) %>% select(-arm,-blacklist) #%>%  #cpg+
  inputMatrix<-inputMatrix %>% filter(chrom!="chrY") 
  #check if raw_medians
  if ("raw_medians" %in% colnames(inputMatrix))
  {
    inputMatrix <- inputMatrix %>% select(-raw_medians)
  }
  #mutate
  #view(inputMatrix$cpg)
  #inputMatrix<-inputMatrix %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(log(cpg,base=2))))
  print(nrow(inputMatrix))
  chromList<-unique(inputMatrix$chrom)
  print(chromList)
  dm_per_cell_vals<-data.frame(cellName=colnames(inputMatrix  %>% select(-chrom,-pos,-cpg)),stringsAsFactors=FALSE)
  alteration_list=c()
  last_coords=c(0,0)
  alteration_delta=c()
  #-0.25
  #0.00007
  
  #1e-4,5e-5
  #diffThreshold<-0
  #chromList=c("chr9")
  #chromList=c("chr5")
  #may not need IQR, may work with mean/median instead
  for (targetChrom in chromList)
  {
    #message(targetChrom)
    IQRv<-apply(inputMatrix  %>% filter(chrom==targetChrom) %>% select(ends_with(scCNVCaller$cellSuffix)),1,quantile,0.4,na.rm=TRUE)
    IQRs<-scale(IQRv,center=TRUE,scale=TRUE)
    
    cm<-cpt.meanvar(data=as.vector(IQRs),test.stat="Normal", penalty="AIC",method = "PELT",minseglen = minSeg)
    print(plot(cm))
    cptlist<-t(rbind(cm@param.est$mean,cm@cpts))
    colnames(cptlist)<-c("Mean","Point")
    cptlist<-as_tibble(cptlist) %>% mutate(Diff = Point - lag(Point))
    print(cptlist)
    if (is.na(cptlist$Diff[1]))
    {
      #  message(cptlist$Diff[1])
      cptlist$Diff[1]=cptlist$Point[1]-1
      #  message(cptlist$Diff[1])
    }
    # print(cptlist)
    #set up list of double minutes; cutoff 0.002
    #1.05e-4
    coord_list<-vector(mode = "list")
    d_loss<-cptlist %>% filter(Mean<(lossCutoff))
    d_lossb<-d_loss
    # print(d_loss)
    if (nrow(d_loss)>0)
    {
      tmp_coords<-inputMatrix  %>% filter(chrom==targetChrom) %>% select(chrom,pos)
      d_lossb<-d_loss %>% mutate(startCoord=as.numeric(tmp_coords$pos[Point-Diff]),endCoord=as.numeric(tmp_coords$pos[Point]),startChrom=tmp_coords$chrom[Point])
      # print(d_lossb)
      #d_lossb<-d_lossb %>% filter((endCoord-startCoord)>=minLength) 
      tmpb_merged<-d_lossb %>% mutate(touching=(startCoord==lag(endCoord)))
      tmpb_merged$touching[1]<-FALSE #otherwise NA
      tmpb_merged <- tmpb_merged %>% mutate(startCoord=if_else(touching==TRUE,lag(startCoord),startCoord),Diff=if_else(touching==TRUE,Diff+lag(Diff),Diff))
      print(length(which(tmpb_merged$touching==TRUE)))
      #  print(tmpb_merged)
      # print(which(tmpb_merged$touching==TRUE)-1)
      if (length(which(tmpb_merged$touching==TRUE)) > 0)
      {
        tmp_merged_final<-tmpb_merged %>% filter(!(row_number() %in% (which(tmpb_merged$touching==TRUE)-1))) %>% select(-touching)
        d_lossb<-tmp_merged_final
      }
      d_lossb <- d_lossb %>% arrange(startChrom,as.numeric(startCoord))
      print(d_lossb)
    }
    #
    #coord_list<-append(coord_list,list(c(last_coords[1],last_coords[2])))
    #change this to loop
    if (nrow(d_lossb)>=1)
    {
      for (i in 1:nrow(d_lossb))
      {
        last_coords[1]<-d_lossb$startCoord[i]
        last_coords[2]<-d_lossb$endCoord[i]
        alterationName<-str_c(targetChrom,sprintf(fmt="%d",last_coords[1]),sprintf(fmt="%d",last_coords[2]),sep="_")
        # message(seg_length)
        
        posList=seq(from=last_coords[1],to=last_coords[2],by = 1e6)
        posList<-sprintf(fmt="%d",posList)
        #posList
        #or max
        t2<-inputMatrix %>% filter(chrom==targetChrom,pos %in% posList) %>% select(-cpg,-pos,-chrom) %>% summarise_if(is.numeric,sum) # %>% select(-cpg)
        t2d<-t2 %>% gather(Cell,Value)
        print(t2d)
        # print(inputMatrix %>% filter(chrom==targetChrom,pos %in% posList) %>% select(-cpg))
        # print(targetChrom)
        # print(posList)
        fit1<-Mclust(-1*log2(t2d$Value+signalBoost),G=2:3,modelNames=c("V"),initialization = list(noise=TRUE))
        if (length(fit1)==0)
        {
          message("broken")
          next()
        }
        message(fit1$G)
        #plot(fit1)
        #if (!is.na(fit1))
        #{
        #print(plot(fit1))
        if (fit1$G>2)
        {
          #collapse clusters
          print(fit1$parameters)
          fit1$classification[fit1$classification==2]<-1
          fit1$classification[fit1$classification==3]<-2
          clust1_mean<-mean(fit1$parameters$mean[1:2])
          clust2_mean<-fit1$parameters$mean[3]
          delta_mean = abs(clust2_mean-clust1_mean)
          #      delta_mean
          #      message(delta_mean)
          if (delta_mean>diffThreshold)
          {
            alteration_list<-cbind(alteration_list,alterationName)
            alteration_delta<-cbind(alteration_delta,delta_mean)
            fit1$classification[fit1$uncertainty>uncertaintyCutLoss]<-0
            #message(fit1$parameters$mean)
            dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification)
          }
          #      #message(fit1$parameters)
        }
        if (fit1$G==2)
        {
          #check for and collapse cluster assignments
          #cells in 2 with val < 1 and vice versa
          fit1$classification[fit1$data<fit1$parameters$mean[1] & fit1$classification==2]<-1
          fit1$classification[fit1$data>fit1$parameters$mean[2] & fit1$classification==1]<-2
          delta_mean = abs(fit1$parameters$mean[2]-fit1$parameters$mean[1])
          message(str_c(delta_mean,"MEAN"))
          if (delta_mean>diffThreshold)
          {
            alteration_list<-cbind(alteration_list,alterationName)
            alteration_delta<-cbind(alteration_delta,delta_mean)
            #add delta mean
            #message(fit1$parameters$mean)
            #uncertainty less than 0.2
            fit1$classification[fit1$uncertainty>uncertaintyCutLoss]<-0
            dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification)
          }
          # }
          #  }
        }
      }
    }
    #cleanup
    
  }
  #final cleanup
  if (length(alteration_list)>0)
  {
    colnames(dm_per_cell_vals)[2:ncol(dm_per_cell_vals)]<-alteration_list
    
    #cutoff for minimum # of cells
    #count # of cells in each cluster
    #rowSums(.[2:ncol(cellQuality)]>0)
    loss_cluster_counts<-dm_per_cell_vals %>% gather(Alteration,Clust,2:ncol(dm_per_cell_vals)) %>% group_by(Alteration) %>% count(Alteration,Clust)
    print(loss_cluster_counts)
    min_loss<-loss_cluster_counts %>% spread(Clust,n) %>% select(-'0') %>% mutate(min=min(`1`,`2`))
    print(min_loss)
    min_loss[is.na(min_loss)]<-0
    #min_loss$Alteration[min_loss$min>lossCutoffCells]
    #cut appropriately
    dm_per_cell_vals <- dm_per_cell_vals %>% select(cellName,min_loss$Alteration[min_loss$min>lossCutoffCells])
    return(list(dm_per_cell_vals,alteration_list,alteration_delta))
  }
  return(list())
}

#' readInputTable
#' Read input counts fragment matrix
#' @param inputFile   input count fragment matrix
#' @param sep   separator (defaults to tab)
#' @keywords import
#' @keywords CNV
#' @export
readInputTable<-function(inputFile,sep="\t")
{
  scData<-fread(inputFile,header=TRUE,stringsAsFactors = FALSE,sep="\t",data.table = FALSE)
  #TODO: Add masking
  #scData<-scData %>% mutate(Cell_id = paste(Cell_id,scCNVCaller$cellSuffix,sep="")) 
  #scData$Cell_id
  colnames(scData)[1]<-"Cell_id"
  
  rownames(scData)<-scData$Cell_id
  scData<-scData %>% select(-Cell_id)
  return(scData)
}
