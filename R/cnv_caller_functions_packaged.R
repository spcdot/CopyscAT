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
  
  #initialize chromosome data
  scCNVCaller$chrom_sizes <- read.delim(genomeFile,stringsAsFactors=FALSE,header=FALSE)
  colnames(scCNVCaller$chrom_sizes) <- c("chrom","length")
  
  scCNVCaller$cpg_data<-read.table(cpgFile,stringsAsFactors = FALSE,header=FALSE)
  colnames(scCNVCaller$cpg_data) <- c("chrom","start","end","cpg_density")
  
  #initialize cytoband data
  scCNVCaller$cytoband_data<-read.table(cytobandFile,stringsAsFactors = FALSE,header=FALSE)
  
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
  scCNVCaller$isMale<-FALSE
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
  if (str_ends(outputDir,"/",negate=TRUE))
  {
    scCNVCaller$locPrefix=str_c(scCNVCaller$locPrefix,"/",sep="")
  }
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
  write.table(summaryStats,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,".tsv"),quote=FALSE)
  return(summaryStats)
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
#' @param upperFilterQuantile - exclude cells with reads above certain percentile (to avoid doublets)
#' @export
normalizeMatrixN <- function(inputMatrix,logNorm=FALSE,maxZero=2000,imputeZeros=FALSE,blacklistProp=0.8, priorCount=1,blacklistCutoff=100,dividingFactor=1e6,upperFilterQuantile=0.95)
{
  sc_t<-data.table(t(inputMatrix))
  #sc_t
  cellReads<-transpose(sc_t[,lapply(.SD,sum)],keep.names="Cell")
  pdf(str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_signal_distribution.pdf"),width=6,height=4)
  hist(cellReads$V1, breaks=50,main=scCNVCaller$outPrefix,xlab = "Signal")
  abline(v=quantile(cellReads$V1,upperFilterQuantile),col=c("red"),lty=2)
  dev.off()
  
  readsCells=cellReads[,mean(V1)]
  #apply quantile filter
  sc_t[,cellReads[cellReads[,V1>quantile(cellReads$V1,upperFilterQuantile)]]$Cell:=NULL]
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
  print(sprintf("Blaclist regions: %d",length(blacklistRegions)))
  print(sprintf("Total bins: %d",nrow(sc_pos)))
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
  scData_k<- scData_nc_split %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),list(~ (./ dividingFactor)))
  return(scData_k)
}

#CALCULATE DMs
#' getDoubleMinutes
#'
#' This function allows you to express your love of cats.
#' @keywords double minutes
#' @param minThreshold set threshold of Z score to run changepoint analysis; higher value = less sensitive but faster
#' @export
getDoubleMinutes = function(inputMatrix,targetCell,doPlot=FALSE,penalty_type="SIC",doLosses=FALSE,peakCutoff=5,lossCutoff=-1,minThreshold=4)
{
  medianDensity<-median(inputMatrix[,targetCell])
  if (medianDensity>0){
    tmp<-inputMatrix[,targetCell]
    #tmp[tmp>4e-6]<-medianDensity
    #consider using IQR Here
    tmp<-as.numeric(scale(tmp,center=median(tmp),scale=IQR(tmp)),stringsAsFactors=FALSE)
    #adjusted minseglength 
    #print(max(tmp))
    if (max(tmp)>minThreshold)
    {
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
    else  {
      return(NULL)
    }
  }
  else {
    return(NULL)
  }
}
#' Test overlap
#'
#' Tests overlap between two intervals
#' Internal use
#' @param interval.1 (text interval)
#' @param interval.2 (text interval)
#' @export
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
#' @param inputMatrix normalized input matrix
#' @param minCells minimum cells to retain a DM
#' @param qualityCutoff2 minimum cells to retain DM (filter 2)
#' @param peakCutoff Z score minimum to retain a putative DM/amplification
#' @param doPlots  Plot occasional cells as QC
#' @param imageNumber  Plot every Nth cell (use in conjunction with doPlots = TRUE)
#' @param minThreshold
#' @export
identifyDoubleMinutes<-function(inputMatrix,minCells=100,qualityCutoff2=100,peakCutoff=5,lossCutoff=-1,doPlots=FALSE,imageNumber=1000,logTrans=FALSE,cpgTransform=FALSE,doLosses=FALSE,minThreshold=4)
{
  dm_per_cell<-data.frame(cellName=colnames(inputMatrix %>% select(ends_with(scCNVCaller$cellSuffix))),stringsAsFactors=FALSE)
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
    #NEW:
    scData_1h <- data.frame(scData_k_cpg %>% filter(chrom==currentChrom ) %>% mutate(loc1=paste(chrom,pos,sep="_")) %>% select(-chrom,-pos,-cpg))
  #  print(scData_1h)
    #seleect peak cutoff
    #unlistValues <- unlist((scData_1h %>% select(-loc1)))
    #peakCutoff = quantile(scale(unlistValues,center=median(unlistValues),scale=IQR(unlistValues)),.99)
    
    #scData_1h[,2:ncol(scData_1h)]<-scale(scData_1h[,2:ncol(scData_1h)],center=TRUE,scale=TRUE)
    dm_list=data.frame(dm="null",stringsAsFactors=FALSE)
    if (m>1)
    {
      #rm(dm_per_cell_clean)
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
    pb<-txtProgressBar(min=0,max=nrow(dm_per_cell)-1,style=3)
    for (cellNum in 1:(nrow(dm_per_cell)-1))
    {
      setTxtProgressBar(pb,cellNum)
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
        pdf(str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_cell",cellNum,"_",currentChrom,".pdf"),width=6,height=4)
        temp_dm<-getDoubleMinutes(scData_1h,cellNum,doLosses=doLosses,peakCutoff = peakCutoff,lossCutoff = lossCutoff,doPlot=TRUE,minThreshold=minThreshold)
        dev.off()
        plottableCell=FALSE
      }
      else
      {
        temp_dm<-getDoubleMinutes(scData_1h,cellNum,doLosses=doLosses,peakCutoff = peakCutoff,lossCutoff = lossCutoff,doPlot=FALSE,minThreshold=minThreshold)
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
    close(pb)
    names(dm_per_cell_temp)[2:nrow(dm_list)]<-dm_list$dm[2:nrow(dm_list)]
    #dm_per_cell
    #dm_per_cell %>% rename_at(2:(nrow(dm_list)-1),~ dm_list$dm)
    #test for overlaps - start index 2
    #start one smaller
    
    dm_list_test2=dm_list %>% mutate(newIndex=1)
    
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
      #TODO: fudge indices to empty blanks
      tmp_fact<-factor(dm_list_test2$newIndex)
      levels(tmp_fact)<-seq(1,length(levels(tmp_fact)))
      dm_list_test2$newIndex<-as.numeric(levels(tmp_fact)[tmp_fact])
      #option: can clean up dataframe here and replace NA with FALSE, and then proceed instead of making a giant one up front
      dm_per_cell_temp[is.na(dm_per_cell_temp)]<-FALSE
      #OK that worked, now to amalgamate cell-wide data
      for (j3 in 2:nrow(dm_list_test2))
      {
        #combine and rename columns appropriately
        dm_per_cell_temp[,dm_list_test2$newIndex[j3]]<-(dm_per_cell_temp[,j3] | dm_per_cell_temp[,dm_list_test2$newIndex[j3]])
      }
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
  if (firstRepeatIndex>2)
  {
  dm_per_cell_clean<-dm_per_cell[,1:(firstRepeatIndex-1)]
#  dm_per_cell_clean %>% summarize_if(is.logical,sum)
  
  #can cutoff here too (def 50-75)
  secondThreshold=qualityCutoff2
  keepAlterations <- dm_per_cell_clean %>% summarize_if(is.logical,sum) %>% gather(Chrom,Value) %>% filter(Value>secondThreshold)
  dm_per_cell_clean <- dm_per_cell_clean %>% select(cellName, keepAlterations$Chrom)
  return(dm_per_cell_clean)
  }
  else
  {
    print("No double minutes identified")
  }
  return(NULL)
}

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
#' Scale matrix bins
#'
#' Merges together adjacent windows
#' @param inputMatrix matrix of input stuff
#' @param binSizeRatio ratio to expand by (2 = double bin size)
#' @param inColumn
#' @keywords cats
#' @examples
#' cat_function()
scaleMatrixBins <- function(inputMatrix, binSizeRatio,inColumn)
{
  newMatrix<-inputMatrix %>% mutate(pos_b=floor(as.numeric(!!rlang::sym(inColumn))/(scCNVCaller$binSize*binSizeRatio))*(scCNVCaller$binSize*binSizeRatio))
  newMatrix<-newMatrix %>% select(-inColumn) %>% group_by(pos_b,chrom) %>% summarise_if(is.numeric,list(sum)) %>% arrange(chrom,pos_b)
  
  return(newMatrix)
}

#' collapseChrom3N
#'
#' This function summarises signal across chromosomes and normalises using the CpG content
#' @param inputMatrix: normalised matrix
#' @param summaryFunction: function to use to summarise signal
#' @param logTrans: whether matrix is log-transformed
#' @param binExpand: whether to merge adjacent bins, if so how many (merges if > 1)
#' @param minimumChromValue: cutoff to keep a chromosome 
#' @param tssEnrich multiply CpG by weighting factor (default 1)
#' @param logBase  log base to apply if doing a log weighting of the CpG scaling factor
#' @param minCPG keep only regions with at least CpG score of X (300)
#' @param powVal scaling factor adjustment (cpg^powVal, 0.7-0.75 works well for hg38 and hg19)
#' @keywords CNV
#' @export
collapseChrom3N<-function(inputMatrix,minimumSegments=40,summaryFunction=cutAverage,logTrans=FALSE,binExpand=1,minimumChromValue=2,tssEnrich=5,logBase=2,minCPG=300,powVal=0.73)
{
  #remove centromeric bins
  #these are ordered in the same fashion so should be accurate
  #divide signal by cpg density + 1 (To avoid dividing by zero)
  #remove blacklist regions
  # %>% filter(cpg!=0)
  scData_k_norm <- inputMatrix %>% mutate(chromArm=scCNVCaller$cytoband_data$V4,cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(chrom=str_c(chrom,chromArm)) %>% filter(cpg>minCPG)
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
  chromXName="chrXq"
  chromYName="chrYq"
  x_tmp<-(tail(scCNVCaller$cytoband_data[which(scCNVCaller$cytoband_data$V1=="chrX"),c(1,4)],n=1))
  chromXName<-paste(x_tmp$V1,x_tmp$V4,sep="")
  y_tmp<-(tail(scCNVCaller$cytoband_data[which(scCNVCaller$cytoband_data$V1=="chrY"),c(1,4)],n=1))
  chromYName<-paste(y_tmp$V1,y_tmp$V4,sep="")
  sckn<-sckn[c(chromXName,chromYName),lapply(.SD,quantile,probs=0.8),by="chrom"]
  nrow(sckn)
  total_cpg[chrom %in% c(chromXName,chromYName)]
  sckn[,cpg:=total_cpg[chrom %in% c(chromXName,chromYName)]$V1]
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
  #print(xy_signal)
  diffXY<-xy_signal$chrXq-xy_signal$chrYq
  # print(abs(diffXY))
  
  if (abs(diffXY)<sexCutoff)
  {
    scCNVCaller$isMale=TRUE
  }
  #top part works awesome
  
  #now to convert rest
  sckn<-data.table(scData_k_norm)
  #print(nrow(sckn))
  sckn <- sckn[chromArm!="cen"][blacklist==FALSE]
  #print(nrow(sckn))
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
  #print(median_chrom_signal)
  
  
  #median_chromosome_density[,1:5]
  #]
  #median_chromosome_density<-scData_k_norm %>%
  #  group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) 
  #message(nrow(scData_k_norm))
  #now append 
  #summaryFunction
  
  
  
  #TODO: cpg of 200 and medianNorm cutoffs only work for 1e6 bins; need to change for others
  if (logTrans==TRUE)
  {
    #5.1e-6 old 
    #5.1e-06, totalcpg 130 for 200k
    #THIS ONE
    #transpose(sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL]
    sckn<-sckn[,lapply(.SD,"/",1+log(tssEnrich*cpg^powVal,base=logBase)),by="chrom"][,cpg:=NULL][,lapply(.SD,summaryFunction),by="chrom"]
    #scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+log(cpg,base=2)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    #  select(-cpg) %>% filter(chrom %in% median_chrom_signal$chrom) %>% 
    #  mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>minimumChromValue) 
    #+log(totalcpg,base=2)
    #was 200
  }
  else
  {
    #THIS ONE
    sckn<-sckn[,lapply(.SD,"/",1+tssEnrich*cpg^powVal),by="chrom"][,cpg:=NULL][,lapply(.SD,summaryFunction),by="chrom"]
    #sckn<-sckn[,lapply(.SD,"/",1+tssEnrich*log(cpg,base=logBase)),by="chrom"][,cpg:=NULL][,lapply(.SD,sum),by="chrom"]
    #scData_prechrom<-scData_k_norm %>% mutate_at(vars(ends_with(scCNVCaller$cellSuffix)),funs(./(1+sqrt(cpg)))) %>% group_by(chrom) %>% summarise_if(is.numeric,list(summaryFunction)) %>% 
    #  select(-cpg) %>% filter(chrom %in% median_chrom_signal$chrom) %>%
    #  mutate(medianNorm=median_chrom_signal$Value) %>% filter(medianNorm>minimumChromValue) 
  }
  
  sckn<-sckn[,medianNorm:=median_chrom_signal$V1]
  sckn<-sckn[which(sckn[,medianNorm>minimumChromValue])]
  #return the appropriate matrices
  return(as.data.frame(sckn))
}

#STEP 1:  bin signal per chromosome
#' filterCells
#'
#' This function filters the collapsed chromosome signal density matrix based on predefined parameters (minimum non-zero chromosome segments)
#' @param inputMatrix - collapsed, normalised matrix
#' @param minimumSegments  the minimum segments required to keep a cell
#' @param minDensity  minimum segment density to consider useful (defaults to zero)
#' @param signalSDcut Number of standard deviations average signal of a cell must be from an average per-cell signal to be considered an outlier (to identify putative doublets; default: 2)
#' @keywords filter
#' @keywords CNV
#' @export
#' @examples
#' filterCells(inputMatrix,minimumSegments=41)
filterCells <- function(inputMatrix,minimumSegments=40,minDensity=0,signalSDcut=2)
{
  tmp<-inputMatrix
  initialCells<-ncol(tmp)-1
  cellQuality<-tmp %>% 
    gather(Cell,Density,ends_with(scCNVCaller$cellSuffix)) %>% select(chrom,Cell,Density) %>%
    spread(chrom,Density) 
  #print(head(cellQuality))
  cellQuality <- cellQuality %>% mutate(count=rowSums(.[2:ncol(cellQuality)]>minDensity))
  #print(cellQuality$count)
  barcodesPassing <- (cellQuality %>% filter(count>=minimumSegments))$Cell
  print(sprintf("%d of %d barcodes passed quality check",length(barcodesPassing),initialCells))
  
  #filter scData_prechrom for quality checks
  tmp <- tmp %>% select(c("chrom",barcodesPassing))
  
  #average total signal per cell
  #remove outlier cells
  aveSignal <- tmp %>% summarize_at(vars(ends_with(scCNVCaller$cellSuffix)),list(sum)) %>% gather(Cell,Average)
  #se_aveSignal <- sd(aveSignal$Average)/sqrt(length(aveSignal$Average))
  sd_aveSignal <- sd(aveSignal$Average)
  mean_aveSignal <- mean(aveSignal$Average)
  #sd_aveSignal
  #mean_aveSignal - signalSDcut*sd_aveSignal
  #select cells
  cellsPassingB <- aveSignal %>% filter(between(Average,mean_aveSignal - signalSDcut*sd_aveSignal,mean_aveSignal + signalSDcut*sd_aveSignal)) %>% select(Cell)
  #print(str_c(length(cellsPassingB$Cell)," barcodes passed quality check 2"))
  print(sprintf("%d of %d barcodes passed quality check",length(cellsPassingB$Cell),initialCells))
  #filter scData_prechrom for quality checks
  scData_chrom_filtered <- tmp %>% select(c("chrom",cellsPassingB$Cell)) 
  #average total signal per cell
 #scData_chrom_filtered <- tmp 
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
 # print(median_list)
 # print(iqr_list)
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
#' @param inputMatrix
#' @param median_iqr matrix of medians and IQRs
#' @param useDummyCells use pseudodiploid control
#' @param propDummy pseudodiploid control as proportion of normal cells (0.25 = 25 percent)
#' @param deltaMean minimum difference between clusters
#' @param subsetSize subset of cells to use for initial computation (useful in large samples)
#' @param fakeCellSD SD of dummy cells as proportion of range of normal cells (0.10 is default, 0.10-0.20 works best in most cases)
#' @param summaryFunction summarising function for chromosome arms (default is trimmed average)
#' @param median_iqr matrix of medians and IQRs
#' @param median_iqr matrix of medians and IQRs
#' @param median_iqr matrix of medians and IQRs
#' @export
identifyCNVClusters <- function(inputMatrix, median_iqr, useDummyCells = TRUE,propDummy=0.25, deltaMean=0.03, subsetSize=500, fakeCellSD=0.1, summaryFunction=cutAverage, minDiff=0.25,deltaBIC2=0.25, bicMinimum=0.1,minMix=0.3,uncertaintyCutoff=0.55, mergeCutoff=3,summarySuffix="",shrinky=0,maxClust=4,IQRCutoff=0.25,medianQuantileCutoff=0.3,normalCells=NULL)
{
  pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_cnv_summary",summarySuffix,".pdf"),width=8,height=6)
  
  scData_chrom<-inputMatrix
  median_chrom_signal<-median_iqr[[1]]
  median_chrom_signal_filter <- median_chrom_signal %>% filter(chrom %in% scData_chrom$chrom) %>% filter(chrom!="chrX",chrom!="chrY")
  IQR_chrom_signal<-median_iqr[[2]]
  IQR_chrom_signal_filter <- IQR_chrom_signal %>% filter(chrom %in% median_chrom_signal_filter$chrom)
  #
  if (medianQuantileCutoff!=-1)
  {
     median_chrom_signal_filter$Value<-median(median_chrom_signal_filter$Value)
  }
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
  
  med_iqr2<-list()
  if (medianQuantileCutoff!=-1)
  {
     med_iqr2<-computeCenters(inputMatrix = scData_chrom2,summaryFun = summaryFunction)
  }
  else
  {
    med_iqr2<-computeCenters(inputMatrix = scData_chrom2 %>% select(chrom,normalCells),summaryFun = summaryFunction) 
  }
  median_chrom_signal_filter <- med_iqr2[[1]] %>% filter(Value>0)
  IQR_chrom_signal_filter <- med_iqr2[[2]] %>% filter(chrom %in% median_chrom_signal_filter$chrom)
  #print(median_chrom_signal_filter)
  #clean up again
  #changed from median
  if (medianQuantileCutoff!=-1)
  {
    median_chrom_signal_filter$Value<-quantile(median_chrom_signal_filter$Value,medianQuantileCutoff)
    IQR_chrom_signal_filter$Value<-quantile(IQR_chrom_signal_filter$Value,IQRCutoff)
  }
  #print(median_chrom_signal_filter)
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
  print(ggplot(scData_chrom_spread,aes(chrom,Density))+geom_violin(scale="width") + theme(axis.text.x=element_text(angle=-90,hjust = 0,vjust=0.5),
                                                                                          
                                                                                          axis.text.y=element_text(color="#000000"),
                                                                                          axis.line.x.bottom = element_line(colour="#000000"),
                                                                                          axis.line.y.left =  element_line(colour="#000000"),
                                                                                          panel.background = element_rect(fill="#ffffff"),
                                                                                          plot.background = element_rect(color="#ffffff",fill="#ffffff")) +
          xlab("Segment")) 
  
  
  
  #initialize variables for clustering
  cell_ids=unique(scData_chrom_spread$Cell)
  chrom_clusters<-data.frame(Chrom=median_chrom_signal_filter$chrom,V1=0,V2=0,V3=0,V4=0,V5=0,V6=0,stringsAsFactors=FALSE)
  cell_assignments<-data.frame(Cells=cell_ids,stringsAsFactors = FALSE)
  
  #initialization
  set.seed(0)
  defSubset=sample(1:length(unique(scData_chrom_spread$Cell)),subsetSize)
  print(ggplot(scData_chrom_spread %>% filter(Cell %in% unique(scData_chrom_spread$Cell)[defSubset]),aes(chrom,Density))+geom_violin(scale="width") + theme(axis.text.x=element_text(angle=-90,hjust = 0,vjust=0.5),
                                                                                          
                                                                                          axis.text.y=element_text(color="#000000"),
                                                                                          axis.line.x.bottom = element_line(colour="#000000"),
                                                                                          axis.line.y.left =  element_line(colour="#000000"),
                                                                                          panel.background = element_rect(fill="#ffffff"),
                                                                                          plot.background = element_rect(color="#ffffff",fill="#ffffff")) +
          xlab("Segment")) 
  
  #can cut down if the mixing proportion minimum is less than X% (e.g. 20)
  for (m5 in 1:(length(median_chrom_signal_filter$chrom)))
  {
    set.seed(0)
    
    #hc doesn't wo
    dens1<-scData_chrom_spread %>% filter(chrom==median_chrom_signal_filter$chrom[m5]) %>% select(Density)
    initParams=list(noise=TRUE,subset=defSubset) #hcPairs=randomPairs(dens1$Density,modelName = "V")) #,noise=sample(1:length(dens1$Density),10)) 
    #message(median_chrom_signal_filter$chrom[m5])
    #
    print(median_chrom_signal_filter$chrom[m5])
    fit = Mclust(dens1$Density,G=1:(maxClust-1),model=c("V"), prior = priorControl(shrinkage=shrinky), initialization = initParams) 
    #print(plot(fit,what="classification")+title(median_chrom_signal_filter$chrom[m5]))
    #print(summary(fit))
    #retype if greater than two components
    if (fit$G>2)
    {
      #print((fit$BIC[fit$G]-fit$BIC[fit$G-1]))
      if ((fit$BIC[fit$G]-fit$BIC[fit$G-1])<deltaBIC2)
      {
        fit = Mclust(dens1$Density,G=1:2,prior = priorControl(),model="V",initialization = initParams)
      }
      else
      {
        range_lists<-list()
        #now do the matching for 3
        overlapping_clusters<-list()
        clust_list1<-c()
        for (c1 in 1:fit$G)
        {
          if (length(which(fit$classification==c1))>0)
          {
            range_lists[[c1]]<-range(fit$data[fit$classification==c1])
          #ignore empty clusters
            clust_list1<-append(clust_list1,c1)
          }
        }
        for (c1 in clust_list1)
        {
          for (c2 in clust_list1)
          {
            if (c1!=c2)
            {
              #print(range_lists[[c1]])
              #print(range_lists[[c2]])
              isOverlap = (range_lists[[c1]][1]<range_lists[[c2]][1]) & (range_lists[[c1]][2]>range_lists[[c2]][2])
              #print(isOverlap)
              if (isOverlap==TRUE)
              {
                #with c1 being the mama cluster
                overlapping_clusters<-c(overlapping_clusters,list(c(c1,c2)))
              }
            }
          }
        }
    #    print(range_lists)
        #process overlaps
    #    print(overlapping_clusters)
        if (length(overlapping_clusters)>0)
        {
        for (k0 in 1:length(overlapping_clusters))
        {
          overlap_lists<-overlapping_clusters[[k0]]
          if (overlap_lists[1]!=overlap_lists[2])
          {
          #print(overlap_lists)
          #test difference in means
          diff_mean<-fit$parameters$mean[overlap_lists[2]]-fit$parameters$mean[overlap_lists[1]]
          #print(diff_mean)
          if (abs(diff_mean) > deltaMean)
          {
            #test quantiles of larger list
            #if meets criteria then divide as previously and recompute means
            if (diff_mean<0)
            {
              #transfer small 2 to cluster 1
              fit$classification[fit$classification==overlap_lists[1] & fit$data<=fit$parameters$mean[overlap_lists[2]]]<-overlap_lists[2]
    #          message(sprintf("transfer %d into %d",overlap_lists[1],overlap_lists[2]))
   #           print(range_lists[[overlap_lists[1]]])
  #            print(range_lists[[overlap_lists[2]]])
              
            }
            else
            {
     #         message(sprintf("transfer %d into %d",overlap_lists[1],overlap_lists[2]))
   #          print(range_lists[[overlap_lists[1]]])
  #            print(range_lists[[overlap_lists[2]]])
              #2 is shifted right of one so right fill the cluster
              fit$classification[fit$classification==overlap_lists[1] & fit$data>=fit$parameters$mean[overlap_lists[2]]]<-overlap_lists[2]
              #relabel others
            }
           #  
            #recompute means
            fit$parameters$mean[overlap_lists[1]]<-mean(fit$data[fit$classification==overlap_lists[1]])
            fit$parameters$mean[overlap_lists[2]]<-mean(fit$data[fit$classification==overlap_lists[2]])
          }
          else
          {
            #merge small into big
            message(sprintf("transfer %d into %d",overlap_lists[2],overlap_lists[1]))
         #   print(range_lists[[overlap_lists[1]]])
        #    print(range_lists[[overlap_lists[2]]])
            fit$classification[fit$classification==overlap_lists[2]]<-overlap_lists[1]
            for (k1 in 1:length(overlapping_clusters))
                  {
                    print(sprintf("overlap %d; list %d",overlapping_clusters[[k1]][1],overlap_lists[1]))
                    if (overlapping_clusters[[k1]][1]==overlap_lists[2])
                    {
                      #message("relabel")
                      overlapping_clusters[[k1]][1]=overlap_lists[1]
                    }
                  }
          }
          }
          
        }
        }
        #check if two 
        if (length(which(fit$classification==2))==0)
        {
          #push pop 3 into 2
          fit$classification[fit$classification==3]<-2
          fit$parameters$mean[2]<-fit$parameters$mean[3]
        }
        #range1<-range(fit$data[fit$classification==1])
        #range2<-range(fit$data[fit$classification==2])
        #range3<-range(fit$data[fit$classification==3])
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
      #print(str_c("deltaBIC",deltaBIC))
      if ((deltaBIC < bicMinimum) | (minProp < propCutoff))
      {
        #message(median_chrom_signal_filter$chrom[m5])
        
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
       # print(varChange)
      #  print(diffMean)
      #  print("Range")
      #  print(range(fit$data[fit$classification==1]))
      #  print(range(fit$data[fit$classification==2]))
        #if the larger slops greater than smaller then we need to fix it
       # if ((diffMean < deltaMean)) # | ((varChange > 4) & (diffMean < deltaMean * 4)))
       # {
      #    message(median_chrom_signal_filter$chrom[m5])
          
      #    fit = Mclust(dens1$Density,G=1,model="V",prior = priorControl(),initialization = initParams)
      #  }
       # else
       # {
          if ((varChange > 2.0) & (diffMean<mergeCutoff))
          {
      #      message(str_c("cutting dist", median_chrom_signal_filter$chrom[m5]))
      #      print(fit$parameters$mean)
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
            #recompute means
            fit$parameters$mean[1]<-mean(fit$data[fit$classification==1])
            fit$parameters$mean[2]<-mean(fit$data[fit$classification==2])
          }
     #   }
      }
    }
    #test uncertainty
    #print(fit$uncertainty)
    fit$classification[fit$uncertainty>uncertaintyCutoff]<-0
    #print plot
    print(plot(fit,what="classification")+title(median_chrom_signal_filter$chrom[m5]))
    print(summary(fit))
    #can do adaptations here to check for overfitting
    chrom_clusters[chrom_clusters$Chrom==median_chrom_signal_filter$chrom[m5],2:(2+length(fit$parameters$mean)-1)]<-fit$parameters$mean
    cell_assignments <- cbind.data.frame(cell_assignments,fit$classification,stringsAsFactors=FALSE)
    
    names(cell_assignments)[m5+1]=median_chrom_signal_filter$chrom[m5]
    #TODO: compute LRT statistic
  }
  names(cell_assignments)[2:(length(median_chrom_signal_filter$chrom)+1)]<-median_chrom_signal_filter$chrom
  #try to assign groups based on imputed clusters
 # chrom_clusters
  dev.off()
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
#' @param initialResultList output list generated by identifyCNVClusters
#' @param medianIQR median and IQR list of sample
#' @param maxClust maximum # of clusters fed into identifyCNVClusters
#' @param minDiff minimum Z-score difference between clusters
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
    #print(clust_list[2:(2+numClusts-1)])
    clust_list[2:(2+numClusts-1)]<-activeClusters[order(activeClusters)]
  #  print(clust_list[2:(2+numClusts-1)])
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
    #newOrder
    #TODO: what about zeros
    levels(cell_assignments[,m5+1])<-newOrder
    #tmp_fact<-factor(dm_list_test2$newIndex)
    #
    cell_assignments[,m5+1]<-as.numeric(levels(cell_assignments[,m5+1])[cell_assignments[,m5+1]])
  #  print("lags")
    minDiffA=minDiff-0.00001
    while (minDiffA>0 & minDiffA<minDiff)
    {
      tmpa<-data.frame(init=seq(2,1+numClusts))
      tmpa<-tmpa %>% mutate(init2=lag(init,default=2)) %>% mutate(clust_diff=(as.double(chrom_clusters_final[m5,init])-as.double(chrom_clusters_final[m5,init2]))) %>% filter(clust_diff!=0) %>% arrange(clust_diff)
      if (nrow(tmpa)==0)
      {
        minDiffA=0
      }
      else
      {
        minDiffA=min(tmpa$clust_diff)
        ##print("min diff")
        ##print(minDiffA)
        ##print(tmpa)
        #now iterate
        for (j1 in 1:nrow(tmpa))
        {
          #print(tmpa$clust_diff[j1])
          if (tmpa$clust_diff[j1]<minDiff)
          {
            #merge these indices
            #print(chrom_clusters_final[m5,])
            #print(str(chrom_clusters_final[m5,]))
            #print(tmpa$init[j1])
            #print(chrom_clusters_final[m5,tmpa$init[j1]])
            #weighted average
            number_a<-length(which(cell_assignments[,m5+1]==(tmpa$init[j1]-1)))
            number_b<-length(which(cell_assignments[,m5+1]==(tmpa$init2[j1]-1)))
           # print(number_a)
          #  print(number_b)
          #  print(range(as.numeric(cell_assignments[,m5+1])))
            if ((number_a + number_b)!=0)
            {
              
            new_mean<-(chrom_clusters_final[m5,tmpa$init[j1]]*number_a+chrom_clusters_final[m5,tmpa$init2[j1]]*number_b)/(number_a+number_b)
            #print(new_mean)
            chrom_clusters_final[m5,tmpa$init[j1]]<-new_mean
            chrom_clusters_final[m5,tmpa$init2[j1]]<-new_mean
            #       message(str_c(chrom_clusters_final[m5,i1],lastValue,i1,sep=","))
            cell_assignments[cell_assignments[,m5+1]==(tmpa$init2[j1]-1),m5+1]=tmpa$init[j1]-1
            break()
            }
            else
            {
              #cluster is empty thus not assignable
              minDiffA = 0
            }
          }
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
    uniq_clusts<-unique(as.numeric(chrom_clusters_final[m5,2:ncol(chrom_clusters_final)]))
    levels(cell_assignments[,m5+1])<-seq(minCellAssign,length(uniq_clusts))
    #tmp_fact<-factor(dm_list_test2$newIndex)
    
    
    chrom_clusters_final[m5,2:ncol(chrom_clusters_final)]<-0
    chrom_clusters_final[m5,2:(1+length(uniq_clusts))]<-uniq_clusts
    #
    # message(levels(cell_assignments[,m5+1]))
    cell_assignments[,m5+1]<-as.numeric(levels(cell_assignments[,m5+1])[cell_assignments[,m5+1]])
    #TODO: merge cluster values and medians into a list which we can merge (e.g. chr0_0,0.5))
  }
 # chrom_clusters_final
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
  
  for (suspCNV in possCNV$Chrom)
  {
    #message(possCNV)
    clusterList<-unique(chrom_clusters_final %>% filter(Chrom==suspCNV) %>% select(-Chrom))
    #print(clusterList)
    if (possCNV$maxClust[possCNV$Chrom==suspCNV]==2)
    {
      abn<-min(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]]))
      #first item is normal one
      #print(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]]))
      #print(abn)
      #print(which(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]])==abn))
      abn_index<-which(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]])==abn)
      #print(order(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]])))
      #compare to median autosomal value
      #if value is less than autosomal, suspect a loss, otherwise suspect a gain to be abnormal
     # possCNV
    #  print(as.numeric(clusterList[abn_index]))
      possCNV$normCluster[possCNV$Chrom==suspCNV] =abn_index
     # print(possCNV)
    }
    else
    {
      #assume middle cluster is normal
      possCNV$normCluster[possCNV$Chrom==suspCNV] =2
     # print(possCNV)
    }
  }
  #merge cell assignments and annotate
  v1<-rownames_to_column(cell_assignments) %>% filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% filter(Chrom %in% possCNV$Chrom)
  v2<-left_join(v1,possCNV,by="Chrom") %>% mutate(clust_text=if_else(clust>normCluster,3,if_else(clust<normCluster,1,2))) %>% mutate(clust_text=if_else(clust==0,2,clust_text))
  v2 <- v2 %>% select(rowname,Chrom,clust_text) %>% spread(Chrom,clust_text) %>% arrange(rowname)
  if(saveOutput==TRUE)
  {
    write.table(x=v2,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv_anno_binary.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  colnames(possCNV)[1]<-"chrom"
  
  final_cnv_binarized<-inner_join(consensus_CNV_clusters %>% gather(chrom,clust,2:ncol(consensus_CNV_clusters)),possCNV,by="chrom") %>% mutate(clust=if_else(clust==normCluster,0,1)) %>% select(rowname,chrom,clust) %>% spread(chrom,clust)
  return(list(v2,possCNV,final_cnv_binarized))
}
#'
#' annotateCNV4
#' 
#' Annotates the filtered CNV calls and saves the results, and estimates absolute copy numbers
#' @param cnvResults   output of filtered CNVs (set of lists)
#' @param saveOutput    save output to a file on disk
#' @param maxClust2    maximum # of clusters
#' @param outputSuffix   output suffix to save files
#' @param sdCNV   sd for CNV copy number estimation
#' @param filterResults filter results to remove likely unaltered chromosomes
#' @param filterRange  copy number range minimum
#' @param minAlteredCells filtering parameter - remove all alterations with a smallest group less than X cells (default is 40)
#' @keywords CNV
#' @keywords output
#' @export
#' @examples
#' annotateCNV4(cleanCNV,saveOutput=TRUE,outputSuffix="_nolog",sdCNV=0.5)
annotateCNV4 <- function(cnvResults,saveOutput=TRUE,maxClust2=4,outputSuffix="_1",sdCNV=0.6,filterResults=TRUE,filterRange=0.8,minAlteredCells=40)
{
  cell_assignments<-cnvResults[[1]]
  cell_assignments<-column_to_rownames(cell_assignments,var = "Cells")
  chrom_clusters_final<-cnvResults[[2]]
  #colnames(chrom_clusters_final)=c("Chrom",seq(1:(maxClust2-2)),"0")
  shift_val=0
  if (length(which(chrom_clusters_final$V2==0))>1)
  {
    shift_val = mean(chrom_clusters_final$V1[which(chrom_clusters_final$V2==0)])
  }
  #print(shift_val)
  #identify possible CNVs - can use bins to test
  possCNV <- cell_assignments %>% summarise_if(is.numeric,list(max)) %>% gather(Chrom,maxClust) %>% dplyr::filter(maxClust>1)
  possCNV <- possCNV %>% mutate(Type="Unknown")
 # print(cell_assignments %>% summarise_if(is.numeric,list(max)))
  #%>% gather(Chrom,maxClust))
 # print(possCNV)
  #label by Z scores
  chrom_clusters_final<-chrom_clusters_final %>% gather("cluster","zscore",starts_with("V")) %>% mutate(cluster=str_remove(cluster,"V")) 
  #print(chrom_clusters_final)
  #print(which(is.nan(chrom_clusters_final$zscore)))
  chrom_clusters_final$zscore<-sapply(sapply(chrom_clusters_final$zscore-shift_val,pnorm,log.p=TRUE),qnorm,2,sdCNV,log.p=TRUE)
  #replace each column
  #print(chrom_clusters_final)
  trimmedCNV<-vector()
  thresholdVal=filterRange
  for (chrom in possCNV$Chrom)
  {
  #  print(chrom)
    zscores<-chrom_clusters_final %>% dplyr::filter(Chrom==chrom) %>% dplyr::select(zscore)
  #  print(zscores)
    #print(head(cell_assignments[chrom]))
    #print(as.character(factor(cell_assignments[,chrom],levels=seq(from=0,to=6),labels = c("2",zscores$zscore))))
    cell_assignments[,chrom]<-as.numeric(as.character(factor(cell_assignments[,chrom],levels=seq(from=0,to=6),labels = c("2",zscores$zscore))))
   # print(zscores)
 #  print(chrom)
#    print(range(zscores))
    if (diff(range(zscores))>=thresholdVal)
    {
 #     print(chrom)
      trimmedCNV<-append(trimmedCNV,chrom)
    }
  }
  #print(trimmedCNV)
  #output CNV list by cluster
  #cell_assignments %>% filter()
  #WITHOUT BLANKS
  #consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  
  #WITH BLANKS
  if (!filterResults)
  {
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  else
  {
    #identify chromosomes not used
    unused = cnvResults[[1]] %>% dplyr::filter(str_detect(Cells,"X",negate=TRUE)) %>% gather(chrom,cluster,starts_with("chr")) %>% group_by(chrom,cluster) %>% filter(cluster!=0) %>% summarise(num=n()) %>% group_by(chrom) %>% summarise(min=min(num)) %>% dplyr::filter(min<minAlteredCells)
   # print(unused)
    #print(trimmedCNV)
    #print(trimmedCNV %in% unused$chrom)
    trimmedCNV<-trimmedCNV[!(trimmedCNV %in% unused$chrom)]
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% trimmedCNV) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  if (saveOutput==TRUE)
  {
    write.table(x=consensus_CNV_clusters,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv_scores.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  return(list(chrom_clusters_final,possCNV,consensus_CNV_clusters))
}
#' annotateCNV4B
#' 
#' Annotates the filtered CNV calls and saves the results, and estimates absolute copy numbers
#' unlike the regular version, this uses predicted normal cells to set the normal values to 2
#' @param cnvResults   output of filtered CNVs (set of lists)
#' @param expectedNormals   list of normal barcodes
#' @param saveOutput    save output to a file on disk
#' @param maxClust2    maximum # of clusters
#' @param outputSuffix   output suffix to save files
#' @param sdCNV   sd for CNV copy number estimation
#' @param filterResults filter results to remove likely unaltered chromosomes
#' @param filterRange  copy number range minimum
#' @param minAlteredCellProp filtering parameter - proportion of normal cell population that should be unaltered, default is 0.75
#' @keywords CNV
#' @keywords output
#' @export
#' @examples
annotateCNV4B <- function(cnvResults,expectedNormals,saveOutput=TRUE,maxClust2=4,outputSuffix="_1",sdCNV=0.6,filterResults=TRUE,filterRange=0.8,minAlteredCellProp=0.75)
{
  cell_assignments<-cnvResults[[1]]
  normal_clusters<-cnvResults[[1]] %>% filter(Cells %in% expectedNormals) %>% summarise_if(is.numeric,mean) %>% mutate_all(round)
  cell_assignments<-column_to_rownames(cell_assignments,var = "Cells")
  chrom_clusters_final<-cnvResults[[2]]
  chrom_clusters_final<-chrom_clusters_final %>% mutate(norm=t(normal_clusters)) %>% mutate(zoffset=if_else(norm=="1",V1,V2))
  #no offset for X or Y
  chrom_clusters_final$zoffset[chrom_clusters_final$Chrom %in% c("chrXp","chrXq","chrYp","chrYq")]<-0
  #colnames(chrom_clusters_final)=c("Chrom",seq(1:(maxClust2-2)),"0")
  shift_val=0
  if (length(which(chrom_clusters_final$V2==0))>1)
  {
    shift_val = mean(chrom_clusters_final$V1[which(chrom_clusters_final$V2==0)])
  }
  #print(shift_val)
  #identify possible CNVs - can use bins to test
  possCNV <- cell_assignments %>% summarise_if(is.numeric,list(max)) %>% gather(Chrom,maxClust) %>% dplyr::filter(maxClust>1)
  possCNV <- possCNV %>% mutate(Type="Unknown")
  # print(cell_assignments %>% summarise_if(is.numeric,list(max)))
  #%>% gather(Chrom,maxClust))
  # print(possCNV)
  #label by Z scores
  chrom_clusters_final<-chrom_clusters_final %>% gather("cluster","zscore",starts_with("V")) %>% mutate(cluster=str_remove(cluster,"V")) 
  #print(chrom_clusters_final)
  #print(which(is.nan(chrom_clusters_final$zscore)))
  chrom_clusters_final$zscore<-sapply(sapply(chrom_clusters_final$zscore-chrom_clusters_final$zoffset,pnorm,log.p=TRUE),qnorm,2,sdCNV,log.p=TRUE)
  #replace each column
  #print(chrom_clusters_final)
  trimmedCNV<-vector()
  thresholdVal=filterRange
  for (chrom in possCNV$Chrom)
  {
    #  print(chrom)
    zscores<-chrom_clusters_final %>% dplyr::filter(Chrom==chrom) %>% dplyr::select(zscore)
    #  print(zscores)
    #print(head(cell_assignments[chrom]))
    #print(as.character(factor(cell_assignments[,chrom],levels=seq(from=0,to=6),labels = c("2",zscores$zscore))))
    cell_assignments[,chrom]<-as.numeric(as.character(factor(cell_assignments[,chrom],levels=seq(from=0,to=6),labels = c("2",zscores$zscore))))
    # print(zscores)
    #  print(chrom)
    #    print(range(zscores))
    if (diff(range(zscores))>=thresholdVal)
    {
      #     print(chrom)
      trimmedCNV<-append(trimmedCNV,chrom)
    }
  }
  #print(trimmedCNV)
  #output CNV list by cluster
  #cell_assignments %>% filter()
  #WITHOUT BLANKS
  #consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  
  #WITH BLANKS
  if (!filterResults)
  {
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  else
  {
    #identify chromosomes not used
    #filter those that have fewer altered than expected # of normal cells
    unused = cnvResults[[1]] %>% dplyr::filter(str_detect(Cells,"X",negate=TRUE)) %>% gather(chrom,cluster,starts_with("chr")) %>% group_by(chrom,cluster) %>% filter(cluster!=0) %>% summarise(num=n()) %>% group_by(chrom) %>% summarise(min=min(num)) %>% dplyr::filter(min<(minAlteredCellProp*length(expectedNormals)))
    # print(unused)
    #print(trimmedCNV)
    #print(trimmedCNV %in% unused$chrom)
    trimmedCNV<-trimmedCNV[!(trimmedCNV %in% unused$chrom)]
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% trimmedCNV) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  if (saveOutput==TRUE)
  {
    write.table(x=consensus_CNV_clusters,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv_scores.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  return(list(chrom_clusters_final,possCNV,consensus_CNV_clusters))
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
  print(ggplot(inputMatrix  %>% gather(Cell,Density,ends_with(scCNVCaller$cellSuffix)),aes(chrom,Density))+geom_violin(scale="width",trim=FALSE) +
          theme(axis.text.x=element_text(angle=-90,vjust = 0.5,hjust = 0, color = "#000000"),
                axis.text.y=element_text(color="#000000"),
                axis.line.x.bottom = element_line(colour="#000000"),
                axis.line.y.left =  element_line(colour="#000000"),
                panel.background = element_rect(fill="#ffffff"),
                plot.background = element_rect(color="#ffffff",fill="#ffffff")) +
          xlab("Segment"))
  dev.off()
}
#LOSS CNV
#' getLOHRegions
#'
#' Calls putative regions of LOH
#' @param inputMatrix - normalised input fragment count matrix
#' @param lossCutoff - minimum Z score for losses (default is -0.25)
#' @param uncertaintyCutLoss - maximum uncertainty to accept a cell assignment
#' @param diffThreshold - difference threshold between bins
#' @param minLength - minimum length of a LOH region
#' @param minSeg - minimum number of bins for changepoint analysis
#' @param lossCutoffCells - minimum number of cells to accept a LOH (default is 100)
#' @param targetFun - target function to identify LOH regions
#' @param signalBoost - value to add prior to logtransforming signal
#' @param lossCutoffReads minimum number of cells showing alteration to keep in final list (filter)
#' @param quantileLimit quantile of signal per cell population to use to identify losses in changepoints (default is 0.3)
#' @keywords LOH
#' @keywords CNV
#' @export
getLOHRegions <- function(inputMatrixIn,lossCutoff=(-0.25), uncertaintyCutLoss=0.5, diffThreshold=0.9, minLength=3e6, minSeg=3, lossCutoffCells=100,targetFun=IQR,quantileLimit=0.3,cpgCutoff=0,meanThreshold=4,dummyQuantile=0.5,dummyPercentile=0.2,dummySd=0.2)
{
  #c("E","V")
  pdf(str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_LOH.pdf"),width=6,height=4)
  inputMatrix<-data.table(inputMatrixIn)
  inputMatrix<-inputMatrix[,cpg:=scCNVCaller$cpg_data$cpg_density]
  inputMatrix<-inputMatrix[,arm:=scCNVCaller$cytoband_data$V4]
  inputMatrix<-inputMatrix[which(inputMatrix[,arm!="cen" & blacklist==0 & cpg>cpgCutoff & chrom!="chrY"])]
  inputMatrix<-inputMatrix[,c("raw_medians","blacklist","cpg","arm"):=NULL]
  #inputMatrix<-input
  #%>% mutate(cpg=cpg_data$cpg_density) %>% mutate(arm=cytoband_data$V4) %>% filter(arm!="cen", blacklist==0, cpg>0) %>% select(-arm,-blacklist) #%>%  #cpg+
  #check if raw_medians
  #mutate
  #view(inputMatrix$cpg)
  #inputMatrix<-inputMatrix %>% mutate_at(vars(ends_with("-1")),funs(./(log(cpg,base=2))))
  #print(nrow(inputMatrix))
  chromList<-unique(inputMatrix$chrom)
  # chromList<-c("chr9","chr10")
  #print(chromList)
  dm_per_cell_vals<-data.frame(cellName=colnames(inputMatrix  %>% select(ends_with(scCNVCaller$cellSuffix))),stringsAsFactors=FALSE)
  alteration_list=c()
  last_coords=c(0,0)
  alteration_delta=c()
  sliceList<-function(fit)
  {
    #reassign items
    #assume there are only two
    #recompute means
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
    #recompute means
     mean1a<-mean(fit$data[fit$classification==1])
     mean2b<-mean(fit$data[fit$classification==2])
    
    
    return(fit$classification)
  }
  #may not need IQR, may work with mean/median instead
  for (targetChrom in chromList)
  {
    inputMatrixK<-inputMatrix[which(inputMatrix[,chrom==targetChrom])]
    tmp_coords<-inputMatrixK[,c("chrom","pos")] 
    # print(tail(colnames(inputMatrixK)))
    #print(transpose(inputMatrixK[,chrom:=NULL],make.names = "pos"))
    #print(str(transpose(inputMatrixK[,chrom:=NULL],make.names = "pos")))
    #print(tail(colnames((transpose(inputMatrixK[,chrom:=NULL],make.names = "pos")))))
    #message(targetChrom)
    IQRv<-transpose(inputMatrixK[,chrom:=NULL],make.names = "pos")[,lapply(.SD,quantile,quantileLimit,na.rm=TRUE)]
    # print(IQRv)
    #IQRv<-apply(inputMatrix  %>% filter(chrom==targetChrom) %>% select(scCNVCaller$cellSuffix),1,quantile,0.3,na.rm=TRUE)
    IQRs<-scale(t(IQRv),center=TRUE,scale=TRUE)
    if (is.nan(as.vector(IQRs)[1])) {
      next
    }
    #expectedSignal<-inputMatrix %>% filter(chrom==targetChrom) %>% select(-cpg,-pos,-chrom) %>% summarise_if(is.numeric,median) %>% gather(barcode,value) %>% summarise_if(is.numeric,mean)
    #[,chrom:=NULL]
    expectedSignalN<-transpose(inputMatrixK[,lapply(.SD,max),by="pos"],make.names="pos")[,lapply(.SD,mean)]
    #print(as.vector(IQRs))
    cm<-cpt.meanvar(data=as.vector(IQRs),test.stat="Normal", penalty="AIC",method = "PELT",minseglen = minSeg)
    cptlist<-t(rbind(cm@param.est$mean,cm@cpts))
    colnames(cptlist)<-c("Mean","Point")
    cptlist<-as_tibble(cptlist) %>% mutate(Diff = Point - lag(Point))
    plot(cm@data.set,xlab="Chromosome bin",ylim=c(-5,5),ylab="Z-score")
    p1<-as_tibble(cptlist) %>% mutate(Start=lag(Point))
    p1$Start[1]<-0
    p1<-p1 %>% select(Start,Mean,Point,Mean)
    for (a in 1:nrow(p1))
    {
      #print(segments(p1$Start[a],p1$Mean[a],p1$Point[a],p1$Mean[a],col="red"))
      segments(p1$Start[a],p1$Mean[a],p1$Point[a],p1$Mean[a],col="red")
    }
    
    #print(cptlist)
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
      #%>% filter(chrom==targetChrom) %>% select(chrom,pos)
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
        set.seed(12354)
        last_coords[1]<-d_lossb$startCoord[i]
        last_coords[2]<-d_lossb$endCoord[i]
        alterationName<-str_c(targetChrom,sprintf(fmt="%d",last_coords[1]),sprintf(fmt="%d",last_coords[2]),sep="_")
        # message(seg_length)
        
        posList=seq(from=last_coords[1],to=last_coords[2],by = 1e6)
        posList<-sprintf(fmt="%d",posList)
        #posList
        #or max
        #print(head(expectedSignal,n=1))
        t2d<-inputMatrixK[which(inputMatrixK[,pos %in% posList])][,lapply(.SD,max),.SDcols=-c("pos")] # %>% select(-cpg,-pos,-chrom) %>% summarise_if(is.numeric,max) # %>% select(-cpg)
        #t2d<-t2 %>% gather(Cell,Value)
        # print(-1*log2(t2d+signalBoost))
        #print(str(-1*log2(t2d+signalBoost)))
        #generate random vals for zeros
        #        print(expectedSignalN)
        deadVal=unlist(quantile(expectedSignalN,dummyQuantile))
        #print(deadVal)
        #print(str(deadVal))
        # print(head(t2d,n=10))
        if (length(t2d)==length(which(t2d==0)))
        {
          next
        }
        #  print(which(t2d==0))
        #  print(length(which(t2d==0)))
        #print(deadVal)
        t2d<-t(t2d)
        # print(head(t2d))
        #ignore the zeros, just make fake cells here, since we generate this already
        
        
        nCells=nrow(t2d)
        valsToFill<-data.frame(vals=rnorm(n=floor(dummyPercentile * nCells),mean=deadVal,sd = dummySd*deadVal))
        rownames(valsToFill)<-paste("X",seq(from=1,to=floor(dummyPercentile * nCells)),"-1",sep="")
        #print(valsToFill)
        # print(valsToFill)
        #t2d[which(t2d==0)]<-valsToFill
        #print(str(t2d))
        t2d<-rbind(t2d,as.matrix(valsToFill))
        #print(tail(t2d))
        #print(tail(rbind(t2d,valsToFill$vals)))
        #  print(head(t2d,n=10))
        t2d_trans<-as.vector(t(sqrt(t2d+deadVal)))
        
        #  print(as.vector(t(-1*log2(t2d+signalBoost))))
        #print(str(as.vector(-1*log2(t2d+signalBoost))))
        hist(t2d_trans,main=alterationName,breaks = 50)
        # print(inputMatrix %>% filter(chrom==targetChrom,pos %in% posList) %>% select(-cpg))
        # print(targetChrom)
        # print(posList)
        fit1<-Mclust(t2d_trans,G=2:3,modelNames=c("V"),initialization = list(noise=TRUE))
        if (length(fit1)==0)
        {
          #message("broken, retry")
          fit1<-Mclust(t2d_trans,G=1:2,modelNames=c("V"),initialization = list(noise=TRUE))
          if (length(fit1)==0)
          {
            #message("broken 2")
            next()
          }
        }
        #message(fit1$G)
        #plot(fit1)
        #if (!is.na(fit1))
        #{
        #print(plot(fit1))
        if (fit1$G>2)
        {
          #collapse clusters
          #which diff is bigger
          #print(summary(fit1))
          plot(fit1,what="classification")
          diff_2=abs(fit1$parameters$mean[3]-fit1$parameters$mean[2])
          diff_1=abs(fit1$parameters$mean[2]-fit1$parameters$mean[1])
          clust1_mean<-mean(fit1$parameters$mean[1:2])
          clust2_mean<-fit1$parameters$mean[3]
          if (diff_2>2.5*diff_1)
          {
            fit1$classification[fit1$classification==2]<-1
            fit1$classification[fit1$classification==3]<-2
          }
          else
          {
            #fit1$classification[fit1$classification==2]<-1
            fit1$classification[fit1$classification==3]<-2
            clust1_mean<-fit1$parameters$mean[1]
            clust2_mean<-mean(fit1$parameters$mean[2:3])
          }
          #      print(clust2_mean)
          #deadVal=unlist(quantile(expectedSignalN,0.3))
          dVal2<-unlist(quantile(expectedSignalN,0.6))
          plot(fit1,what="classification")
          #recompute means
          fit1$parameters$mean[1]<-mean(fit1$data[fit1$classification==1])
          fit1$parameters$mean[2]<-mean(fit1$data[fit1$classification==2])
          fit1$classification<-sliceList(fit1)
          plot(fit1,what="classification")
          #  print(dVal2)
          #signalDiff=2^(-1*clust2_mean)-dVal2
          signalDiff=dVal2-clust1_mean^2
          #print(signalDiff)
          #print(fit1$parameters$mean)
          
          #print(alterationName)
          delta_mean = clust2_mean-clust1_mean
          #print(delta_mean)
          #print(clust2_mean-clust1_mean)
          #      delta_mean
          #      message(delta_mean)
          if (abs(delta_mean)>diffThreshold)
          {
            alteration_list<-cbind(alteration_list,alterationName)
            alteration_delta<-cbind(alteration_delta,signalDiff)
            fit1$classification[fit1$uncertainty>uncertaintyCutLoss]<-0
            #message(fit1$parameters$mean)
            dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification[1:nCells])
          }
          #      #message(fit1$parameters)
        }
        if (fit1$G==2)
        {
          plot(fit1,what="classification")
          
          #check for and collapse cluster assignments
          #cells in 2 with val < 1 and vice versa
          #fit1$classification[fit1$data<fit1$parameters$mean[1] & fit1$classification==2]<-1
          #fit1$classification[fit1$data>fit1$parameters$mean[2] & fit1$classification==1]<-2
          fit1$classification<-sliceList(fit1)
          plot(fit1,what="classification")
          delta_mean = fit1$parameters$mean[2]-fit1$parameters$mean[1]
          #message(str_c(delta_mean,"MEAN"))
          #       print(delta_mean)
          dVal2<-unlist(quantile(expectedSignalN,0.6))
          
          signalDiff=dVal2-fit1$parameters$mean[1]^2
          #print(signalDiff)
          #print(delta_mean)
          if (abs(delta_mean)>diffThreshold)
          {
            alteration_list<-cbind(alteration_list,alterationName)
            alteration_delta<-cbind(alteration_delta,signalDiff)
            #add delta mean
            #message(fit1$parameters$mean)
            #uncertainty less than 0.2
            fit1$classification[fit1$uncertainty>uncertaintyCutLoss]<-0
            dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification[1:nCells])
          }
          # }
          #  }
        }
      }
    }
    #cleanup
    
  }
  #final cleanup
  #save images
  dev.off()
  if (length(alteration_list)>0)
  {
    colnames(dm_per_cell_vals)[2:ncol(dm_per_cell_vals)]<-alteration_list
    
    #cutoff for minimum # of cells
    #count # of cells in each cluster
    #rowSums(.[2:ncol(cellQuality)]>0)
    loss_cluster_counts<-dm_per_cell_vals %>% gather(Alteration,Clust,2:ncol(dm_per_cell_vals)) %>% group_by(Alteration) %>% count(Alteration,Clust)
    #  print(loss_cluster_counts)
    min_loss<-loss_cluster_counts %>% spread(Clust,n) %>% select(-'0') %>% mutate(min=min(`1`,`2`))
    #  print(min_loss)
    min_loss[is.na(min_loss)]<-0
    #min_loss$Alteration[min_loss$min>lossCutoffCells]
    #cut appropriately
    dm_per_cell_vals <- dm_per_cell_vals %>% select(cellName,min_loss$Alteration[min_loss$min>lossCutoffCells]) %>% mutate_at(vars(starts_with('chr')),funs(if_else(.==1,-1,.))) %>% mutate_at(vars(starts_with('chr')),funs(if_else(.==2,1,.)))
    return(list(dm_per_cell_vals,alteration_list,alteration_delta))
  }
  return(list())
}
#' AnnotateLosses
#'
#' This function allows you to annotate loss regions
#' @param lossFile the file for LOH putative reginos
#' @keywords Losses
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
  scData<-scData %>% dplyr::select(-Cell_id)
  return(scData)
}

#' smoothClusters
#' Read clusters from file and CNV calls from input or an external file
#' Note: assumes files are comma delimited
#' Returns a smoothed matrix
#' @param inputClusterFile   input cluster matrix
#' @param inputCNVList   list of CNVs with copy number estimates (2 as normal)
#' @param inputCNVClusterFile  alternate option to specify external comma delimited file for CNV calls
#' @param percentPositive   minimum percentage of a cluster to use to call positive for CNV, default is 0.5 (50 percent)
#' @keywords CNV
#' @keywords clusters
#' @export
smoothClusters <- function(inputClusterFile,inputCNVList,inputCNVClusterFile="",percentPositive=0.5,removeEmpty=TRUE) 
{
  inputClusters<-read.table(inputClusterFile,stringsAsFactors=FALSE,header=TRUE,sep=",")
  colnames(inputClusters)[2]<-"clust"
  inputCNV<-""
  if (inputCNVClusterFile!="")
  {
    inputCNV<-read.table(inputCNVClusterFile,stringsAsFactors=FALSE,header=TRUE,sep=",")
    inputCNV <- inputCNV %>% mutate_at(vars(starts_with("chr")),list(as.double))
    #print(inputCNV
    
  }
  else
  {
    #input from CNV calls
    inputCNV<-inputCNVList[[1]]
    #colnames(inputCNV)[1]<-"Barcode"
  }
  colnames(inputCNV)[1]<-"Barcode"
  #print(colnames(inputClusters))
  #print(inputClusters$Barcode)
  
  #colnames(inputClusters)[1]<-"Barcode"
  #print(inputClusters[,1])
  b1<-left_join(inputClusters,inputCNV,by="Barcode")
  #SMOOTH
  b1c<-b1 %>% mutate(clust=inputClusters[,2])
  specRound <- function(x,boost=0.1)
  {
    #  print(x-2)
    # print(sign(x-2)*round(abs(x-2)+boost))
    #print(abs(x-2)+boost)
    return (2+sign(x-2)*round(abs(x-2)+boost,digits=0))
  }
  #b1c %>% mutate_at(vars(starts_with("chr")),funs(if_else(is.na(.),2L,.)))
  #b1c_round<-b1c %>% mutate_at(vars(starts_with("chr")),funs(if_else(is.na(.),2L,.))) %>% group_by(clust) %>% summarise_at(vars(starts_with("chr")),list(mean))  %>% mutate_at(vars(starts_with("chr")),funs(specRound(.,boost=-0.1)))
  b1c_round<-b1c %>% mutate_at(vars(starts_with("chr")),funs(if_else(is.na(.),2,.))) %>% 
    group_by(clust) %>% summarise_at(vars(starts_with("chr")),list(mean))  %>% 
    mutate_at(vars(starts_with("chr")),funs(specRound(.,0.5 - percentPositive)))
  if (removeEmpty)
  {
  isMultiple<-function(x){
    x2<-as.numeric(x)
    if (max(x2)==min(x2)){
      return(TRUE)
    }
    else
    {
      return(FALSE)
    }}
  bb<-data.table(b1c_round)
 # print(bb)
  #print(colnames(bb)[which(bb[,lapply(.SD,isMultiple)]==FALSE),.SDcols=colnames(bb)[which(str_detect(colnames(bb),"chr"))])
  bb[,.SD,.SDcols=colnames(bb)[which(bb[,lapply(.SD,isMultiple),.SDcols=-c("clust")]==FALSE)]]
  b1c_round<-bb
  }
  cluster_clean<-left_join(inputClusters,b1c_round,by="clust")
  return(cluster_clean)
}
#' generateReferences
#' Generate reference files from UCSC for genome of interest
#' Returns CPG, chromosome size and cytoband data into target directory
#' @param genomeObject   BSgenome object of interest
#' @param genomeText   Shorthand for genome in UCSC (defaults to hg38)
#' @param tileWidth  Width of tiles - (default: 1e6)
#' @param outputDir   Output directory to save reference files to (default: ~)
#' @keywords CNV
#' @keywords reference
#' @export
generateReferences <- function(genomeObject,genomeText="hg38",tileWidth=1e6,outputDir="~")
{
  fileSuffixes=str_c(outputDir,"/",genomeText,"_",tileWidth,sep="")
  print(str_c("Output to ",fileSuffixes,sep=""))
  #chrom sizes - 
  chrom_sizes<-rownames_to_column(data.frame(length=seqlengths(genomeObject)),var="chrom") %>% mutate(keeper=str_detect(chrom,pattern="alt|random|chrUn|fix|_|MT",negate=TRUE)) %>% filter(keeper==TRUE) %>% select(-keeper)
  print(chrom_sizes)
  chroms<-GRanges(seqnames=genomeObject@seqinfo)
  #print(which(chroms@seqnames %in% chrom_sizes$chrom))
  #chroms<-chroms[which(chroms@seqnames %in% chrom_sizes$chrom)]
  #chroms<-chroms[which(chroms@seqnames %in% chrom_sizes$chrom)]
  
  #print(chroms)
  #https://bioconductor.org/packages/devel/bioc/vignettes/rtracklayer/inst/doc/rtracklayer.pdf
  #print(seqlengths(genomeObject))
  #a<-tile(chroms,width=1e6)
  #tile the genome
  tiles<-tileGenome(seqlengths(chroms),tilewidth=tileWidth,cut.last.tile.in.chrom=T)
  print(tiles)
  mySession = browserSession("UCSC")
  genome(mySession) <- genomeText
  tbl.cytobands <- getTable(
    ucscTableQuery(mySession, track="cytoBand",
                   table="cytoBand"))
  tbl.cpgIslandExt <- getTable(
    ucscTableQuery(mySession, track="cpgIslandExtUnmasked",
                   table="cpgIslandExtUnmasked"))
  
  tbl.cytobands
  # print(tbl.cytobands)
  tbl.cytobands<-GRanges(tbl.cytobands)
  tbl.cpgIslandExt<-GRanges(tbl.cpgIslandExt)
  #now intersect these to generate the reference files
  #http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/findOverlaps-methods.html
  #MATCHES FOR Cytoband
  centromeres<-which(mcols(tbl.cytobands)$gieStain=="acen")
  mcols(tbl.cytobands)$name[centromeres]<-"cen"
  mcols(tbl.cytobands)$gieStain<-NULL
  print(tbl.cytobands$name)
  tbl.cytobands$name<-as.character(tbl.cytobands$name)
  #pad names
  #print(which(tbl.cytobands$name==""))
  tbl.cytobands$name[which(tbl.cytobands$name=="")]="-"
  tbl.cytobands$name<-str_remove(tbl.cytobands$name,"[0-9.]+")
  matches<-findOverlaps(tbl.cytobands,tiles,select="all",ignore.strand=TRUE)
  #matches
  a3<-tiles[subjectHits(matches)]
 # subjectHits(matches)
  # print(tiles)
  # print(a3)
  
  #print(tbl.cytobands[queryHits(matches)])
  #tbl.cytobands
  mcols(a3) <- cbind.data.frame(
    mcols(a3),
    mcols(tbl.cytobands[queryHits(matches)]))
  #head(a3,n=50)
  #label Cytobands
  #https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/GRanges-class.html
  #ALL CLEAN :)
 #mcols(empties)<-cbind.data.frame(mcols(empties),name="")
  #empties
  empties<-tiles[-unique(queryHits(matches))]
  mcols(empties)<-cbind.data.frame(mcols(empties),name="p")
  a3<-sort.GenomicRanges(append(GRanges(a3),empties))
  a3<-unique(a3)
  
  #add empties to zero (cpgNum = 0)
  
  #a3<-sort.GenomicRanges(append(GRanges(a3),empties))
  
  #overlap stuff again to remove doubles
  
  #label CPGs
  #head(tbl.cpgIslandExt)
  tbl.cpgIslandExt<-sort.GenomicRanges(tbl.cpgIslandExt)
  matches<-findOverlaps(tiles,tbl.cpgIslandExt,select="all",ignore.strand=TRUE)
  #sum up the stuff
  #https://www.rdocumentation.org/packages/S4Vectors/versions/0.4.0/topics/Hits-class
  
  b1<-tbl.cpgIslandExt[subjectHits(matches)]
  mcols(b1) <- cbind.data.frame(
    mcols(b1),
    ranges(tiles[queryHits(matches)]),
    chrom(tiles[queryHits(matches)]))
  #print(b1)
  #print(tiles[queryHits(matches)])
  
  
  #NEED TO FIND TILES THAT AREN'T MATCHES
  length(unique(queryHits(matches)))
  #
  b1t<-as_tibble(b1)
  #subjectHits vs queryHits
  b1t<-b1t %>% mutate(interval=str_c(chrom.tiles.queryHits.matches...,start.1,end.1,sep="-"))
  b1t$interval
  b1t2<-b1t %>% group_by(interval) %>% summarise_at(vars(cpgNum),sum) %>% separate(interval,into=c("chrom","start","end"),sep="-")
  #missing a few here
  #GET EMPTIES
  empties<-tiles[-unique(queryHits(matches))]
  mcols(empties)<-cbind.data.frame(mcols(empties),cpgNum=0)
  #add empties to zero (cpgNum = 0)
  
  cpg_densities<-sort.GenomicRanges(append(GRanges(b1t2),empties))
  #TODO: fill with missing zeros
  
  #REMOVE GARBAGE CHROMOSOMES
  cpg_densities<-keepStandardChromosomes(cpg_densities,pruning.mode = "tidy")
  cytoband_densities<-keepStandardChromosomes(a3,pruning.mode="tidy")
  sort.GenomicRanges(cytoband_densities)
  
  write.table(as_tibble(sort.GenomicRanges(cpg_densities))[,c(1:3,6)],str_c(fileSuffixes,"_cpg_densities.tsv",sep=""),quote=FALSE,sep="\t",col.names = FALSE,row.names=FALSE)
  write.table(as_tibble(sort.GenomicRanges(cytoband_densities))[,c(1:3,6)],str_c(fileSuffixes,"_cytoband_densities_granges.tsv",sep=""),quote=FALSE,sep="\t",col.names = FALSE,row.names=FALSE)
  write.table(as_tibble(chrom_sizes),str_c(outputDir,"/",genomeText,"_chrom_sizes.tsv",sep=""),quote=FALSE,sep="\t",col.names = FALSE,row.names=FALSE)
}
#' identifyNonNeoplastic
#' Detect neoplastic vs neoplastic cells using unsupervised clustering
#' Returns list of barcodes and cluster assignments, and list of non-neoplastic cells
#' @param inputMatrix   Matrix (after filterCells is called)
#' @param estimatedCellularity   Expected cellularity of tumour (default 0.8)
#' @param nmfComponents   # of components to use for NMF/SVD
#' @param outputHeatmap   Output a heatmap of the clustering (default is yes)
#' @param cutHeight  Cutting height as percent of max for dendrogram
#' @param methodHclust  Method for clustering (default is Ward.D)
#' @keywords CNV
#' @keywords clustering
#' @export
identifyNonNeoplastic <- function(inputMatrix,estimatedCellularity=0.8,nmfComponents=5,outputHeatmap=TRUE,cutHeight=0.6,methodHclust="ward.D")
{
  #uses NMF and fast hclust packages
  message("Running NMF")
  res <- nmf(column_to_rownames(inputMatrix,var="chrom"), c(nmfComponents), "brunet", seed="nndsvd",.stop=nmf.stop.threshold(0.1),maxIter=2500)
  message("Computing clusters")
  #dist<-dist(t(coef(res)),method="euclidean")
  tscale<-scale(x=t(coef(res)),center=TRUE)
  dist=dist(tscale,method="euclidean")
  
  clusts<-fastcluster::hclust(dist,method = methodHclust)
  #average or ward
  #clusts<-agnes(t(coef(res)),diss = FALSE,metric="euclidean",method="ward")
  #YEAH THIS WORKS AWESOME
  #some tricksy samples may need 4
  
  #add PDF
  cell_assigns<-cutree(clusts,h=max(clusts$height)*cutHeight)
  
  if (outputHeatmap==TRUE)
  {
    pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_heatmap.pdf"),width=6,height=6)
    #print(cell_assigns)
    #print(factor(cell_assigns))
    heatmap.2(coef(res),Rowv=FALSE,Colv=as.dendrogram(clusts),dendrogram="column",density.info="none",trace="none",scale="none",labCol=FALSE,col=colorRampPalette(viridis(5)),symkey=FALSE,useRaster=TRUE,ColSideColors=viridis(n=length(unique(as.character(cell_assigns))))[cell_assigns])
    legend("topright",fill=viridis(n=3),x.intersp = 0.8,y.intersp=0.8,legend=unique(as.character(cell_assigns)),horiz = FALSE,cex = 0.9,border=TRUE,bty="n")
    dev.off()
  }
  cluster_order<-data.frame(cluster=cell_assigns) %>% group_by(cluster) %>% summarise(number=n()) %>% arrange(number)
  tumor_cell_ids=names(which(cell_assigns!=cluster_order$cluster[1]))
  normalCluster=cluster_order$cluster[1]
  normal_cell_ids=names(which(cell_assigns==cluster_order$cluster[1]))
  #flip if less than 50% tumor
  if (estimatedCellularity<0.5)
  {
    tumor_cell_ids=names(which(cell_assigns!=cluster_order$cluster[nrow(cluster_order)]))
    normal_cell_ids=names(which(cell_assigns==cluster_order$cluster[nrow(cluster_order)]))
    normalCluster=cluster_order$cluster[nrow(cluster_order)]
  }
  return(list(cellAssigns=cell_assigns,normalBarcodes=normal_cell_ids,clusterNormal=normalCluster))
}