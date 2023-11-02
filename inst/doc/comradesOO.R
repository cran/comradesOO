## ----include = FALSE----------------------------------------------------------
library(rmarkdown)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)




## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo = T, results = 'hide',message=FALSE---------------------------------
# Load the comrades-OO Library 
library(comradesOO)

## ----echo = F, results = 'hide',message=FALSE---------------------------------
# Here are the other libraries on which comradesOO relies
library(seqinr)
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(MASS)
library(ggplot2)
library(doParallel)
library(igraph)
library(R4RNA)
library(RColorBrewer)
library(heatmap3)
library(mixtools)
library(TopDom)
library(tidyverse)
library(RRNA)
library(ggrepel)

## ----echo = T, results = 'hide',message=FALSE---------------------------------
# Load the comrades-OO Library 
library(comradesOO)

## -----------------------------------------------------------------------------
# Set up the sample table
sampleTableRow1 = c(system.file("extdata", 
                                "s1.txt", 
                                package="comradesOO"), "s", "1", "s1")
sampleTableRow2 = c(system.file("extdata", 
                                "c1.txt", 
                                package="comradesOO"), "c", "1", "c1")
sampleTable2 = rbind.data.frame(sampleTableRow1, sampleTableRow2)

# add the column names 
colnames(sampleTable2) = c("file", "group", "sample", "sampleName")

sampleTable2

## -----------------------------------------------------------------------------
rna = c("ENSG000000XXXXX_NR003286-2_RN18S1_rRNA")

## -----------------------------------------------------------------------------
path18SFata <- system.file("extdata", 
                           "18S.fasta", 
                           package="comradesOO")

rnaRefs = list()
rnaRefs[[rna]] = read.fasta(path18SFata)


## -----------------------------------------------------------------------------
path18SFata <- system.file("extdata", 
                           "ribovision18S.txt", 
                           package="comradesOO")
known18S = read.table(path18SFata,
                      header = F)

## -----------------------------------------------------------------------------
pathShape <- system.file("extdata",
                         "reactivities.txt", 
                         package="comradesOO")
shape = read.table(pathShape,
                      header = F)

## -----------------------------------------------------------------------------
# runComradesOO(rna,
#                     rnaSize =0 ,
#                     sampleTable,
#                     cores = 3,
#                     stepCount = 2,
#                     clusterCutoff = 20,
#                     clusteredCds,
#                     trimFactor = 2.5, 
#                     clusterCutoff = 1,
#                     rnaRefs,
#                     start,
#                     end,
#                     evCutoff = 1,
#                     ensembl = 50,
#                     constraintNumber = 20,
#                     shape = 0)


## -----------------------------------------------------------------------------
# load the object
cds = comradesDataSet(rnas = rna,
                      sampleTable = sampleTable2)

## -----------------------------------------------------------------------------
# Check status of instance 
cds

## -----------------------------------------------------------------------------

# Returns the size of the RNA
rnaSize(cds)

## -----------------------------------------------------------------------------
# Returns the sample table 
sampleTable(cds)

## -----------------------------------------------------------------------------
# Returns indexes of the samples in the control and not control groups
group(cds)

## -----------------------------------------------------------------------------
# Get the sample names of the instance
sampleNames(cds)

## ----eval=FALSE,echo=FALSE----------------------------------------------------
#  # Return the hybFiles slot
#  hybFiles(cds)

## ----eval=FALSE, echo=FALSE---------------------------------------------------
#  # Return the matrixList slotmatrixList(cds)

## -----------------------------------------------------------------------------
data = getData(cds,       # The object      
        "hybFiles",       # The Type of data to return     
        "original")[[1]]  # The stage of the analysis for the return data
head(data)

## -----------------------------------------------------------------------------

# Returns the RNAs with highest number of assigned reads 
topTranscripts(cds, # The comradesDataSet instance
               2)  # The number of entried to return


## -----------------------------------------------------------------------------
# Returns the RNAs that interact with the RNA of interest
topInteracters(cds, # The comradesDataSet instance
               1)   # The number of entries to return

## -----------------------------------------------------------------------------
# Returns the Interacions with the highest number of assigned reads
topInteractions(cds, # The comradesDataSet instance
                2)   # The number of entries to return


## -----------------------------------------------------------------------------
# Cluster the reads
clusteredCds = clusterComrades(cds = cds,         # The comradesDataSet instance 
                               cores = 1,         # The number of cores
                               stepCount = 2,     # The number of steps in the random walk
                               clusterCutoff = 3) # The minimum number of reads for a cluster to be considered

## -----------------------------------------------------------------------------
# Check status of instance 
clusteredCds

## -----------------------------------------------------------------------------
# Returns the number of clusters in each sample
#clusterNumbers(clusteredCds)

## -----------------------------------------------------------------------------
# Returns the number reads in clusters
readNumbers( clusteredCds)

## -----------------------------------------------------------------------------
getData(clusteredCds,        # The object             
        "clusterTableList",  # The Type of data to return     
        "original")[[1]]     # The stage of the analysis for the return data

## ----eval=FALSE, echo=FALSE---------------------------------------------------
#  getData(clusteredCds,         # The object
#          rna,                  # The rna of interest
#          "clusterGrangesList", # The Type of data to return
#          "original")[[1]]      # The stage of the analysis for the return data

## -----------------------------------------------------------------------------
# Trim the Clusters
trimmedClusters = trimClusters(clusteredCds = clusteredCds, # The comradesDataSet instance 
                               trimFactor = 1,              # The cutoff for cluster trimming (see above)
                               clusterCutoff = 30)          # The minimum number of reads for a cluster to be considered

## -----------------------------------------------------------------------------
# Check status of instance 
trimmedClusters

## ----error=FALSE,eval = FALSE, results=FALSE,message=FALSE--------------------
#  # Fold the RNA in part of whole
#  foldedCds = foldComrades(trimmedClusters,
#                           rna = rna,
#                           rnaRefs = rnaRefs,
#                           start = 1700,
#                           end = 1869,
#                           shape = 0,
#                           ensembl = 40,
#                           constraintNumber  = 5,
#                           evCutoff = 50)

## ----eval = FALSE-------------------------------------------------------------
#  # Check status of instance
#  foldedCds

## -----------------------------------------------------------------------------
# Plot heatmaps for each sample
plotMatrices(cds = cds,         # The comradesDataSet instance 
             type = "original", # The "analysis stage"
             directory = 0,     # The directory for output (0 for standard out)
             a = 1,             # Start coord for x-axis
             b = rnaSize(cds),  # End coord for x-axis
             c = 1,             # Start coord for y-axis
             d = rnaSize(cds),  # End coord for y-axis
             h = 5)             # The hight of the image (if saved)


## -----------------------------------------------------------------------------
# Plot heatmaps for all samples combined and all controls combined
plotMatricesAverage(cds = cds, # The comradesDataSet instance 
             type = "original", # The "analysis stage"
             directory = 0,     # The directory for output (0 for standard out)
             a = 1,             # Start coord for x-axis
             b = rnaSize(cds),  # End coord for x-axis
             c = 1,             # Start coord for y-axis
             d = rnaSize(cds),  # End coord for y-axis
             h = 5)             # The hight of the image (if saved)


## ----eval = FALSE-------------------------------------------------------------
#  domainDF = data.frame()
#  for(j in c(20,30,40,50,60,70)){
#      #for(i in which(sampleTable(cds)$group == "s")){
#  
#      timeMats = as.matrix(getData(x = cds,
#                                   data = "matrixList",
#                                   type = "noHost")[[1]])
#  
#      timeMats = timeMats/ (sum(timeMats)/1000000)
#      tmp = tempfile()
#      write.table(timeMats, file = tmp,quote = F,row.names = F, col.names = F)
#  
#      tdData2 = readHiC(
#          file = tmp,
#          chr = "rna18s",
#          binSize = 10,
#          debug = getOption("TopDom.debug", FALSE)
#      )
#  
#      tdData =  TopDom(
#          tdData2 ,
#          window.size = j,
#          outFile = NULL,
#          statFilter = TRUE,
#          debug = getOption("TopDom.debug", FALSE)
#      )
#  
#      td = tdData$domain
#      td$sample = sampleTable(cds)$sampleName[1]
#      td$window = j
#      domainDF = rbind.data.frame(td, domainDF)
#  
#  }
#  
#  
#  
#  ggplot(domainDF) +
#      geom_segment(aes(x = from.coord/10,
#                       xend = to.coord/10, y = as.factor(sub("s","",sample)),
#                       yend = (as.factor(sub("s","",sample)) ), colour = tag),
#                   size  = 20, alpha = 0.8) +
#      facet_grid(window~.)+
#      theme_bw()
#  

## ----eval = FALSE-------------------------------------------------------------
#  plotEnsemblePCA(foldedCds,
#                  labels = T, # plot labels for structures
#                  split = F)  # split samples over different facets (T/f)
#  

## ----eval = FALSE-------------------------------------------------------------
#  plotComparisonArc(foldedCds = foldedCds,
#                    s1 = "s1",            # The sample of the 1st structure
#                    s2 = "s1",            # The sample of the 2nd structure
#                    n1 = 1,               # The number of the 1st structure
#                    n2 = 2)               # The number of the 2nd structure

## ----eval = F-----------------------------------------------------------------
#  plotStructure(foldedCds = foldedCds,
#                rnaRefs = rnaRefs,
#                s = "s1",          # The sample of the structure
#                n = 1)             # The number of the structure

## ----message=FALSE, echo=F,eval = FALSE, results='hide'-----------------------
#  
#  plotStructure(foldedCds = foldedCds,
#                rnaRefs = rnaRefs,
#                s =  "s1",        # The sample of the structure
#                n = 1)            # The number of the structure
#  

## -----------------------------------------------------------------------------

getInteractions(cds,
                "ENSG00000XXXXXX_NR003287-2_RN28S1_rRNA") %>%
    mutate(sample =sub("\\d$","",sample) )%>%
    group_by(rna,Position,sample)%>%
    summarise(sum =  sum(depth)) %>%
    ggplot()+
    geom_area(aes(x = Position,
                  y = sum, 
                  fill = sample), 
              stat = "identity")+
    facet_grid(sample~.) +
    theme_bw()


## -----------------------------------------------------------------------------
getReverseInteractions(cds,
                       rna) %>%
    mutate(sample =sub("\\d$","",sample) )%>%
    group_by(rna,Position,sample)%>%
    summarise(sum =  sum(depth)) %>%
    ggplot()+
    geom_area(aes(x = Position,
                  y = sum, 
                  fill = sample), 
                    stat = "identity")+
    facet_grid(sample~.)+
    theme_bw()

## -----------------------------------------------------------------------------
expansionSize = 5
knownMat = matrix(0, nrow = rnaSize(cds), ncol = rnaSize(cds))
for(i in 1:nrow(known18S)){
    knownMat[ (known18S$V1[i]-expansionSize):(known18S$V1[i]+expansionSize),
              (known18S$V2[i]-expansionSize):(known18S$V2[i]+expansionSize)] =
        knownMat[(known18S$V1[i]-expansionSize):(known18S$V1[i]+expansionSize),
                 (known18S$V2[i]-expansionSize):(known18S$V2[i]+expansionSize)] +1
}
knownMat = knownMat + t(knownMat)


## -----------------------------------------------------------------------------
# use compare known to gett he known and not know clusters
knowClusteredCds = compareKnown(trimmedClusters, # The comradesDataSet instance 
                                knownMat, # A contact matrix of know interactions
                                "trimmedClusters") # The analysis stage of clustering to compare 

knowClusteredCds

## -----------------------------------------------------------------------------
# Plot heatmaps for all samples combined and all controls combined
plotMatricesAverage(cds = knowClusteredCds, # The comradesDataSet instance 
             type = "KnownAndNovel", # The "analysis stage"
             directory = 0,     # The directory for output (0 for standard out)
             a = 1,             # Start coord for x-axis
             b = rnaSize(cds),  # End coord for x-axis
             c = 1,             # Start coord for y-axis
             d = rnaSize(cds),  # End coord for y-axis
             h = 5)             # The hight of the image (if saved)


## -----------------------------------------------------------------------------
# Get the number of clusters for each analysis Stage
clusterNumbers(knowClusteredCds)

## -----------------------------------------------------------------------------
# Get the number of reads in each cluster for each analysis stage
readNumbers(knowClusteredCds)

## ----eval = FALSE-------------------------------------------------------------
#  head(compareKnownStructures(foldedCds,
#                              known18S)) # the comarison set

## ----eval = FALSE-------------------------------------------------------------
#  ggplot(compareKnownStructures(foldedCds, known18S)) +
#      geom_hline(yintercept = c(0.5,0.25,0.75,0,1),
#                 colour = "grey",
#                 alpha = 0.2)+
#      geom_vline(xintercept = c(0.5,0.25,0.75,0,1),
#                 colour = "grey",
#                 alpha = 0.2)+
#      geom_point(aes(x = sensitivity,
#                     y = precision,
#                     size = as.numeric(as.character(unlist(foldedCds@dgs))),
#                     colour = str_sub(structureID,
#                                      start = 1 ,
#                                      end = 2))) +
#      xlim(0,1)+
#      ylim(0,1)+
#      theme_classic()
#  

