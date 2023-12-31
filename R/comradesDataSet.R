#' @include comradesOO.R


#' comradesDataSet
#'
#' An S4 class to represent a COMRADES dataset
#' @rdname comradesDataSet
#' @export
#'
setClass(
  "comradesDataSet",
  slots = c(
    rnas = "character",
    rnaSize = "numeric",
    sampleTable = "data.frame",
    # meta data
    hybFiles = "list",
    # data tables
    matrixList = "list",
    clusterGrangesList = "list",
    clusterTableList = "list",
    clusterTableFolded = "data.frame",
    interactionTable = "data.frame",
    viennaStructures = "list",
    dgs = "list"
  ),
  prototype = list()
)

setValidity("comradesDataSet", function(object) {
  
})




#' comradesDataSet
#'
#' \code{comradesDataSet} objects are used to store the input meta-data, data and
#' create a framework for the storage of results. Whilst creating the object,
#' the original hyb files are also filtered for the RNA of interest. Check the 
#' package vignette for more information.
#'
#'
#' @param rnas vector - The names of the RNA interest, these must be displayed
#' the same way as in the input Hyb Files.
#' @param rnaSize named list - The sizes (nt) of the RNAs of interest, the list
#'  elements must have same names as the \code{rnas} vector and each each contain
#'  one numeric value.
#' @param sampleTable string - The address of the sample table, the sample table
#'  must have 4 columns, fileName (the full path and file name of the input
#'  hyb file for each sample ), group ("s" - sample or "c" - control),
#'  sample (1,2,3, etc), sampleName (must be unique).
#'
#' @return A comradesDataSet object.
#'
#'
#'
#' @slot clusterTableFolded table - a table similar to the \code{clusterTableList}
#' it contains coordinates of the clusters along with vienna format fold and
#' RNA sequences for each cluster
#' @slot clusterTableList List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. contains a table
#' with coordinates and information about the clusters identified
#' @slot clusterGrangesList List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. contains GRanges
#' objects of the original duplexes with their cluster membership
#' @slot sampleTable table - Column names; fileName, group (s or c),
#'  sample (1,2,3, etc), sampleName (must be unique)
#' @slot rnas string - a single RNA to analyse - must be present in \code{rnas(cdsObject)}
#' @slot rnaSize if set to 0 this will be calculated 
#' @slot matrixList List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{matrixList(cds)[[rna]][[type]][[sample]]}. Contains a set
#' of contact matrices, each cell contains the number of duplexes identified
#' for position x,y.
#' @slot hybFiles List - Follows the pattern for list slots of comradesDataSet
#' objects, \code{hybFiles(cds)[[rna]][[type]][[sample]]}. Contains a set of
#' tables, these are the original Hyb files that were read in.
#' @slot interactionTable Table of interactions discovered in step1 of the folding
#' @slot viennaStructures List of vienna format structures from final prediction
#' @slot dgs List of free energies
#'
#' @name comradesDataSet
#' @docType class
#' @rdname comradesDataSet
#'
#' @examples 
#' # make example input
#' c4 = c(rep("transcript1",100),rep("transcript2",100) )
#' c10 = c(rep("transcript1",200) )
#' c1 = 1:200
#' c2 = rep(paste(rep("A", 40), collapse = ""),200)
#' c3 = rep(".",200)
#' c9 = rep(".",200)
#' c15 = rep(".",200)
#' c5 = rep(1,200)
#' c11 = rep(21,200)
#' c6 = rep(20,200)
#' c12= rep(40,200)
#' # short distance 50
#' c7 = sample(1:5, 50, replace = TRUE)
#' c8 = sample(20:25, 50, replace = TRUE)
#' c13 = sample(20:25, 50, replace = TRUE)
#' c14 = sample(40:45, 50, replace = TRUE)
#' # long distance 50
#' c7 = c(c7,sample(1:5, 50, replace = TRUE))
#' c8 = c(c8,sample(20:25, 50, replace = TRUE))
#' c13 = c(c13,sample(60:70, 50, replace = TRUE))
#' c14 = c(c14,sample(80:83, 50, replace = TRUE))
#' # inter RNA 100
#' c7 = c(c7,sample(1:5, 100, replace = TRUE))
#' c8 = c(c8,sample(20:25, 100, replace = TRUE))
#' c13 = c(c13,sample(1:5, 100, replace = TRUE))
#' c14 = c(c14,sample(20:25, 100, replace = TRUE))
#' exampleInput = data.frame(V1 = c1,
#'                          V2 = c2,
#'                          V3 = c3,
#'                          V4 = c4,
#'                          V5 = as.numeric(c5),
#'                          V6 = as.numeric(c6),
#'                          V7 = as.numeric(c7),
#'                          V8 = as.numeric(c8),
#'                          V9 = c9,
#'                          V10 = c10,
#'                          V11 = as.numeric(c11),
#'                          V12 = as.numeric(c12),
#'                          V13 = as.numeric(c13),
#'                          V14 = as.numeric(c14),
#'                          V15 = c15)
#' file = tempfile()
#' write.table(exampleInput,file = file, 
#' quote = FALSE, 
#' row.names = FALSE, sep = "\t", col.names = FALSE)
#'
#'
#' # Set up the sample table 
#' sampleTabler1 = c(file, "s", "1", "s1")
#' sampleTabler2 = c(file, "c", "1", "c1")
#' # make the sample table 
#' sampleTable2 = rbind.data.frame(sampleTabler1, sampleTabler2)
#' # add the column names 
#' colnames(sampleTable2) = c("file", "group", "sample", "sampleName")
#'
#' # Choose RNA and set up the object ----
#' rna = c("transcript1")
#' 
#' 
#' # load the object
#' cds = comradesDataSet(rnas = rna,
#'                      rnaSize = 0,
#'                      sampleTable = sampleTable2)
#'
#'cds
#'
#' @export
comradesDataSet <- function(rnas,
                            rnaSize = 0 ,
                            sampleTable) {
  ###########################################################
  # Read in the sample table
  ###########################################################
  # check the inputs here, stop if wrong
  message(" ******************************************** ")
  message(" *****            COMRADES-OO          ****** ")
  message(" ******************************************** ")
  message(" *****-------*******************-------****** ")
  message(" *****       Reading SampleTable       ****** ")
  
  
  # Read in sample table

  #check for more than two samples
  # if( nrow(sampleTable) < 2 ){
  #        stop( "The sample Table must contain at least 1
  #              sample and 1 control" )
  #    }
  message(paste(" *****       Detected ",
                nrow(sampleTable), " Samples      ****** "))
  
  
  
  ###########################################################
  #check column names of sampleTable
  colnamesST = c("file", "group", "sample", "sampleName")
  if (all(colnames(sampleTable) != colnamesST)) {
    stop("Column names of metaData table should be :
              file, group, sample, sampleNames")
  }
  
  # Get the comparison groups check group has the c and s
  if (!(unique(as.character(sampleTable$group))[1] %in% c("c", "s") &
        unique(as.character(sampleTable$group))[2] %in% c("c", "s"))) {
    stop("Groups should be c and s")
  }
  
  
  # Make group into a list with control and sample
  group = sampleTable[, "group"]
  group2 = list()
  group2[["c"]] = which(group == "c")
  group2[["s"]] = which(group == "s")
  
  group = group2
  message(paste(
    " *****     detected group c::",
    paste(group[["c"]],
          collapse = " ") ,
    paste(rep(" ",
              (
                length(group[["c"]]) * (3 - length(group[["c"]]))) * 2),
          collapse = ""),
    "   ***** "
  ))
  message(paste(
    " *****     detected group s::",
    paste(group[["s"]],
          collapse = " ") ,
    paste(rep(" ",
              (
                length(group[["s"]]) * (3 - (length(group[["s"]]))) * 2)),
          collapse = ""),
    "   ***** "
  ))
  
  
  
  ###########################################################
  # Get the sampleNames
  sampleNames = c()
  if (is.null(sampleTable$sampleName)) {
    stop("The sample Table must have a column named sampleName")
  } else if (length(unique(sampleTable$sampleName)) !=
             length(sampleTable$sampleName)) {
    stop("Sample names must be unique")
  } else{
    sampleNames = as.character(sampleTable$sampleName)
    spaces =  (length(sampleNames) * (3 - length(sampleNames))) * 2 
    if(spaces < 0 ){spaces = 0}
    message(paste(
      " ****  ",
      paste(rep(" ",
                spaces,
            collapse = ""),
      " Sample Names: ",
      paste(sampleNames, collapse = " "),
      paste(rep(" ",
                spaces,
            collapse = ""),
      " **** "
    ))))
  }
  
  
  
  ###########################################################
  # Read in the  hyb files
  ###########################################################
  #load the files into a list
  message(" *****         Reading Hyb Files        ***** ")
  
  hybFiles = list()
  hybFiles[["all"]] = list()
  hybFiles[["all"]][["all"]] = list()
  
  
  
  # read in the tables
  inputs <- lapply(as.character(sampleTable$file),
                   function(file)
                     read.table(file,colClasses = c("character",
                                                    "character",
                                                    "character",
                                                    "character",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "character",
                                                    "character",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "integer",
                                                    "character")))
  
  #check column names
  if (all(sapply(inputs, function(file)
    ! (identical(
      colnames(file),
      c(
        "V1",
        "V2",
        "V3",
        "V4",
        "V5",
        "V6",
        "V7",
        "V8",
        "V9",
        "V10",
        "V11",
        "V12",
        "V13",
        "V14",
        "V15"
      )
    ))))) {
    stop(
      " The input files do not look they are produced with the
                 nextflow pipeline, please check the documentation. "
    )
  }
  hybFiles[["all"]][["all"]] = inputs
  names(hybFiles[["all"]][["all"]]) = sampleNames
  
  
  #get the maxpos
  getMax = function(file, rna) {
    return(max(file[file$V4 == rnas & file$V10 == rnas, "V14"],
               max(file[file$V4 == rnas & file$V10 == rnas, "V8"])))
  }
  
  maxPos = max(sapply(inputs, getMax, rna = rnas))
  
  
  
  
  
  
  ###########################################################
  message(" *****     Getting RNAs of Interest    ****** ")
  message(" *****    RNA of interest + Host RNA    ***** ")
  hybFiles[[rnas]][["original"]] = swapHybs(hybList = hybFiles[["all"]][["all"]],
                                            rna = rnas)
  hybFiles[[rnas]][["host"]] = swapHybs3(hybList = hybFiles[["all"]][["all"]],
                                         rna = rnas)
  message(" *****      RNA of interest Alone       ***** ")
  hybFiles[[rnas]][["noHost"]] = swapHybs2(hybList = hybFiles[["all"]][["all"]],
                                           rna = rnas)
  
  
  
  ###########################################################
  # Make matrices of the specific RNA without host
  ###########################################################
  message(" *****         Making Matrices         ****** ")
  
  matrixList = list()
  matrixList[[rnas]] = list()
  if (rnaSize > 0) {
    rnaSize2 =   rnaSize
  } else{
    rnaSize2 = maxPos
  }
  
  message(paste(" *****          RNA Size: ",
                rnaSize2, "        ***** "))
  matrixList[[rnas]][["noHost"]] = list()
  matrixList[[rnas]][["noHost"]] = getMatrices(hybFiles[[rnas]][["noHost"]],
                                               rnas, rnaSize2)
  names(matrixList[[rnas]][["noHost"]]) = sampleNames
  matrixList[[rnas]][["original"]] = list()
  matrixList[[rnas]][["original"]] = getMatrices(hybFiles[[rnas]][["original"]],
                                                 rnas, rnaSize2)
  names(matrixList[[rnas]][["original"]]) = sampleNames
  
  
  
  
  ###########################################################
  message(" *****         Creating object          ***** ")
  message(" *****-------*******************-------****** ")
  message(" ******************************************** ")
  message(" ******************************************** ")
  #create comrades dataset object
  object  = new(
    "comradesDataSet",
    rnas = rnas,
    rnaSize = rnaSize2,
    sampleTable = sampleTable,
    hybFiles = hybFiles,
    matrixList = matrixList,
    clusterGrangesList = list(),
    clusterTableList = list(),
    clusterTableFolded = data.frame(),
    interactionTable = data.frame(),
    viennaStructures = list(),
    dgs = list()
  )
  
  return(object)
  
  
}
