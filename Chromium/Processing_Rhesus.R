# Load Raw Counts + pre-process (calculate percent.mito, filter empty droplets, exclude doublets) then merge sample objects)
library(Seurat)
library(CellMembrane)

--------------------------------------------------------
# Load Raw Counts + pre-process (calculate percent.mito, filter empty droplets, exclude doublets)
# Conditions
emptyDropsLower <- 200
emptyDropsFdrThreshold <- 0.001
maxAllowableCells <- 40000
minAllowableCells <- 500
useEmptyDropsCellRanger <- TRUE
nExpectedCells <- 10000
useCellBender <- FALSE
useSoupX <- FALSE

outputPrefix <- 'SingleCell.rawCounts'

seuratObjects <- list()
datasetIdToName <- list()
datasetIdToOutputFileId<- list()
datasetIdToReadset<- list()

seuratObjects[['521595']] <- '521595_RawData'
datasetIdToName[['521595']] <- 'ONPRC025_XX_D52_B_S2'
datasetIdToOutputFileId[['521595']] <- 521595
datasetIdToReadset[['521595']] <- 92312
seuratObjects[['521591']] <- '521591_RawData'
datasetIdToName[['521591']] <- 'ONPRC017_XX_D130_B_S10'
datasetIdToOutputFileId[['521591']] <- 521591
datasetIdToReadset[['521591']] <- 84652
seuratObjects[['521588']] <- '521588_RawData'
datasetIdToName[['521588']] <- 'ONPRC027_XX_D50_B_S7'
datasetIdToOutputFileId[['521588']] <- 521588
datasetIdToReadset[['521588']] <- 92310
seuratObjects[['521586']] <- '521586_RawData'
datasetIdToName[['521586']] <- 'ONPRC018_XX_D50_A_S1'
datasetIdToOutputFileId[['521586']] <- 521586
datasetIdToReadset[['521586']] <- 84655
seuratObjects[['521583']] <- '521583_RawData'
datasetIdToName[['521583']] <- 'ONPRC016_XX_D130_B_S8'
datasetIdToOutputFileId[['521583']] <- 521583
datasetIdToReadset[['521583']] <- 84651
seuratObjects[['521579']] <- '521579_RawData'
datasetIdToName[['521579']] <- 'ONPRC017_XX_D130_A_S9'
datasetIdToOutputFileId[['521579']] <- 521579
datasetIdToReadset[['521579']] <- 84648
seuratObjects[['521576']] <- '521576_RawData'
datasetIdToName[['521576']] <- 'ONPRC025_XX_D52_A_S1'
datasetIdToOutputFileId[['521576']] <- 521576
datasetIdToReadset[['521576']] <- 92309
seuratObjects[['521574']] <- '521574_RawData'
datasetIdToName[['521574']] <- 'ONPRC027_XX_D50_A_S6'
datasetIdToOutputFileId[['521574']] <- 521574
datasetIdToReadset[['521574']] <- 92311
seuratObjects[['521570']] <- '521570_RawData'
datasetIdToName[['521570']] <- ‘ONPRC012_XX_D100_B_S8'
datasetIdToOutputFileId[['521570']] <- 521570
datasetIdToReadset[['521570']] <- 84650
seuratObjects[['521568']] <- '521568_RawData'
datasetIdToName[['521568']] <- 'ONPRC018_XX_D50_B_S2'
datasetIdToOutputFileId[['521568']] <- 521568
datasetIdToReadset[['521568']] <- 84654
seuratObjects[['521563']] <- '521563_RawData'
datasetIdToName[['521563']] <- ‘ONPRC012_XX_D100_A_S7'
datasetIdToOutputFileId[['521563']] <- 521563
datasetIdToReadset[['521563']] <- 84653
seuratObjects[['521560']] <- '521560_RawData'
datasetIdToName[['521560']] <- 'ONPRC016_XX_D130_A_S7'
datasetIdToOutputFileId[['521560']] <- 521560
datasetIdToReadset[['521560']] <- 84646
seuratObjects[['521556']] <- '521556_RawData'
datasetIdToName[['521556']] <- ‘ONPRC001_XX_D130_A_S1'
datasetIdToOutputFileId[['521556']] <- 521556
datasetIdToReadset[['521556']] <- 84649
seuratObjects[['521554']] <- '521554_RawData'
datasetIdToName[['521554']] <- ‘ONPRC001_XX_D130_B_S2'
datasetIdToOutputFileId[['521554']] <- 521554
datasetIdToReadset[['521554']] <- 84647
seuratObjects[['521550']] <- '521550_RawData'
datasetIdToName[['521550']] <- ‘ONPRC006_XX_D100_B_S11'
datasetIdToOutputFileId[['521550']] <- 521550
datasetIdToReadset[['521550']] <- 84645
seuratObjects[['521548']] <- '521548_RawData'
datasetIdToName[['521548']] <- ‘ONPRC007_XX_D100_A_S8'
datasetIdToOutputFileId[['521548']] <- 521548
datasetIdToReadset[['521548']] <- 84643
seuratObjects[['521545']] <- '521545_RawData'
datasetIdToName[['521545']] <- ‘ONPRC006_XX_D100_A_S10'
datasetIdToOutputFileId[['521545']] <- 521545
datasetIdToReadset[['521545']] <- 84644
seuratObjects[['521535']] <- '521535_RawData'
datasetIdToName[['521535']] <- ‘ONPRC007_XX_D100_B_S9'
datasetIdToOutputFileId[['521535']] <- 521535
datasetIdToReadset[['521535']] <- 84642

# Binds arguments from the environment to the target function
bindArgs <- function(fun, seuratObj, allowableArgNames = NULL, disallowedArgNames = NULL) {
  boundArgs <- list()
  boundArgs[['seuratObj']] <- seuratObj
  
  for (name in names(formals(fun))) {
    if (!is.null(disallowedArgNames) && (name %in% disallowedArgNames)) {
      next
    }
    else if (name %in% names(boundArgs)) {
      next
    }
    else if (exists(name)) {
      if (!is.null(allowableArgNames) && !(name %in% allowableArgNames)) {
        next
      }
      
      val <- get(name)
      displayVal <- val
      if (all(is.null(val))) {
        displayVal <- 'NULL'
      } else if (all(is.na(val))) {
        displayVal <- 'NA'
      } else if (is.object(val)) {
        displayVal <- '[object]'
      }
      
      if (length(displayVal) > 1) {
        displayVal <- paste0(displayVal, collapse = ',')
      }
      
      print(paste0('Binding argument: ', name, ': ', displayVal))
      if (all(is.null(val))) {
        boundArgs[name] <- list(NULL)
      } else {
        boundArgs[[name]] <- val
      }
    }
  }
  
  formals(fun)[names(boundArgs)] <- boundArgs
  
  fun
}

clearSeuratCommands <- function(seuratObj, maxSize = 500000) {
  for (commandName in names(seuratObj@commands)) {
    val <- object.size(x = slot(seuratObj@commands[[commandName]], 'call.string'))
    if (val > maxSize) {
      print(paste0('Clearing call.string for: ', commandName, '. size: ', format(val, units = 'auto')))
      slot(seuratObj@commands[[commandName]], 'call.string') <- ''
    }
  }
  
  return(seuratObj)
}

printName <- function(datasetId) {
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  print(paste0('Processing dataset: ', datasetName))
}

savedFiles <- data.frame(datasetId = character(), datasetName = character(), filename = character(), outputFileId = character(), readsetId = character())
if (file.exists('/work/savedSeuratObjects.txt')) {
  print('Deleting pre-existing savedSeuratObjects.txt file')
  unlink('/work/savedSeuratObjects.txt')
}

file.create('/work/savedSeuratObjects.txt')
## [1] TRUE
print(paste0('Total lines in savedSeuratObjects.txt on job start:', length(readLines('savedSeuratObjects.txt'))))
## [1] "Total lines in savedSeuratObjects.txt on job start:0"
saveData <- function(seuratObj, datasetId) {
  print(paste0('Saving dataset: ', datasetId))
  print(seuratObj)
  
  datasetIdForFile <- makeLegalFileName(datasetId)
  fn <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.rds')
  
  message(paste0('Saving RDS file: ', fn, ' with ', ncol(seuratObj), ' cells'))
  barcodeFile <- paste0(outputPrefix, '.', datasetIdForFile, '.cellBarcodes.csv')
  metaFile <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.meta.txt')
  
  saveRDS(seuratObj, file = fn)
  
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  
  # NOTE: this is the ID of the original loupe file. Needed for operations like appending cell hashing or CITE-seq
  outputFileId <- ifelse(datasetId %in% names(datasetIdToOutputFileId), yes = datasetIdToOutputFileId[[datasetId]], no = NA)
  
  readsetId <- ifelse(datasetId %in% names(datasetIdToReadset), yes = datasetIdToReadset[[datasetId]], no = NA)
  print(paste0('readsetId: ', readsetId))
  
  toAppend <- data.frame(datasetId = datasetId, datasetName = datasetName, filename = fn, outputFileId = outputFileId, readsetId = readsetId)
  if (nrow(toAppend) != 1) {
    warning(paste0('Error saving seurat objects, more than one row:'))
    print(toAppend)
    stop('Error saving seurat objects, more than one row!')
  }
  
  write.table(toAppend, file = 'savedSeuratObjects.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)
  print(paste0('Total lines in savedSeuratObjects.txt after save:', length(readLines('savedSeuratObjects.txt'))))
  
# Write cell barcodes and metadata:
metaDf <- seuratObj@meta.data
metaDf$cellbarcode <- colnames(seuratObj)
write.table(metaDf, file = metaFile, quote = T, row.names = F, sep = ',', col.names = T)
write.table(data.frame(CellBarcode = colnames(seuratObj)), file = barcodeFile, quote = F, row.names = F, sep = ',', col.names = F)
}

intermediateFiles <- c()
addIntermediateFile <- function(f) {
  print(paste0('Adding intermediate file: ', f))
  intermediateFiles <<- c(intermediateFiles, f)
}

makeLegalFileName <- function(fn) {
  fn <- gsub(fn, pattern = '\\\\', replacement = '_')
  fn <- gsub(fn, pattern = '[/ ,;]', replacement = '_')
  fn <- gsub(fn, pattern = '\\|', replacement = '_')
  return(fn)
}

errorMessages <- c()
addErrorMessage <- function(f) {
  print(paste0('Adding error: ', f))
  errorMessages <<- c(errorMessages, f)
}

print('Updating future.globals.maxSize')
## [1] "Updating future.globals.maxSize"
options(future.globals.maxSize = Inf)

options('Seurat.memsafe' = TRUE)

if (Sys.getenv('SEURAT_MAX_THREADS') != '') {
  print(paste0('Setting future::plan workers to: ', Sys.getenv('SEURAT_MAX_THREADS')))
  future::plan(strategy='multisession', workers=Sys.getenv('SEURAT_MAX_THREADS'))
}


# Load Raw Counts
for (datasetId in names(seuratObjects)) {
  printName(datasetId)
  rawCountDir <- seuratObjects[[datasetId]]
  
  datasetName <- datasetIdToName[[datasetId]]
  
  previouslyFilteredMatrix <- NULL
  if (useCellBender) {
    # This is the 10x project:
    print('Will use cellbender-adjusted counts instead of raw counts:')
    previouslyFilteredMatrix <- paste0(rawCountDir, '/raw_feature_bc_matrix.cellbender_filtered.h5')
    if (!file.exists(previouslyFilteredMatrix)) {
      stop(paste0('Unable to find file: ', previouslyFilteredMatrix, '. You can re-run cellbender to fix this.'))
    }
  }
  
seuratObj <- CellMembrane::ReadAndFilter10xData(
  dataDir = rawCountDir, 
  datasetId = datasetId, 
  datasetName = datasetName, 
  emptyDropsLower = emptyDropsLower, 
  emptyDropsFdrThreshold = emptyDropsFdrThreshold, 
  useEmptyDropsCellRanger = useEmptyDropsCellRanger, 
  nExpectedCells = nExpectedCells, 
  useSoupX = useSoupX, 
  previouslyFilteredMatrix = previouslyFilteredMatrix)
  
  if (!is.null(maxAllowableCells) && ncol(seuratObj) > maxAllowableCells) {
    addErrorMessage(paste0('The seurat object has ', ncol(seuratObj), ' cells, which is more than the max allowable cells (', maxAllowableCells, '). Please review emptyDrops results as this probably means thresholds were suboptimal.'))
  }
  
  if (!is.null(minAllowableCells) && ncol(seuratObj) < minAllowableCells) {
    addErrorMessage(paste0('The seurat object has ', ncol(seuratObj), ' cells, which is less than the min allowable cells (', minAllowableCells, '). Please review emptyDrops results as this probably means thresholds were suboptimal.'))
  }
  
  saveData(seuratObj, datasetId)
  
# Cleanup
rm(seuratObj)
gc()
}

## Note: CellMembrane would normally run mitoGenesPattern = "^MT-", but "There were no genes matching mitoGenesPattern, so attempting to identify MT genes using name"
## "Total mito features: 13"
## Example code (don't need to run, embedded)
CalculatePercentMito <- function(seuratObj, mitoGenesPattern = "^MT-", annotateMitoFromReferenceIfNoHitsFound = TRUE, outputColName = 'p.mito') {
  mito.features <- grep(pattern = mitoGenesPattern, x = rownames(x = seuratObj), value = TRUE)
  if (length(mito.features) == 0 && annotateMitoFromReferenceIfNoHitsFound) {
    print('There were no genes matching mitoGenesPattern, so attempting to identify MT genes using name')
    mito.features <- c('ATP6','ATP8','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND4L','ND5','ND6')
    mito.features.prefix <- c('MT-', mito.features)
    i1 <- length(intersect(mito.features, rownames(seuratObj)))
    i2 <- length(intersect(mito.features.prefix, rownames(seuratObj)))
    if (i2 > 0 && i2 > i1) {
      print('Selecting using reference gene set with MT- prefix')
      mito.features <- mito.features.prefix
    }
    
    mito.features <- intersect(mito.features, rownames(seuratObj))
  }
  
  print(paste0('Total mito features: ', length(mito.features)))
  mito.features <- intersect(mito.features, rownames(seuratObj))
  print(paste0('Total intersecting with seurat rownames (total: ', length(rownames(seuratObj)),'): ', length(mito.features)))
  
  if (all(is.null(mito.features)) || length(mito.features) == 0) {
    print('No mito features found')
    seuratObj[[outputColName]] <- 0
  } else {
    p.mito <- Matrix::colSums(x = GetAssayData(object = seuratObj, layer = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObj, layer = 'counts'))
    seuratObj[[outputColName]] <- p.mito
  }
  
  return(seuratObj)

--------
# Example output info

## [1] "Saving dataset: 521595"
## An object of class Seurat 
## 35432 features across 13677 samples within 1 assay 
## Active assay: RNA (35432 features, 0 variable features)
## [1] "readsetId: 92312"
## [1] "Total lines in savedSeuratObjects.txt after save:1"
## [1] "Processing dataset: ONPRC017_XX_D130_B_S10"
## [1] "Loading counts and running emptyDrops"

## [1] "Knee: 3429"
## [1] "Inflection: 2123"
## [1] "Performing emptyDropsCellRanger with 10000 iterations"
## [1] "Input cells: 1831254"
## [1] "Cells with non-finite FDR: 1816571"
## [1] "Cells passing FDR (0.001): 13204"
## [1] "Cells failing FDR: 1479"

## [1] "Total rows with non-finite probabilities: 8530"
## [1] "Min UMI count in a droplet called a cell: 501"
## [1] "Max UMI count in a droplet not called a cell: 1955"
## [1] "Replacing underscores with hyphens in feature names"
## [1] "Adding barcode prefix: 521591"
## [1] "There were no genes matching mitoGenesPattern, so attempting to identify MT genes using name"
## [1] "Total mito features: 13"
## [1] "Total intersecting with seurat rownames (total: 35432): 13"
## [1] "Total lines in savedSeuratObjects.txt after save:18"

if (length(intermediateFiles) > 0) {
  write.table(data.frame(file = intermediateFiles), file = 'intermediateFiles.txt', quote = FALSE, delim = '\t', row.names = FALSE, col.names = FALSE)
}

if (length(errorMessages) > 0) {
  print('There were errors:')
  for (msg in errorMessages) {
    print(msg)
  }
  
  write(errorMessages, file = 'seuratErrors.txt')
}

if (file.exists('savedSeuratObjects.txt')) {
  print(paste0('Total lines in savedSeuratObjects.txt:', length(readLines('savedSeuratObjects.txt'))))
} else {
  print('File does not exist: savedSeuratObjects.txt')
}
## [1] "Total lines in savedSeuratObjects.txt:18"

  
------------------------------------------
# Exclude doublets - CellMembrane uses scDblFinder
# Conditions
dropDoublets <- TRUE

outputPrefix <- 'SingleCell.df'

seuratObjects <- list()
datasetIdToName <- list()
datasetIdToOutputFileId<- list()
datasetIdToReadset<- list()
seuratObjects[['521595']] <- 'SingleCell.rawCounts.521595.seurat.rds'
datasetIdToName[['521595']] <- 'ONPRC025_XX_D52_B_S2'
datasetIdToOutputFileId[['521595']] <- 521595
datasetIdToReadset[['521595']] <- 92312
seuratObjects[['521591']] <- 'SingleCell.rawCounts.521591.seurat.rds'
datasetIdToName[['521591']] <- 'ONPRC017_XX_D130_B_S10'
datasetIdToOutputFileId[['521591']] <- 521591
datasetIdToReadset[['521591']] <- 84652
seuratObjects[['521588']] <- 'SingleCell.rawCounts.521588.seurat.rds'
datasetIdToName[['521588']] <- 'ONPRC027_XX_D50_B_S7'
datasetIdToOutputFileId[['521588']] <- 521588
datasetIdToReadset[['521588']] <- 92310
seuratObjects[['521586']] <- 'SingleCell.rawCounts.521586.seurat.rds'
datasetIdToName[['521586']] <- 'ONPRC018_XX_D50_A_S1'
datasetIdToOutputFileId[['521586']] <- 521586
datasetIdToReadset[['521586']] <- 84655
seuratObjects[['521583']] <- 'SingleCell.rawCounts.521583.seurat.rds'
datasetIdToName[['521583']] <- 'ONPRC016_XX_D130_B_S8'
datasetIdToOutputFileId[['521583']] <- 521583
datasetIdToReadset[['521583']] <- 84651
seuratObjects[['521579']] <- 'SingleCell.rawCounts.521579.seurat.rds'
datasetIdToName[['521579']] <- 'ONPRC017_XX_D130_A_S9'
datasetIdToOutputFileId[['521579']] <- 521579
datasetIdToReadset[['521579']] <- 84648
seuratObjects[['521576']] <- 'SingleCell.rawCounts.521576.seurat.rds'
datasetIdToName[['521576']] <- 'ONPRC025_XX_D52_A_S1'
datasetIdToOutputFileId[['521576']] <- 521576
datasetIdToReadset[['521576']] <- 92309
seuratObjects[['521574']] <- 'SingleCell.rawCounts.521574.seurat.rds'
datasetIdToName[['521574']] <- 'ONPRC027_XX_D50_A_S6'
datasetIdToOutputFileId[['521574']] <- 521574
datasetIdToReadset[['521574']] <- 92311
seuratObjects[['521570']] <- 'SingleCell.rawCounts.521570.seurat.rds'
datasetIdToName[['521570']] <- 'ONPRC012_XX_D100_B_S8'
datasetIdToOutputFileId[['521570']] <- 521570
datasetIdToReadset[['521570']] <- 84650
seuratObjects[['521568']] <- 'SingleCell.rawCounts.521568.seurat.rds'
datasetIdToName[['521568']] <- 'ONPRC018_XX_D50_B_S2'
datasetIdToOutputFileId[['521568']] <- 521568
datasetIdToReadset[['521568']] <- 84654
seuratObjects[['521563']] <- 'SingleCell.rawCounts.521563.seurat.rds'
datasetIdToName[['521563']] <- 'ONPRC012_XX_D100_A_S7'
datasetIdToOutputFileId[['521563']] <- 521563
datasetIdToReadset[['521563']] <- 84653
seuratObjects[['521560']] <- 'SingleCell.rawCounts.521560.seurat.rds'
datasetIdToName[['521560']] <- 'ONPRC016_XX_D130_A_S7'
datasetIdToOutputFileId[['521560']] <- 521560
datasetIdToReadset[['521560']] <- 84646
seuratObjects[['521556']] <- 'SingleCell.rawCounts.521556.seurat.rds'
datasetIdToName[['521556']] <- 'ONPRC001_XX_D130_A_S1'
datasetIdToOutputFileId[['521556']] <- 521556
datasetIdToReadset[['521556']] <- 84649
seuratObjects[['521554']] <- 'SingleCell.rawCounts.521554.seurat.rds'
datasetIdToName[['521554']] <- 'ONPRC001_XX_D130_B_S2'
datasetIdToOutputFileId[['521554']] <- 521554
datasetIdToReadset[['521554']] <- 84647
seuratObjects[['521550']] <- 'SingleCell.rawCounts.521550.seurat.rds'
datasetIdToName[['521550']] <- 'ONPRC006_XX_D100_B_S11'
datasetIdToOutputFileId[['521550']] <- 521550
datasetIdToReadset[['521550']] <- 84645
seuratObjects[['521548']] <- 'SingleCell.rawCounts.521548.seurat.rds'
datasetIdToName[['521548']] <- 'ONPRC007_XX_D100_A_S8'
datasetIdToOutputFileId[['521548']] <- 521548
datasetIdToReadset[['521548']] <- 84643
seuratObjects[['521545']] <- 'SingleCell.rawCounts.521545.seurat.rds'
datasetIdToName[['521545']] <- 'ONPRC006_XX_D100_A_S10'
datasetIdToOutputFileId[['521545']] <- 521545
datasetIdToReadset[['521545']] <- 84644
seuratObjects[['521535']] <- 'SingleCell.rawCounts.521535.seurat.rds'
datasetIdToName[['521535']] <- 'ONPRC007_XX_D100_B_S9'
datasetIdToOutputFileId[['521535']] <- 521535
datasetIdToReadset[['521535']] <- 84642

# Binds arguments from the environment to the target function
bindArgs <- function(fun, seuratObj, allowableArgNames = NULL, disallowedArgNames = NULL) {
  boundArgs <- list()
  boundArgs[['seuratObj']] <- seuratObj
  
  for (name in names(formals(fun))) {
    if (!is.null(disallowedArgNames) && (name %in% disallowedArgNames)) {
      next
    }
    else if (name %in% names(boundArgs)) {
      next
    }
    else if (exists(name)) {
      if (!is.null(allowableArgNames) && !(name %in% allowableArgNames)) {
        next
      }
      
      val <- get(name)
      displayVal <- val
      if (all(is.null(val))) {
        displayVal <- 'NULL'
      } else if (all(is.na(val))) {
        displayVal <- 'NA'
      } else if (is.object(val)) {
        displayVal <- '[object]'
      }
      
      if (length(displayVal) > 1) {
        displayVal <- paste0(displayVal, collapse = ',')
      }
      
      print(paste0('Binding argument: ', name, ': ', displayVal))
      if (all(is.null(val))) {
        boundArgs[name] <- list(NULL)
      } else {
        boundArgs[[name]] <- val
      }
    }
  }
  
  formals(fun)[names(boundArgs)] <- boundArgs
  
  fun
}

clearSeuratCommands <- function(seuratObj, maxSize = 500000) {
  for (commandName in names(seuratObj@commands)) {
    val <- object.size(x = slot(seuratObj@commands[[commandName]], 'call.string'))
    if (val > maxSize) {
      print(paste0('Clearing call.string for: ', commandName, '. size: ', format(val, units = 'auto')))
      slot(seuratObj@commands[[commandName]], 'call.string') <- ''
    }
  }
  
  return(seuratObj)
}

printName <- function(datasetId) {
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  print(paste0('Processing dataset: ', datasetName))
}

savedFiles <- data.frame(datasetId = character(), datasetName = character(), filename = character(), outputFileId = character(), readsetId = character())
if (file.exists('/work/savedSeuratObjects.txt')) {
  print('Deleting pre-existing savedSeuratObjects.txt file')
  unlink('/work/savedSeuratObjects.txt')
}

file.create('/work/savedSeuratObjects.txt')
## [1] TRUE
print(paste0('Total lines in savedSeuratObjects.txt on job start:', length(readLines('savedSeuratObjects.txt'))))
## [1] "Total lines in savedSeuratObjects.txt on job start:0"
saveData <- function(seuratObj, datasetId) {
  print(paste0('Saving dataset: ', datasetId))
  print(seuratObj)
  
  datasetIdForFile <- makeLegalFileName(datasetId)
  fn <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.rds')
  
  message(paste0('Saving RDS file: ', fn, ' with ', ncol(seuratObj), ' cells'))
  barcodeFile <- paste0(outputPrefix, '.', datasetIdForFile, '.cellBarcodes.csv')
  metaFile <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.meta.txt')
  
  saveRDS(seuratObj, file = fn)
  
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  
  # NOTE: this is the ID of the original loupe file. Needed for operations like appending cell hashing or CITE-seq
  outputFileId <- ifelse(datasetId %in% names(datasetIdToOutputFileId), yes = datasetIdToOutputFileId[[datasetId]], no = NA)
  
  readsetId <- ifelse(datasetId %in% names(datasetIdToReadset), yes = datasetIdToReadset[[datasetId]], no = NA)
  print(paste0('readsetId: ', readsetId))
  
  toAppend <- data.frame(datasetId = datasetId, datasetName = datasetName, filename = fn, outputFileId = outputFileId, readsetId = readsetId)
  if (nrow(toAppend) != 1) {
    warning(paste0('Error saving seurat objects, more than one row:'))
    print(toAppend)
    stop('Error saving seurat objects, more than one row!')
  }
  
  write.table(toAppend, file = 'savedSeuratObjects.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)
  print(paste0('Total lines in savedSeuratObjects.txt after save:', length(readLines('savedSeuratObjects.txt'))))
  
# Write cell barcodes and metadata:
metaDf <- seuratObj@meta.data
metaDf$cellbarcode <- colnames(seuratObj)
write.table(metaDf, file = metaFile, quote = T, row.names = F, sep = ',', col.names = T)
write.table(data.frame(CellBarcode = colnames(seuratObj)), file = barcodeFile, quote = F, row.names = F, sep = ',', col.names = F)
}

intermediateFiles <- c()
addIntermediateFile <- function(f) {
  print(paste0('Adding intermediate file: ', f))
  intermediateFiles <<- c(intermediateFiles, f)
}

makeLegalFileName <- function(fn) {
  fn <- gsub(fn, pattern = '\\\\', replacement = '_')
  fn <- gsub(fn, pattern = '[/ ,;]', replacement = '_')
  fn <- gsub(fn, pattern = '\\|', replacement = '_')
  return(fn)
}

errorMessages <- c()
addErrorMessage <- function(f) {
  print(paste0('Adding error: ', f))
  errorMessages <<- c(errorMessages, f)
}

print('Updating future.globals.maxSize')
## [1] "Updating future.globals.maxSize"
options(future.globals.maxSize = Inf)

options('Seurat.memsafe' = TRUE)

if (Sys.getenv('SEURAT_MAX_THREADS') != '') {
  print(paste0('Setting future::plan workers to: ', Sys.getenv('SEURAT_MAX_THREADS')))
  future::plan(strategy='multisession', workers=Sys.getenv('SEURAT_MAX_THREADS'))
}

# Exclude doublets - CellMembrane uses scDblFinder
DoubletFinder
for (datasetId in names(seuratObjects)) {
  printName(datasetId)
  seuratObj <- readRDS(seuratObjects[[datasetId]])
  
  seuratObj <- CellMembrane::FindDoublets(seuratObj, dropDoublets = dropDoublets)
  
  saveData(seuratObj, datasetId)
  
# Cleanup
rm(seuratObj)
gc()
}

# Example output info
## [1] "Processing dataset: ONPRC025_XX_D52_B_S2"
## [1] "Finding doublets using scDblFinder"

## [1] "No reductions calculated, cannot plot tSNE/UMAP"
## [1] "Dropping 2615 doublets"
## [1] "Saving dataset: 521595"
## An object of class Seurat 
## 35432 features across 11062 samples within 1 assay 
## Active assay: RNA (35432 features, 0 variable features)
## [1] "readsetId: 92312"
## [1] "Total lines in savedSeuratObjects.txt after save:1"

## [1] "Total lines in savedSeuratObjects.txt after save:18"
if (length(intermediateFiles) > 0) {
  write.table(data.frame(file = intermediateFiles), file = 'intermediateFiles.txt', quote = FALSE, delim = '\t', row.names = FALSE, col.names = FALSE)
}

if (length(errorMessages) > 0) {
  print('There were errors:')
  for (msg in errorMessages) {
    print(msg)
  }
  
  write(errorMessages, file = 'seuratErrors.txt')
}

if (file.exists('savedSeuratObjects.txt')) {
  print(paste0('Total lines in savedSeuratObjects.txt:', length(readLines('savedSeuratObjects.txt'))))
} else {
  print('File does not exist: savedSeuratObjects.txt')
}
## [1] "Total lines in savedSeuratObjects.txt:18"

-------------------------------------------------
# Filter Raw Counts
# Conditions
nCount_RNA.low <- 0
nCount_RNA.high <- 30000
nFeature.low <- 200
nFeature.high <- 6000
pMito.low <- 0
pMito.high <- 0.4

outputPrefix <- 'SingleCell.df.frc'

seuratObjects <- list()
datasetIdToName <- list()
datasetIdToOutputFileId<- list()
datasetIdToReadset<- list()
seuratObjects[['521595']] <- 'SingleCell.df.521595.seurat.rds'
datasetIdToName[['521595']] <- 'ONPRC025_XX_D52_B_S2'
datasetIdToOutputFileId[['521595']] <- 521595
datasetIdToReadset[['521595']] <- 92312
seuratObjects[['521591']] <- 'SingleCell.df.521591.seurat.rds'
datasetIdToName[['521591']] <- 'ONPRC017_XX_D130_B_S10'
datasetIdToOutputFileId[['521591']] <- 521591
datasetIdToReadset[['521591']] <- 84652
seuratObjects[['521588']] <- 'SingleCell.df.521588.seurat.rds'
datasetIdToName[['521588']] <- 'ONPRC027_XX_D50_B_S7'
datasetIdToOutputFileId[['521588']] <- 521588
datasetIdToReadset[['521588']] <- 92310
seuratObjects[['521586']] <- 'SingleCell.df.521586.seurat.rds'
datasetIdToName[['521586']] <- 'ONPRC018_XX_D50_A_S1'
datasetIdToOutputFileId[['521586']] <- 521586
datasetIdToReadset[['521586']] <- 84655
seuratObjects[['521583']] <- 'SingleCell.df.521583.seurat.rds'
datasetIdToName[['521583']] <- 'ONPRC016_XX_D130_B_S8'
datasetIdToOutputFileId[['521583']] <- 521583
datasetIdToReadset[['521583']] <- 84651
seuratObjects[['521579']] <- 'SingleCell.df.521579.seurat.rds'
datasetIdToName[['521579']] <- 'ONPRC017_XX_D130_A_S9'
datasetIdToOutputFileId[['521579']] <- 521579
datasetIdToReadset[['521579']] <- 84648
seuratObjects[['521576']] <- 'SingleCell.df.521576.seurat.rds'
datasetIdToName[['521576']] <- 'ONPRC025_XX_D52_A_S1'
datasetIdToOutputFileId[['521576']] <- 521576
datasetIdToReadset[['521576']] <- 92309
seuratObjects[['521574']] <- 'SingleCell.df.521574.seurat.rds'
datasetIdToName[['521574']] <- 'ONPRC027_XX_D50_A_S6'
datasetIdToOutputFileId[['521574']] <- 521574
datasetIdToReadset[['521574']] <- 92311
seuratObjects[['521570']] <- 'SingleCell.df.521570.seurat.rds'
datasetIdToName[['521570']] <- 'ONPRC012_XX_D100_B_S8'
datasetIdToOutputFileId[['521570']] <- 521570
datasetIdToReadset[['521570']] <- 84650
seuratObjects[['521568']] <- 'SingleCell.df.521568.seurat.rds'
datasetIdToName[['521568']] <- 'ONPRC018_XX_D50_B_S2'
datasetIdToOutputFileId[['521568']] <- 521568
datasetIdToReadset[['521568']] <- 84654
seuratObjects[['521563']] <- 'SingleCell.df.521563.seurat.rds'
datasetIdToName[['521563']] <- 'ONPRC012_XX_D100_A_S7'
datasetIdToOutputFileId[['521563']] <- 521563
datasetIdToReadset[['521563']] <- 84653
seuratObjects[['521560']] <- 'SingleCell.df.521560.seurat.rds'
datasetIdToName[['521560']] <- 'ONPRC016_XX_D130_A_S7'
datasetIdToOutputFileId[['521560']] <- 521560
datasetIdToReadset[['521560']] <- 84646
seuratObjects[['521556']] <- 'SingleCell.df.521556.seurat.rds'
datasetIdToName[['521556']] <- 'ONPRC001_XX_D130_A_S1'
datasetIdToOutputFileId[['521556']] <- 521556
datasetIdToReadset[['521556']] <- 84649
seuratObjects[['521554']] <- 'SingleCell.df.521554.seurat.rds'
datasetIdToName[['521554']] <- 'ONPRC001_XX_D130_B_S2'
datasetIdToOutputFileId[['521554']] <- 521554
datasetIdToReadset[['521554']] <- 84647
seuratObjects[['521550']] <- 'SingleCell.df.521550.seurat.rds'
datasetIdToName[['521550']] <- 'ONPRC006_XX_D100_B_S11'
datasetIdToOutputFileId[['521550']] <- 521550
datasetIdToReadset[['521550']] <- 84645
seuratObjects[['521548']] <- 'SingleCell.df.521548.seurat.rds'
datasetIdToName[['521548']] <- 'ONPRC007_XX_D100_A_S8'
datasetIdToOutputFileId[['521548']] <- 521548
datasetIdToReadset[['521548']] <- 84643
seuratObjects[['521545']] <- 'SingleCell.df.521545.seurat.rds'
datasetIdToName[['521545']] <- 'ONPRC006_XX_D100_A_S10'
datasetIdToOutputFileId[['521545']] <- 521545
datasetIdToReadset[['521545']] <- 84644
seuratObjects[['521535']] <- 'SingleCell.df.521535.seurat.rds'
datasetIdToName[['521535']] <- 'ONPRC007_XX_D100_B_S9'
datasetIdToOutputFileId[['521535']] <- 521535
datasetIdToReadset[['521535']] <- 84642

# Binds arguments from the environment to the target function
bindArgs <- function(fun, seuratObj, allowableArgNames = NULL, disallowedArgNames = NULL) {
  boundArgs <- list()
  boundArgs[['seuratObj']] <- seuratObj
  
  for (name in names(formals(fun))) {
    if (!is.null(disallowedArgNames) && (name %in% disallowedArgNames)) {
      next
    }
    else if (name %in% names(boundArgs)) {
      next
    }
    else if (exists(name)) {
      if (!is.null(allowableArgNames) && !(name %in% allowableArgNames)) {
        next
      }
      
      val <- get(name)
      displayVal <- val
      if (all(is.null(val))) {
        displayVal <- 'NULL'
      } else if (all(is.na(val))) {
        displayVal <- 'NA'
      } else if (is.object(val)) {
        displayVal <- '[object]'
      }
      
      if (length(displayVal) > 1) {
        displayVal <- paste0(displayVal, collapse = ',')
      }
      
      print(paste0('Binding argument: ', name, ': ', displayVal))
      if (all(is.null(val))) {
        boundArgs[name] <- list(NULL)
      } else {
        boundArgs[[name]] <- val
      }
    }
  }
  
  formals(fun)[names(boundArgs)] <- boundArgs
  
  fun
}

clearSeuratCommands <- function(seuratObj, maxSize = 500000) {
  for (commandName in names(seuratObj@commands)) {
    val <- object.size(x = slot(seuratObj@commands[[commandName]], 'call.string'))
    if (val > maxSize) {
      print(paste0('Clearing call.string for: ', commandName, '. size: ', format(val, units = 'auto')))
      slot(seuratObj@commands[[commandName]], 'call.string') <- ''
    }
  }
  
  return(seuratObj)
}

printName <- function(datasetId) {
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  print(paste0('Processing dataset: ', datasetName))
}

savedFiles <- data.frame(datasetId = character(), datasetName = character(), filename = character(), outputFileId = character(), readsetId = character())
if (file.exists('/work/savedSeuratObjects.txt')) {
  print('Deleting pre-existing savedSeuratObjects.txt file')
  unlink('/work/savedSeuratObjects.txt')
}

file.create('/work/savedSeuratObjects.txt')
## [1] TRUE
print(paste0('Total lines in savedSeuratObjects.txt on job start:', length(readLines('savedSeuratObjects.txt'))))
## [1] "Total lines in savedSeuratObjects.txt on job start:0"
saveData <- function(seuratObj, datasetId) {
  print(paste0('Saving dataset: ', datasetId))
  print(seuratObj)
  
  datasetIdForFile <- makeLegalFileName(datasetId)
  fn <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.rds')
  
  message(paste0('Saving RDS file: ', fn, ' with ', ncol(seuratObj), ' cells'))
  barcodeFile <- paste0(outputPrefix, '.', datasetIdForFile, '.cellBarcodes.csv')
  metaFile <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.meta.txt')
  
  saveRDS(seuratObj, file = fn)
  
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  
  # NOTE: this is the ID of the original loupe file. Needed for operations like appending cell hashing or CITE-seq
  outputFileId <- ifelse(datasetId %in% names(datasetIdToOutputFileId), yes = datasetIdToOutputFileId[[datasetId]], no = NA)
  
  readsetId <- ifelse(datasetId %in% names(datasetIdToReadset), yes = datasetIdToReadset[[datasetId]], no = NA)
  print(paste0('readsetId: ', readsetId))
  
  toAppend <- data.frame(datasetId = datasetId, datasetName = datasetName, filename = fn, outputFileId = outputFileId, readsetId = readsetId)
  if (nrow(toAppend) != 1) {
    warning(paste0('Error saving seurat objects, more than one row:'))
    print(toAppend)
    stop('Error saving seurat objects, more than one row!')
  }
  
  write.table(toAppend, file = 'savedSeuratObjects.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)
  print(paste0('Total lines in savedSeuratObjects.txt after save:', length(readLines('savedSeuratObjects.txt'))))
  
# Write cell barcodes and metadata:
metaDf <- seuratObj@meta.data
metaDf$cellbarcode <- colnames(seuratObj)
write.table(metaDf, file = metaFile, quote = T, row.names = F, sep = ',', col.names = T)
write.table(data.frame(CellBarcode = colnames(seuratObj)), file = barcodeFile, quote = F, row.names = F, sep = ',', col.names = F)
}

intermediateFiles <- c()
addIntermediateFile <- function(f) {
  print(paste0('Adding intermediate file: ', f))
  intermediateFiles <<- c(intermediateFiles, f)
}

makeLegalFileName <- function(fn) {
  fn <- gsub(fn, pattern = '\\\\', replacement = '_')
  fn <- gsub(fn, pattern = '[/ ,;]', replacement = '_')
  fn <- gsub(fn, pattern = '\\|', replacement = '_')
  return(fn)
}

errorMessages <- c()
addErrorMessage <- function(f) {
  print(paste0('Adding error: ', f))
  errorMessages <<- c(errorMessages, f)
}

print('Updating future.globals.maxSize')
## [1] "Updating future.globals.maxSize"
options(future.globals.maxSize = Inf)

options('Seurat.memsafe' = TRUE)

if (Sys.getenv('SEURAT_MAX_THREADS') != '') {
  print(paste0('Setting future::plan workers to: ', Sys.getenv('SEURAT_MAX_THREADS')))
  future::plan(strategy='multisession', workers=Sys.getenv('SEURAT_MAX_THREADS'))
}

# Filter Raw Counts
Filter Raw Counts
for (datasetId in names(seuratObjects)) {
  printName(datasetId)
  seuratObj <- readRDS(seuratObjects[[datasetId]])
  
  seuratObj <- bindArgs(CellMembrane::FilterRawCounts, seuratObj)()
  
  saveData(seuratObj, datasetId)
  
# Cleanup
rm(seuratObj)
gc()
}

# Example output info
## [1] "Processing dataset: ONPRC025_XX_D52_B_S2"
## [1] "Binding argument: nCount_RNA.high: 30000"
## [1] "Binding argument: nCount_RNA.low: 0"
## [1] "Binding argument: nFeature.high: 6000"
## [1] "Binding argument: nFeature.low: 200"
## [1] "Binding argument: pMito.high: 0.4"
## [1] "Binding argument: pMito.low: 0"
## [1] "Filtering Cells..."
## [1] "Initial cells: 11062"

## [1] "After nCount_RNA filter: 10980"
## [1] "After nFeature_RNA filter: 10821"
## [1] "After p.mito filter: 10815"
## [1] "Final cells: 10815"
## [1] "Saving dataset: 521595"
## An object of class Seurat 
## 35432 features across 10815 samples within 1 assay 
## Active assay: RNA (35432 features, 0 variable features)
## [1] "readsetId: 92312"
## [1] "Total lines in savedSeuratObjects.txt after save:1"

## [1] "Total lines in savedSeuratObjects.txt after save:18"
if (length(intermediateFiles) > 0) {
  write.table(data.frame(file = intermediateFiles), file = 'intermediateFiles.txt', quote = FALSE, delim = '\t', row.names = FALSE, col.names = FALSE)
}

if (length(errorMessages) > 0) {
  print('There were errors:')
  for (msg in errorMessages) {
    print(msg)
  }
  
  write(errorMessages, file = 'seuratErrors.txt')
}

if (file.exists('savedSeuratObjects.txt')) {
  print(paste0('Total lines in savedSeuratObjects.txt:', length(readLines('savedSeuratObjects.txt'))))
} else {
  print('File does not exist: savedSeuratObjects.txt')
}
## [1] "Total lines in savedSeuratObjects.txt:18"


--------------------------------
# Append Cell Saturation
# Conditions
outputPrefix <- 'SingleCell.df.frc.saturation'

seuratObjects <- list()
datasetIdToName <- list()
datasetIdToOutputFileId<- list()
datasetIdToReadset<- list()
seuratObjects[['521595']] <- 'SingleCell.df.frc.521595.seurat.rds'
datasetIdToName[['521595']] <- 'ONPRC025_XX_D52_B_S2'
datasetIdToOutputFileId[['521595']] <- 521595
datasetIdToReadset[['521595']] <- 92312
seuratObjects[['521591']] <- 'SingleCell.df.frc.521591.seurat.rds'
datasetIdToName[['521591']] <- 'ONPRC017_XX_D130_B_S10'
datasetIdToOutputFileId[['521591']] <- 521591
datasetIdToReadset[['521591']] <- 84652
seuratObjects[['521588']] <- 'SingleCell.df.frc.521588.seurat.rds'
datasetIdToName[['521588']] <- 'ONPRC027_XX_D50_B_S7'
datasetIdToOutputFileId[['521588']] <- 521588
datasetIdToReadset[['521588']] <- 92310
seuratObjects[['521586']] <- 'SingleCell.df.frc.521586.seurat.rds'
datasetIdToName[['521586']] <- 'ONPRC018_XX_D50_A_S1'
datasetIdToOutputFileId[['521586']] <- 521586
datasetIdToReadset[['521586']] <- 84655
seuratObjects[['521583']] <- 'SingleCell.df.frc.521583.seurat.rds'
datasetIdToName[['521583']] <- 'ONPRC016_XX_D130_B_S8'
datasetIdToOutputFileId[['521583']] <- 521583
datasetIdToReadset[['521583']] <- 84651
seuratObjects[['521579']] <- 'SingleCell.df.frc.521579.seurat.rds'
datasetIdToName[['521579']] <- 'ONPRC017_XX_D130_A_S9'
datasetIdToOutputFileId[['521579']] <- 521579
datasetIdToReadset[['521579']] <- 84648
seuratObjects[['521576']] <- 'SingleCell.df.frc.521576.seurat.rds'
datasetIdToName[['521576']] <- 'ONPRC025_XX_D52_A_S1'
datasetIdToOutputFileId[['521576']] <- 521576
datasetIdToReadset[['521576']] <- 92309
seuratObjects[['521574']] <- 'SingleCell.df.frc.521574.seurat.rds'
datasetIdToName[['521574']] <- 'ONPRC027_XX_D50_A_S6'
datasetIdToOutputFileId[['521574']] <- 521574
datasetIdToReadset[['521574']] <- 92311
seuratObjects[['521570']] <- 'SingleCell.df.frc.521570.seurat.rds'
datasetIdToName[['521570']] <- 'ONPRC012_XX_D100_B_S8'
datasetIdToOutputFileId[['521570']] <- 521570
datasetIdToReadset[['521570']] <- 84650
seuratObjects[['521568']] <- 'SingleCell.df.frc.521568.seurat.rds'
datasetIdToName[['521568']] <- 'ONPRC018_XX_D50_B_S2'
datasetIdToOutputFileId[['521568']] <- 521568
datasetIdToReadset[['521568']] <- 84654
seuratObjects[['521563']] <- 'SingleCell.df.frc.521563.seurat.rds'
datasetIdToName[['521563']] <- 'ONPRC012_XX_D100_A_S7'
datasetIdToOutputFileId[['521563']] <- 521563
datasetIdToReadset[['521563']] <- 84653
seuratObjects[['521560']] <- 'SingleCell.df.frc.521560.seurat.rds'
datasetIdToName[['521560']] <- 'ONPRC016_XX_D130_A_S7'
datasetIdToOutputFileId[['521560']] <- 521560
datasetIdToReadset[['521560']] <- 84646
seuratObjects[['521556']] <- 'SingleCell.df.frc.521556.seurat.rds'
datasetIdToName[['521556']] <- 'ONPRC001_XX_D130_A_S1'
datasetIdToOutputFileId[['521556']] <- 521556
datasetIdToReadset[['521556']] <- 84649
seuratObjects[['521554']] <- 'SingleCell.df.frc.521554.seurat.rds'
datasetIdToName[['521554']] <- 'ONPRC001_XX_D130_B_S2'
datasetIdToOutputFileId[['521554']] <- 521554
datasetIdToReadset[['521554']] <- 84647
seuratObjects[['521550']] <- 'SingleCell.df.frc.521550.seurat.rds'
datasetIdToName[['521550']] <- 'ONPRC006_XX_D100_B_S11'
datasetIdToOutputFileId[['521550']] <- 521550
datasetIdToReadset[['521550']] <- 84645
seuratObjects[['521548']] <- 'SingleCell.df.frc.521548.seurat.rds'
datasetIdToName[['521548']] <- 'ONPRC007_XX_D100_A_S8'
datasetIdToOutputFileId[['521548']] <- 521548
datasetIdToReadset[['521548']] <- 84643
seuratObjects[['521545']] <- 'SingleCell.df.frc.521545.seurat.rds'
datasetIdToName[['521545']] <- 'ONPRC006_XX_D100_A_S10'
datasetIdToOutputFileId[['521545']] <- 521545
datasetIdToReadset[['521545']] <- 84644
seuratObjects[['521535']] <- 'SingleCell.df.frc.521535.seurat.rds'
datasetIdToName[['521535']] <- 'ONPRC007_XX_D100_B_S9'
datasetIdToOutputFileId[['521535']] <- 521535
datasetIdToReadset[['521535']] <- 84642

# Binds arguments from the environment to the target function
bindArgs <- function(fun, seuratObj, allowableArgNames = NULL, disallowedArgNames = NULL) {
  boundArgs <- list()
  boundArgs[['seuratObj']] <- seuratObj
  
  for (name in names(formals(fun))) {
    if (!is.null(disallowedArgNames) && (name %in% disallowedArgNames)) {
      next
    }
    else if (name %in% names(boundArgs)) {
      next
    }
    else if (exists(name)) {
      if (!is.null(allowableArgNames) && !(name %in% allowableArgNames)) {
        next
      }
      
      val <- get(name)
      displayVal <- val
      if (all(is.null(val))) {
        displayVal <- 'NULL'
      } else if (all(is.na(val))) {
        displayVal <- 'NA'
      } else if (is.object(val)) {
        displayVal <- '[object]'
      }
      
      if (length(displayVal) > 1) {
        displayVal <- paste0(displayVal, collapse = ',')
      }
      
      print(paste0('Binding argument: ', name, ': ', displayVal))
      if (all(is.null(val))) {
        boundArgs[name] <- list(NULL)
      } else {
        boundArgs[[name]] <- val
      }
    }
  }
  
  formals(fun)[names(boundArgs)] <- boundArgs
  
  fun
}

clearSeuratCommands <- function(seuratObj, maxSize = 500000) {
  for (commandName in names(seuratObj@commands)) {
    val <- object.size(x = slot(seuratObj@commands[[commandName]], 'call.string'))
    if (val > maxSize) {
      print(paste0('Clearing call.string for: ', commandName, '. size: ', format(val, units = 'auto')))
      slot(seuratObj@commands[[commandName]], 'call.string') <- ''
    }
  }
  
  return(seuratObj)
}

printName <- function(datasetId) {
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  print(paste0('Processing dataset: ', datasetName))
}

savedFiles <- data.frame(datasetId = character(), datasetName = character(), filename = character(), outputFileId = character(), readsetId = character())
if (file.exists('/work/savedSeuratObjects.txt')) {
  print('Deleting pre-existing savedSeuratObjects.txt file')
  unlink('/work/savedSeuratObjects.txt')
}

file.create('/work/savedSeuratObjects.txt')
## [1] TRUE
print(paste0('Total lines in savedSeuratObjects.txt on job start:', length(readLines('savedSeuratObjects.txt'))))
## [1] "Total lines in savedSeuratObjects.txt on job start:0"
saveData <- function(seuratObj, datasetId) {
  print(paste0('Saving dataset: ', datasetId))
  print(seuratObj)
  
  datasetIdForFile <- makeLegalFileName(datasetId)
  fn <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.rds')
  
  message(paste0('Saving RDS file: ', fn, ' with ', ncol(seuratObj), ' cells'))
  barcodeFile <- paste0(outputPrefix, '.', datasetIdForFile, '.cellBarcodes.csv')
  metaFile <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.meta.txt')
  
  saveRDS(seuratObj, file = fn)
  
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  
  # NOTE: this is the ID of the original loupe file. Needed for operations like appending cell hashing or CITE-seq
  outputFileId <- ifelse(datasetId %in% names(datasetIdToOutputFileId), yes = datasetIdToOutputFileId[[datasetId]], no = NA)
  
  readsetId <- ifelse(datasetId %in% names(datasetIdToReadset), yes = datasetIdToReadset[[datasetId]], no = NA)
  print(paste0('readsetId: ', readsetId))
  
  toAppend <- data.frame(datasetId = datasetId, datasetName = datasetName, filename = fn, outputFileId = outputFileId, readsetId = readsetId)
  if (nrow(toAppend) != 1) {
    warning(paste0('Error saving seurat objects, more than one row:'))
    print(toAppend)
    stop('Error saving seurat objects, more than one row!')
  }
  
  write.table(toAppend, file = 'savedSeuratObjects.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)
  print(paste0('Total lines in savedSeuratObjects.txt after save:', length(readLines('savedSeuratObjects.txt'))))
  
# Write cell barcodes and metadata:
metaDf <- seuratObj@meta.data
metaDf$cellbarcode <- colnames(seuratObj)
write.table(metaDf, file = metaFile, quote = T, row.names = F, sep = ',', col.names = T)
write.table(data.frame(CellBarcode = colnames(seuratObj)), file = barcodeFile, quote = F, row.names = F, sep = ',', col.names = F)
}

intermediateFiles <- c()
addIntermediateFile <- function(f) {
  print(paste0('Adding intermediate file: ', f))
  intermediateFiles <<- c(intermediateFiles, f)
}

makeLegalFileName <- function(fn) {
  fn <- gsub(fn, pattern = '\\\\', replacement = '_')
  fn <- gsub(fn, pattern = '[/ ,;]', replacement = '_')
  fn <- gsub(fn, pattern = '\\|', replacement = '_')
  return(fn)
}

errorMessages <- c()
addErrorMessage <- function(f) {
  print(paste0('Adding error: ', f))
  errorMessages <<- c(errorMessages, f)
}

print('Updating future.globals.maxSize')
## [1] "Updating future.globals.maxSize"
options(future.globals.maxSize = Inf)

options('Seurat.memsafe' = TRUE)

if (Sys.getenv('SEURAT_MAX_THREADS') != '') {
  print(paste0('Setting future::plan workers to: ', Sys.getenv('SEURAT_MAX_THREADS')))
  future::plan(strategy='multisession', workers=Sys.getenv('SEURAT_MAX_THREADS'))
}


molInfoFiles <- list(
  '521574-RNA' = '521574.RNA.molInfo.h5',
  '521570-RNA' = '521570.RNA.molInfo.h5',
  '521568-RNA' = '521568.RNA.molInfo.h5',
  '521583-RNA' = '521583.RNA.molInfo.h5',
  '521550-RNA' = '521550.RNA.molInfo.h5',
  '521548-RNA' = '521548.RNA.molInfo.h5',
  '521579-RNA' = '521579.RNA.molInfo.h5',
  '521545-RNA' = '521545.RNA.molInfo.h5',
  '521576-RNA' = '521576.RNA.molInfo.h5',
  '521591-RNA' = '521591.RNA.molInfo.h5',
  '521588-RNA' = '521588.RNA.molInfo.h5',
  '521556-RNA' = '521556.RNA.molInfo.h5',
  '521586-RNA' = '521586.RNA.molInfo.h5',
  '521554-RNA' = '521554.RNA.molInfo.h5',
  '521535-RNA' = '521535.RNA.molInfo.h5',
  '521595-RNA' = '521595.RNA.molInfo.h5',
  '521563-RNA' = '521563.RNA.molInfo.h5',
  '521560-RNA' = '521560.RNA.molInfo.h5'
)
Append Saturation
for (datasetId in names(seuratObjects)) {
  printName(datasetId)
  seuratObj <- readRDS(seuratObjects[[datasetId]])
  
  if (!'DatasetId' %in% names(seuratObj@meta.data)) {
    stop('Seurat object lacks a DatasetId field!')
  }
  
  seuratObj <- CellMembrane::AppendPerCellSaturationInBulk(seuratObj, molInfoList = molInfoFiles)
  saveData(seuratObj, datasetId)
  
  # Cleanup
  rm(seuratObj)
  gc()
}

# Example output info
## [1] "Processing dataset: ONPRC025_XX_D52_B_S2"
## [1] "Adding saturation"
## [1] "Skipping dataset: 521574-RNA"
## [1] "Skipping dataset: 521570-RNA"
## [1] "Skipping dataset: 521568-RNA"
## [1] "Skipping dataset: 521583-RNA"
## [1] "Skipping dataset: 521550-RNA"
## [1] "Skipping dataset: 521548-RNA"
## [1] "Skipping dataset: 521579-RNA"
## [1] "Skipping dataset: 521545-RNA"
## [1] "Skipping dataset: 521576-RNA"
## [1] "Skipping dataset: 521591-RNA"
## [1] "Skipping dataset: 521588-RNA"
## [1] "Skipping dataset: 521556-RNA"
## [1] "Skipping dataset: 521586-RNA"
## [1] "Skipping dataset: 521554-RNA"
## [1] "Skipping dataset: 521535-RNA"
## [1] "Calculating saturation: 521595-RNA"
## [1] "Skipping dataset: 521563-RNA"
## [1] "Skipping dataset: 521560-RNA"
## [1] "Adding saturation to assay RNA for 10815 cells"

## [1] "Saving dataset: 521595"
## An object of class Seurat 
## 35432 features across 10815 samples within 1 assay 
## Active assay: RNA (35432 features, 0 variable features)
## [1] "readsetId: 92312"
## [1] "Total lines in savedSeuratObjects.txt after save:1"

## [1] "Total lines in savedSeuratObjects.txt after save:18"
if (length(intermediateFiles) > 0) {
  write.table(data.frame(file = intermediateFiles), file = 'intermediateFiles.txt', quote = FALSE, delim = '\t', row.names = FALSE, col.names = FALSE)
}

if (length(errorMessages) > 0) {
  print('There were errors:')
  for (msg in errorMessages) {
    print(msg)
  }
  
  write(errorMessages, file = 'seuratErrors.txt')
}

if (file.exists('savedSeuratObjects.txt')) {
  print(paste0('Total lines in savedSeuratObjects.txt:', length(readLines('savedSeuratObjects.txt'))))
} else {
  print('File does not exist: savedSeuratObjects.txt')
}
## [1] "Total lines in savedSeuratObjects.txt:18"


-----------------------------
# Merge Seurat objects
# Conditions
projectName <- 'all18ov_v3'
dietSeurat <- FALSE
errorOnBarcodeSuffix <- TRUE
assaysToDrop <- c('RNA.orig')

outputPrefix <- 'SingleCell.df.frc.saturation.merge'

seuratObjects <- list()
datasetIdToName <- list()
datasetIdToOutputFileId<- list()
datasetIdToReadset<- list()
seuratObjects[['521595']] <- 'SingleCell.df.frc.saturation.521595.seurat.rds'
datasetIdToName[['521595']] <- 'ONPRC025_XX_D52_B_S2'
datasetIdToOutputFileId[['521595']] <- 521595
datasetIdToReadset[['521595']] <- 92312
seuratObjects[['521591']] <- 'SingleCell.df.frc.saturation.521591.seurat.rds'
datasetIdToName[['521591']] <- 'ONPRC017_XX_D130_B_S10'
datasetIdToOutputFileId[['521591']] <- 521591
datasetIdToReadset[['521591']] <- 84652
seuratObjects[['521588']] <- 'SingleCell.df.frc.saturation.521588.seurat.rds'
datasetIdToName[['521588']] <- 'ONPRC027_XX_D50_B_S7'
datasetIdToOutputFileId[['521588']] <- 521588
datasetIdToReadset[['521588']] <- 92310
seuratObjects[['521586']] <- 'SingleCell.df.frc.saturation.521586.seurat.rds'
datasetIdToName[['521586']] <- 'ONPRC018_XX_D50_A_S1'
datasetIdToOutputFileId[['521586']] <- 521586
datasetIdToReadset[['521586']] <- 84655
seuratObjects[['521583']] <- 'SingleCell.df.frc.saturation.521583.seurat.rds'
datasetIdToName[['521583']] <- 'ONPRC016_XX_D130_B_S8'
datasetIdToOutputFileId[['521583']] <- 521583
datasetIdToReadset[['521583']] <- 84651
seuratObjects[['521579']] <- 'SingleCell.df.frc.saturation.521579.seurat.rds'
datasetIdToName[['521579']] <- 'ONPRC017_XX_D130_A_S9'
datasetIdToOutputFileId[['521579']] <- 521579
datasetIdToReadset[['521579']] <- 84648
seuratObjects[['521576']] <- 'SingleCell.df.frc.saturation.521576.seurat.rds'
datasetIdToName[['521576']] <- 'ONPRC025_XX_D52_A_S1'
datasetIdToOutputFileId[['521576']] <- 521576
datasetIdToReadset[['521576']] <- 92309
seuratObjects[['521574']] <- 'SingleCell.df.frc.saturation.521574.seurat.rds'
datasetIdToName[['521574']] <- 'ONPRC027_XX_D50_A_S6'
datasetIdToOutputFileId[['521574']] <- 521574
datasetIdToReadset[['521574']] <- 92311
seuratObjects[['521570']] <- 'SingleCell.df.frc.saturation.521570.seurat.rds'
datasetIdToName[['521570']] <- 'ONPRC012_XX_D100_B_S8'
datasetIdToOutputFileId[['521570']] <- 521570
datasetIdToReadset[['521570']] <- 84650
seuratObjects[['521568']] <- 'SingleCell.df.frc.saturation.521568.seurat.rds'
datasetIdToName[['521568']] <- 'ONPRC018_XX_D50_B_S2'
datasetIdToOutputFileId[['521568']] <- 521568
datasetIdToReadset[['521568']] <- 84654
seuratObjects[['521563']] <- 'SingleCell.df.frc.saturation.521563.seurat.rds'
datasetIdToName[['521563']] <- 'ONPRC012_XX_D100_A_S7'
datasetIdToOutputFileId[['521563']] <- 521563
datasetIdToReadset[['521563']] <- 84653
seuratObjects[['521560']] <- 'SingleCell.df.frc.saturation.521560.seurat.rds'
datasetIdToName[['521560']] <- 'ONPRC016_XX_D130_A_S7'
datasetIdToOutputFileId[['521560']] <- 521560
datasetIdToReadset[['521560']] <- 84646
seuratObjects[['521556']] <- 'SingleCell.df.frc.saturation.521556.seurat.rds'
datasetIdToName[['521556']] <- 'ONPRC001_XX_D130_A_S1'
datasetIdToOutputFileId[['521556']] <- 521556
datasetIdToReadset[['521556']] <- 84649
seuratObjects[['521554']] <- 'SingleCell.df.frc.saturation.521554.seurat.rds'
datasetIdToName[['521554']] <- 'ONPRC001_XX_D130_B_S2'
datasetIdToOutputFileId[['521554']] <- 521554
datasetIdToReadset[['521554']] <- 84647
seuratObjects[['521550']] <- 'SingleCell.df.frc.saturation.521550.seurat.rds'
datasetIdToName[['521550']] <- 'ONPRC006_XX_D100_B_S11'
datasetIdToOutputFileId[['521550']] <- 521550
datasetIdToReadset[['521550']] <- 84645
seuratObjects[['521548']] <- 'SingleCell.df.frc.saturation.521548.seurat.rds'
datasetIdToName[['521548']] <- 'ONPRC007_XX_D100_A_S8'
datasetIdToOutputFileId[['521548']] <- 521548
datasetIdToReadset[['521548']] <- 84643
seuratObjects[['521545']] <- 'SingleCell.df.frc.saturation.521545.seurat.rds'
datasetIdToName[['521545']] <- 'ONPRC006_XX_D100_A_S10'
datasetIdToOutputFileId[['521545']] <- 521545
datasetIdToReadset[['521545']] <- 84644
seuratObjects[['521535']] <- 'SingleCell.df.frc.saturation.521535.seurat.rds'
datasetIdToName[['521535']] <- 'ONPRC007_XX_D100_B_S9'
datasetIdToOutputFileId[['521535']] <- 521535
datasetIdToReadset[['521535']] <- 84642

# Binds arguments from the environment to the target function
bindArgs <- function(fun, seuratObj, allowableArgNames = NULL, disallowedArgNames = NULL) {
  boundArgs <- list()
  boundArgs[['seuratObj']] <- seuratObj
  
  for (name in names(formals(fun))) {
    if (!is.null(disallowedArgNames) && (name %in% disallowedArgNames)) {
      next
    }
    else if (name %in% names(boundArgs)) {
      next
    }
    else if (exists(name)) {
      if (!is.null(allowableArgNames) && !(name %in% allowableArgNames)) {
        next
      }
      
      val <- get(name)
      displayVal <- val
      if (all(is.null(val))) {
        displayVal <- 'NULL'
      } else if (all(is.na(val))) {
        displayVal <- 'NA'
      } else if (is.object(val)) {
        displayVal <- '[object]'
      }
      
      if (length(displayVal) > 1) {
        displayVal <- paste0(displayVal, collapse = ',')
      }
      
      print(paste0('Binding argument: ', name, ': ', displayVal))
      if (all(is.null(val))) {
        boundArgs[name] <- list(NULL)
      } else {
        boundArgs[[name]] <- val
      }
    }
  }
  
  formals(fun)[names(boundArgs)] <- boundArgs
  
  fun
}

clearSeuratCommands <- function(seuratObj, maxSize = 500000) {
  for (commandName in names(seuratObj@commands)) {
    val <- object.size(x = slot(seuratObj@commands[[commandName]], 'call.string'))
    if (val > maxSize) {
      print(paste0('Clearing call.string for: ', commandName, '. size: ', format(val, units = 'auto')))
      slot(seuratObj@commands[[commandName]], 'call.string') <- ''
    }
  }
  
  return(seuratObj)
}

printName <- function(datasetId) {
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  print(paste0('Processing dataset: ', datasetName))
}

savedFiles <- data.frame(datasetId = character(), datasetName = character(), filename = character(), outputFileId = character(), readsetId = character())
if (file.exists('/work/savedSeuratObjects.txt')) {
  print('Deleting pre-existing savedSeuratObjects.txt file')
  unlink('/work/savedSeuratObjects.txt')
}

file.create('/work/savedSeuratObjects.txt')
## [1] TRUE
print(paste0('Total lines in savedSeuratObjects.txt on job start:', length(readLines('savedSeuratObjects.txt'))))
## [1] "Total lines in savedSeuratObjects.txt on job start:0"
saveData <- function(seuratObj, datasetId) {
  print(paste0('Saving dataset: ', datasetId))
  print(seuratObj)
  
  datasetIdForFile <- makeLegalFileName(datasetId)
  fn <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.rds')
  
  message(paste0('Saving RDS file: ', fn, ' with ', ncol(seuratObj), ' cells'))
  barcodeFile <- paste0(outputPrefix, '.', datasetIdForFile, '.cellBarcodes.csv')
  metaFile <- paste0(outputPrefix, '.', datasetIdForFile, '.seurat.meta.txt')
  
  saveRDS(seuratObj, file = fn)
  
  datasetName <- ifelse(datasetId %in% names(datasetIdToName), yes = datasetIdToName[[datasetId]], no = datasetId)
  
  # NOTE: this is the ID of the original loupe file. Needed for operations like appending cell hashing or CITE-seq
  outputFileId <- ifelse(datasetId %in% names(datasetIdToOutputFileId), yes = datasetIdToOutputFileId[[datasetId]], no = NA)
  
  readsetId <- ifelse(datasetId %in% names(datasetIdToReadset), yes = datasetIdToReadset[[datasetId]], no = NA)
  print(paste0('readsetId: ', readsetId))
  
  toAppend <- data.frame(datasetId = datasetId, datasetName = datasetName, filename = fn, outputFileId = outputFileId, readsetId = readsetId)
  if (nrow(toAppend) != 1) {
    warning(paste0('Error saving seurat objects, more than one row:'))
    print(toAppend)
    stop('Error saving seurat objects, more than one row!')
  }
  
  write.table(toAppend, file = 'savedSeuratObjects.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)
  print(paste0('Total lines in savedSeuratObjects.txt after save:', length(readLines('savedSeuratObjects.txt'))))
  
# Write cell barcodes and metadata:
metaDf <- seuratObj@meta.data
metaDf$cellbarcode <- colnames(seuratObj)
write.table(metaDf, file = metaFile, quote = T, row.names = F, sep = ',', col.names = T)
write.table(data.frame(CellBarcode = colnames(seuratObj)), file = barcodeFile, quote = F, row.names = F, sep = ',', col.names = F)
}

intermediateFiles <- c()
addIntermediateFile <- function(f) {
  print(paste0('Adding intermediate file: ', f))
  intermediateFiles <<- c(intermediateFiles, f)
}

makeLegalFileName <- function(fn) {
  fn <- gsub(fn, pattern = '\\\\', replacement = '_')
  fn <- gsub(fn, pattern = '[/ ,;]', replacement = '_')
  fn <- gsub(fn, pattern = '\\|', replacement = '_')
  return(fn)
}

errorMessages <- c()
addErrorMessage <- function(f) {
  print(paste0('Adding error: ', f))
  errorMessages <<- c(errorMessages, f)
}

print('Updating future.globals.maxSize')
## [1] "Updating future.globals.maxSize"
options(future.globals.maxSize = Inf)

options('Seurat.memsafe' = TRUE)

if (Sys.getenv('SEURAT_MAX_THREADS') != '') {
  print(paste0('Setting future::plan workers to: ', Sys.getenv('SEURAT_MAX_THREADS')))
  future::plan(strategy='multisession', workers=Sys.getenv('SEURAT_MAX_THREADS'))
}

# Merge Seurat objects
Merge Seurat Objects
doDiet <- exists('doDiet') && doDiet

mergeBatch <- function(dat) {
  toMerge <- list()
  for (datasetId in names(dat)) {
    message(paste0('Loading: ', datasetId))
    if (doDiet) {
      toMerge[[datasetId]] <- Seurat::DietSeurat(readRDS(dat[[datasetId]]))
      gc()
    } else {
      toMerge[[datasetId]] <- readRDS(dat[[datasetId]])
    }
  }
  
  if (!is.null(assaysToDrop)) {
    for (assayName in assaysToDrop) {
      print(paste0('Dropping assay: ', assayName))
      for (datasetId in names(toMerge)) {
        if (assayName %in% names(toMerge[[datasetId]]@assays)) {
          toMerge[[datasetId]]@assays[[assayName]] <- NULL
        }
      }
    }
  }
  
  seuratObj <- CellMembrane::MergeSeuratObjs(toMerge, projectName = projectName, doGC = doDiet, errorOnBarcodeSuffix = errorOnBarcodeSuffix)
  return(seuratObj)
}

if (length(seuratObjects) == 1) {
  print('There is only one seurat object, no need to merge')
  datasetId <- names(seuratObjects)[[1]]
  saveData(seuratObjects[[datasetId]], datasetId)
} else {
  batchSize <- 20
  numBatches <- ceiling(length(seuratObjects) / batchSize)
  mergedObjects <- list()
  for (i in 1:numBatches) {
    message(paste0('Merging batch ', i, ' of ', numBatches))
    start <- 1 + (i-1)*batchSize
    end <- min(start+batchSize-1, length(seuratObjects))
    message(paste0('processing: ', start, ' to ', end, ' of ', length(seuratObjects)))
    
    mergedObjects[[i]] <- mergeBatch(seuratObjects[start:end])
    gc()
  }
  
  message('Done with batches')
  if (length(mergedObjects) == 1) {
    seuratObj <- mergedObjects[[1]]
  } else {
    message('performing final merge')
    seuratObj <- merge(x = mergedObjects[[1]], y = mergedObjects[2:length(mergedObjects)], project = mergedObjects[[1]]@project.name)
  }
  
  rm(mergedObjects)
  gc()
  
  saveData(seuratObj, projectName)
}

# Example output info 
## [1] "Dropping assay: RNA.orig"
## [1] "Adding dataset: 521595"
## [1] "Adding dataset: 521591"
## [1] "Adding dataset: 521588"
## [1] "Adding dataset: 521586"
## [1] "Adding dataset: 521583"
## [1] "Adding dataset: 521579"
## [1] "Adding dataset: 521576"
## [1] "Adding dataset: 521574"
## [1] "Adding dataset: 521570"
## [1] "Adding dataset: 521568"
## [1] "Adding dataset: 521563"
## [1] "Adding dataset: 521560"
## [1] "Adding dataset: 521556"
## [1] "Adding dataset: 521554"
## [1] "Adding dataset: 521550"
## [1] "Adding dataset: 521548"
## [1] "Adding dataset: 521545"
## [1] "Adding dataset: 521535"
## [1] "Saving dataset: all18ov_v3"
## An object of class Seurat 
## 35432 features across 135335 samples within 1 assay 
## Active assay: RNA (35432 features, 0 variable features)
## [1] "readsetId: NA"
## [1] "Total lines in savedSeuratObjects.txt after save:1"
if (length(intermediateFiles) > 0) {
  write.table(data.frame(file = intermediateFiles), file = 'intermediateFiles.txt', quote = FALSE, delim = '\t', row.names = FALSE, col.names = FALSE)
}

if (length(errorMessages) > 0) {
  print('There were errors:')
  for (msg in errorMessages) {
    print(msg)
  }
  
  write(errorMessages, file = 'seuratErrors.txt')
}

if (file.exists('savedSeuratObjects.txt')) {
  print(paste0('Total lines in savedSeuratObjects.txt:', length(readLines('savedSeuratObjects.txt'))))
} else {
  print('File does not exist: savedSeuratObjects.txt')
}
## [1] "Total lines in savedSeuratObjects.txt:1"


-----------------------------------------------
Session Info
sessionInfo()
## R version 4.2.1 (2022-06-23)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.6 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
##  [1] bookdown_0.33    digest_0.6.31    R6_2.5.1         jsonlite_1.8.4  
##  [5] evaluate_0.20    rmdformats_1.0.4 cachem_1.0.7     rlang_1.1.0     
##  [9] cli_3.6.1        jquerylib_0.1.4  bslib_0.4.2      rmarkdown_2.20  
## [13] tools_4.2.1      xfun_0.37        yaml_2.3.7       fastmap_1.1.1   
## [17] compiler_4.2.1   htmltools_0.5.5  knitr_1.42       sass_0.4.5
