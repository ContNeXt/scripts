## need to install these packages in advance: 
#BiocManager::install("hgu133a2cdf")

setwd("/Users/rfigueiredo/Documents/contNeXt/R")

library("GEOquery")
library("affy")
library("oligo")

# to get the publivation information
#citation("GEOquery")

### dir definitions
data_dir = "/Users/rfigueiredo/contnext_data/data/"
rawData_path = paste0(data_dir, "raw_data/") 
exprsData_path = paste0(data_dir, "exprs_data/") 

## load all data
#metadataFinal = read.csv(paste0(data_dir, "metadata/final_metadata.tsv"), sep="\t")

datasets = as.character(unique(metadataFinal$dataset))
#datasets = datasets[1:100]

datasetsCausingProblems = readLines("datasetsCausingProblems.txt") #create an empty file, add to it as you run into problems
datasetsNoData = readLines("datasetsNoData.txt") #create an empty file, add to it as you run into problems
pb <- txtProgressBar(min = 0, max = length(datasets)-1, style = 3)

suppressMessages(
  for (i in 1:length(datasets)) {
    if (datasets[i]%in%datasetsCausingProblems) {
      next
    }
    fname_exprsSet = paste0(exprsData_path, datasets[i], ".Rdata")
    if (file.exists(fname_exprsSet)) {
      next
    }
    raw_data_file = getGEOSuppFiles(datasets[i], filter_regex = "*_RAW.tar", fetch= FALSE, baseDir = rawData_path)
    getGEOSuppFiles(datasets[i], filter_regex = "*_RAW.tar", baseDir = rawData_path)
    
    if (!is.null(raw_data_file)) {
      if (!(length(raw_data_file$fname)==0)) {
        if (length(raw_data_file$fname)>1) {
          print("Problem1")
        } else {
          untar(paste0(rawData_path, datasets[i], "/", raw_data_file$fname), exdir = paste0(rawData_path, datasets[i]))
          cel_files = list.files(paste0(rawData_path, datasets[i]), pattern = "*.CEL|*.cel")
          rawData <- ReadAffy(filenames=paste0(rawData_path, datasets[i], "/", cel_files))
          normed_data = affy::rma(rawData)
          save(normed_data, file = fname_exprsSet)
        }
      } else {
        #print(datasets[i])
        datasetsNoData = c(datasetsNoData, datasets[i])
        #print(datasetsNoData)
      }
    } else {
      datasetsNoData = c(datasetsNoData, datasets[i])
    }
    setTxtProgressBar(pb, i-1)
    rm(list=setdiff(ls(), c("datasets", "datasetsCausingProblems", "datasetsNoData", "rawData_path", "exprsData_path", "pb", "i")))
    filesToRemove = paste0(tempdir(), "/", list.files(tempdir()))
    filesToRemove = setdiff(filesToRemove, paste0(tempdir(), "/", "A-AFFY-44.adf.txt"))
    #filesToRemove = setdiff(filesToRemove, '/tmp/RtmpVpoO8X/rs-graphics-1b408c82-c603-41ef-b5a3-26723f71f4eb')
    file.remove(filesToRemove)
    gc()
  }
)
close(pb)

write(datasetsNoData, ncolumns = 1, file="datasetsNoData.txt")

datasets[i]
i

#####------------- remove metadata of experiments that couldn't be loaded -----------------------
## load all data
rm(list=setdiff(ls(), "exprsData_path"))
metadataFinal = read.csv(paste0(data_dir, "metadata/final_metadata.tsv"), sep="\t")

datasetsALL = unique(metadataFinal$dataset)
#datasetsLoaded = list.files("Z:/NO BACKUP/exprsSets_rma/")
datasetsLoaded = list.files(exprsData_path)
datasetsLoaded = unlist(strsplit(datasetsLoaded, ".Rdata"))

missingDatasets = setdiff(datasetsALL, datasetsLoaded)
#datasetsCausingProblems = readLines("datasetsCausingProblems.txt")
#datasetsToDelete = readLines("datasetsToDelete.txt")

#missingDatasets = setdiff(missingDatasets, datasetsCausingProblems)
#missingDatasets = setdiff(missingDatasets, datasetsToDelete)
length(missingDatasets)==0

metadataFinal = subset(metadataFinal, dataset%in%datasetsLoaded)
dim(metadataFinal)
#xxxx samples, yyy experiments

save(metadataFinal, file="metadataFinal_afterDataLoading.RData")
write.table(metadataFinal, file="metadataFinal_afterDataLoading.tsv", sep="\t", quote=F, row.names=F, col.names=T)
