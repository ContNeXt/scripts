data_dir = "/Users/rfigueiredo/contnext_data/data/"

# path
setwd(paste0(data_dir,"data_for_coexp_network_construction"))

options(stringsAsFactors=FALSE)
rm(list=ls())
library("ArrayExpress")
library("sva")
library("jsonlite")
library("data.table")
library("dplyr")
library("tidyr")
library("genefilter")
library("hgu133plus2.db")

group = "organism.part"
#group = "cell.type"
#group = "cell.line"

#copied from group definition script
organism part
groupIDs = c('0000178', '0002113', '0002371', '0002097', '0001155', '0001379', 
             '0002048', '0000310', '0004802', '0002107', '0001013', '0000029',
             '0000955', '0002367', '0001134', '0000992', '0002046', '0001295',
             '0001225', '0001264', '0012168', '0001296', '0012652', '0002018',
             '0001052', '0000945', '0003729', '0001836', '0000947', '0002331',
             '0009835', '0002037', '0001507', '0004911', '0001876', '0000002',
             '0001235', '0000317', '0001987', '0016529', '0001377', '0001158',
             '0002038', '0001891', '0002316', '0000173')
 
# # cell type
# groupIDs = c('2000001', '0000842', '0000236', '0000738', '0000576', '0000624',
#              '0000235', '0002328', '0000084', '0000542', '0000451', '0000583',
#              '0000182', '0000057', '0000094', '0002322', '0002327', '0000775',
#              '0000066', '0000037', '0000540', '0000192', '0000625', '0000623',
#              '0002540', '0002620', '0000034', '0000127', '0000134', '0000115')
# 
# 
# # cell line
# groupIDs = c('0007606', '0003684', '0001582', '0009989', '0009015', '0001601',
#              '0007634', '0037116', '0004307', '0007986', '0003665', '0001230',
#              '0009357', '0009348', '0037163', '0037295', '0003704', '0001571',
#              '0037117', '0009456', '0002699', '0002172')
 


skippedIDs = c() 

for(ID in groupIDs) {
  
  print(paste0("Next group: ", ID))
  
  if (ID%in%skippedIDs) {
    next
  }
  
  datasets = readLines(paste0(group, "/", ID, "/datasets.txt"))
  
  #load the correct metadata file (the one where you have also defined your per disease datasets from)
  #metadataFinal = read.table(file = paste0(data_dir, "metadata/", "metadataFinal_afterDataLoading.tsv"), sep = '\t', header = TRUE, quote = "")
  load(file=paste0(data_dir, "metadata/", "metadataFinal_afterDataLoading.RData"))
  
  #load the metadata file that has only samples for a specific sub context
  metadata = read.table(file = paste0(group, "/", ID,"/metadata.tsv"), sep = '\t', header = TRUE, quote = "", fill = TRUE)
  
  source("/Users/rfigueiredo/Documents/ContNeXt/R/getCelFile.R")
  
  datadir = paste0("../build_coexp_network/", group, "/", ID, "/")  # this is where I want the files for the next script to go
  dl_datadir = "/Users/rfigueiredo/contnext_data/exprs_data/" #this is where the downloaded data is
  if (!dir.exists(datadir)) {
    dir.create(datadir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  }
  
  if(file.exists(paste0(datadir, "data_annotated.RData"))) {
    next
  }
  
  ###### Download data
  datasetFiles = paste0(dl_datadir, datasets, ".Rdata")
  
  
  if (!(all(file.exists(paste0(datadir, datasets, ".Rdata"))))) {
    print("Start downloading datasets!")

    pb <- txtProgressBar(min = 0, max = length(datasets)-1, style = 3)

    suppressMessages(
      for (i in 1:length(datasets)) {
        #file.exists(paste0(datadir, datasets[7], ".Rdata"))
        fname_exprsSet = paste0(datadir, datasets[i], ".Rdata")
        if (file.exists(fname_exprsSet)) {
          next
        }
        rawset = ArrayExpress(datasets[i])
        if (!is.null(rawset)) {
          if (length(rawset)>1) {
            indx = which(names(rawset)=="A-AFFY-44")
            rawData = rawset[[indx]]
            AEsetnorm = oligo::rma(rawData)
            save(AEsetnorm, file = fname_exprsSet)
          } else {
            AEsetnorm = oligo::rma(rawset)
            save(AEsetnorm, file = fname_exprsSet)
          }
        }
        setTxtProgressBar(pb, i-1)
      }
    )

    close(pb)
    print("Finished downloading datasets!")

  } else {
    print("Download was already done!")
  }
  
  datasets = datasets[file.exists(paste0(dl_datadir, datasets, ".Rdata"))]
  metadataFinal = subset(metadataFinal, dataset%in%datasets)
  metadata = subset(metadata, dataset%in%datasets)
  
  
  ###### Merge data
  if (!file.exists(paste0(datadir, "exprsALL.RData"))) {
    
    getSampleIDs = function() {
      tryCatch({
        sampleIDs = sampleData_dataset$sampleID[sapply(sampleNames, FUN = grep, x = sampleData_dataset$acc_sample)]
        return(sampleIDs)
      }, error = function(msg) {
        sampleData_dataset$celFiles = getCelFileNamesFromAE(datasets[i])
        sampleNames = colnames(exprsData)
        sampleIDs = sampleData_dataset$sampleID[match(sampleNames,sampleData_dataset$celFiles)]
        return(sampleIDs)
      })
    }
    
    print("Start merging datasets!")
    
    pb <- txtProgressBar(min = 0, max = length(datasets)-1, style = 3)
    
    exprsALL = matrix(data = , nrow = 54675, ncol = 0) 
    
    batchALL = c()
    suppressMessages(
      #for (i in 1:length(datasets)) {
      for (i in 700:length(datasets)) {
        fname_exprsSet = paste0(dl_datadir, datasets[i], ".Rdata")
        load(fname_exprsSet) #normed_data
        exprsData = exprs(normed_data)
        if (all(is.na(exprsData))) {
          datasets = datasets[-i]
          next
        }
        sampleNames = colnames(exprsData)
        sampleNames = unlist(strsplit(sampleNames, ".CEL.gz|.cel.gz|.CEL|.cel|_.*"))
        sampleData_dataset = subset(metadataFinal, subset=dataset==datasets[i])
        colnames(exprsData) = sampleNames
        
        if (all(rownames(exprsALL)==rownames(exprsData))) {
          exprsALL = cbind(exprsALL, exprsData)
          batchALL = c(batchALL, rep(datasets[i], ncol(exprsData)))
        } else {
          print("Error: unequal probe identifier")
          break
        }
        setTxtProgressBar(pb, i-1)
      }
    )
    
    save(exprsALL, batchALL, file=paste0(datadir, "exprsALL.RData"))
    print("Finished merging datasets!")
  } else {
    print("Merging was already done!")
    load(file=paste0(datadir, "exprsALL.RData"))
  }
  
  dim(exprsALL)
  
  metadataFinal = subset(metadataFinal, sample_id%in%colnames(exprsALL))
  metadata = subset(metadata, sample_id%in%colnames(exprsALL))
  
  
  ###### Batch correct data
  if (!file.exists(paste0(datadir, "data_batchCorrected.RData"))) {
    print("Start batch correcting data!")
    #metadata_subset = subset(metadataFinal, dataset%in%datasets)
    modcombat = model.matrix(~1, data=as.data.frame(t(exprsALL)))
    #metadata_subset)
    #combat_edata = ComBat(dat=exprsALL, batch=metadata_subset$dataset, mod=modcombat,
    #                      par.prior=TRUE, prior.plots=FALSE)()
    combat_edata = ComBat(dat=exprsALL, batch= batchALL,#as.factor(as.character(metadata_subset$dataset)),
                          mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
    save(combat_edata, file=paste0(datadir, "data_batchCorrected.RData"))
    print("Finished batch correcting data!")
  } else {
    print("Batch correction was already done!")
    load(file=paste0(datadir, "data_batchCorrected.RData"))
  }
  
  
  ###### Filter for needed samples only
  if (!file.exists(paste0(datadir, "data_batchCorrected_filtered.RData"))) {
    print("Start filtering datasets for needed samples!")
    sampleIDs = metadata$sample_id
    edata_filtered = subset(combat_edata, select = colnames(combat_edata)%in%sampleIDs)
    #in the case of duplicate column names:
    #colnames(edata_filtered)[colnames(edata_filtered) == "GSM578906"] = c("GSM578906_1", "GSM578906_2")
    save(edata_filtered, file=paste0(datadir, "data_batchCorrected_filtered.RData"))
    print("Finished filtering datasets for needed samples!")
  } else {
    print("Filtering datasets for needed samples was already done!")
    load(file=paste0(datadir, "data_batchCorrected_filtered.RData"))
  }
  
  
  ###### Probe filtering and annotation
  if (!file.exists(paste0(datadir, "data_annotated.RData"))) {
    print("Start probe filtering and annotation datasets!")
    eset_combat = ExpressionSet(edata_filtered, annotation = "hgu133plus2") 
    eset_combat_filtered = featureFilter(eset_combat, require.entrez=TRUE, remove.dupEntrez=TRUE,
                                         feature.exclude="^AFFX")
    
    mapGeneNames_probe2HGNC = function(probeIDs) {
      genes_HGNC = unlist(mapIds(hgu133plus2.db,
                                 keys=probeIDs,
                                 column="SYMBOL",
                                 keytype="PROBEID",
                                 multiVals="first"))
      return(genes_HGNC)
    }
    
    combat_edata_filtered = exprs(eset_combat_filtered)
    
    symbols = mapGeneNames_probe2HGNC(rownames(combat_edata_filtered))
    combat_edata_annotated = combat_edata_filtered[!is.na(symbols),]
    symbols_filtered = symbols[!is.na(symbols)]
    rownames(combat_edata_annotated) = symbols_filtered
    
    save(combat_edata_annotated, file=paste0(datadir, "data_annotated.RData"))
    print("Finished probe filtering and annotation datasets!")
    print("Everything is prepared for WGCNA.")
  } else {
    print("Probe filtering and annotation was already done!")
    print("Everything is prepared for WGCNA. Nothing to do anymore.")
  }
  
  
  
}
