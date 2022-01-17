data_dir = "/Users/rfigueiredo/contnext_data/data/"

# path
setwd(paste0(data_dir,"data_for_coexp_network_construction"))

library(WGCNA)

options(stringsAsFactors = FALSE)

group = "tissue"
#group = "cell_type"
#group = "cell_line"


#copied from group definition script
# organism part
groupIDs = c('0000178', '0002113', '0002371', '0002097', '0001155', '0001379', 
             '0002048', '0000310', '0004802', '0002107', '0001013', '0000029',
             '0000955', '0002367', '0001134', '0000992', '0002046', '0001295',
             '0001225', '0001264', '0012168', '0001296', '0012652', '0002018',
             '0001052', '0000945', '0003729', '0001836', '0000947', '0002331',
             '0009835', '0002037', '0001507', '0004911', '0001876', '0000002',
             '0001235', '0000317', '0001987', '0016529', '0001377', '0001158',
             '0002038', '0001891', '0002316', '0000173')

# cell type
# groupIDs = c('2000001', '0000842', '0000236', '0000738', '0000576', '0000624',
#              '0000235', '0002328', '0000084', '0000542', '0000451', '0000583',
#              '0000182', '0000057', '0000094', '0002322', '0002327', '0000775',
#              '0000066', '0000037', '0000540', '0000192', '0000625', '0000623',
#              '0002540', '0002620', '0000034', '0000127', '0000134', '0000115')


# cell line
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
  
  datadir = paste0("../coexp_networks/", group, "/", ID, "/")
  
  if(file.exists(paste0(datadir, "coexp_network_edges.tsv"))) {
    next
  }
  
  data = get(load(paste0(datadir, "data_annotated.RData"))) # or however the files will be named
  powers = c(seq(from = 2, to=12, by=2))                                              
  sft = pickSoftThreshold(t(data), powerVector = powers, verbose = 5)                 #
  power = sft$powerEstimate                                                           #
  if (is.na(power)) {                                                                 #
    if ((max(sft$fitIndices[,2]))>0.7) {                                              #
      power = sft$fitIndices[which(sft$fitIndices[,2]==max(sft$fitIndices[,2])),1]    #
    } else {
      if (ID=="9538") {
        power = 30
      } else if (ID=="9952") {
        power = 30
      } else if (ID=="1909") {
        power = 30
      } else if (ID=="0050745") {
        power = 30
      } else {
        powers = c(seq(from = 12, to=30, by=2))                                              
        sft = pickSoftThreshold(t(data), powerVector = powers, verbose = 5)                 #
        power = sft$powerEstimate
      }
    }
  }
  print(power)
  if (!file.exists(paste0(datadir, "coexp_network.rda"))) {
    net = blockwiseModules(t(data), power = power, maxBlockSize = nrow(data),
                           minModuleSize = 20, verbose = 4, #numericLabels = TRUE, 
                           saveTOMs = TRUE, loadTOM = FALSE,
                           saveTOMFileBase = paste0(datadir, "coexp_network"))
    save(net, file = paste0(datadir, "coexp_network.rda"))
  } else {
    load(file=paste0(datadir, "coexp_network.rda"))
  }
  
  # TOM <- TOMsimilarityFromExpr(datExpr, power = power)
  load(paste0(datadir, "coexp_network-block.1.RData"))
  #add names to TOM
  
  vis <- exportNetworkToVisANT(TOM,
                               #file = paste0(datadir, "VisANTInput-", ID),
                               weighted = TRUE,
                               threshold = 0,
                               probeToGene = data.frame(c(seq(from = 1, to=length(rownames(data)))), rownames(data)) )

  threshold = quantile(vis$weight,0.99)
  print(threshold)
  vis2 = subset(vis, weight>=threshold)
  
  write.table(net$colors, file=paste0(datadir, "coexp_network_modules.tsv"), sep="\t", quote=F, row.names=T, col.names=F)
  write.table(vis2, file=paste0(datadir, "coexp_network_edges.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  
  rm(list=setdiff(ls(), c("groupIDs", "ID", "skippedIDs", "group")))
}
