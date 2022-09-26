library(Seurat)
library(SCORPION)
library(SuperCell)
library(Matrix)
library(pbapply)

TF <- read.table('~/SCORPION/hg38_TF.tsv', sep = '\t', header = TRUE)
PPI <- read.table('~/SCORPION/hg38_PPI.tsv', sep = '\t', header = TRUE)

load('../Results/PreprocessedData/ALL.RData')
allIdents <- paste0(ALL$orig.ident, '-', ALL$CT)
ALL <- ALL@assays$RNA@counts

pbsapply(names(table(allIdents)[table(allIdents) >= 30]), function(ID){
  outFile <- paste0('../Results/Networks/',ID, '.RData')
  if(!file.exists(outFile)){
    idData <- ALL[,allIdents %in% ID]
    if(ncol(idData) >= 30){
      idData <- idData[rowSums(idData != 0) > 5,]
      N <- scorpion(tfMotifs = TF, gexMatrix = idData, ppiNet = PPI)[[1]]
      N <- round(N, 3)
      gc()
      save(N, file = outFile)
    }
  }
})
