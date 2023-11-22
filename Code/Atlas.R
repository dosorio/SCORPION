library(Seurat)
library(Matrix)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

# GSE132465
GSE132465 <- read.table('../Data/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz', header = TRUE, sep = '\t', row.names = 1)
gc()
GSE132465 <- as.matrix(GSE132465)
gc()
GSE132465 <- as(GSE132465, 'dgCMatrix')
gc()
newBC <- colnames(GSE132465)
newBC <- do.call(rbind.data.frame, strsplit(newBC, '\\.|\\_'))
colnames(newBC) <- c('Patient', 'Type', 'Barcode')
colnames(GSE132465) <- paste0(newBC[,1], '-', newBC[,2], '_', newBC[,3])
writeMM(GSE132465, '../Results/PreprocessedData/GSE132465_matrix.mtx')
writeLines(rownames(GSE132465), '../Results/PreprocessedData/GSE132465_genes.txt')
writeLines(colnames(GSE132465), '../Results/PreprocessedData/GSE132465_barcodes.txt')

# EMTAB8107
EMTAB8107 <- readMM('../Data/2098_matrix.mtx')
colnames(EMTAB8107) <- readLines('../Data/2098_barcodes.tsv')
rownames(EMTAB8107) <- read.table('../Data/2098_genes.tsv')[,1]
mdEMTAB8107 <- read.csv('../Data/2099-Colorectalcancer_metadata.csv.gz')
newBC <- colnames(EMTAB8107)
newBC <- do.call(rbind.data.frame, strsplit(newBC, '_'))
newBC[,3] <- newBC[,2]
newBC[,2] <- mdEMTAB8107$TumorSite
newBC[,2] <- gsub('C', 'T', newBC[,2])
newBC[,1] <- paste0('P', mdEMTAB8107$PatientNumber)
colnames(newBC) <- c('Patient', 'Type', 'Barcode')
colnames(EMTAB8107) <- paste0(newBC[,1], '-', newBC[,2], '_', newBC[,3])
writeMM(EMTAB8107, '../Results/PreprocessedData/EMTAB8107_matrix.mtx')
writeLines(rownames(EMTAB8107), '../Results/PreprocessedData/EMTAB8107_genes.txt')
writeLines(colnames(EMTAB8107), '../Results/PreprocessedData/EMTAB8107_barcodes.txt')

# GSE144735
GSE144735 <- read.table('../Data/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz', header = TRUE, sep = '\t', row.names = 1)
gc()
GSE144735 <- as.matrix(GSE144735)
gc()
GSE144735 <- as(GSE144735, 'dgCMatrix')
gc()
newBC <- colnames(GSE144735)
newBC <- do.call(rbind.data.frame, strsplit(newBC, '\\.|\\_'))
colnames(newBC) <- c('Patient', 'Type', 'Barcode')
colnames(GSE144735) <- paste0(newBC[,1], '-', newBC[,2], '_', newBC[,3])
writeMM(GSE144735, '../Results/PreprocessedData/GSE144735_matrix.mtx')
writeLines(rownames(GSE144735), '../Results/PreprocessedData/GSE144735_genes.mtx')
writeLines(colnames(GSE144735), '../Results/PreprocessedData/GSE144735_barcodes.mtx')

# PMID35365629
PMID35365629 <- readRDS('../Data/41467_2022_29366_MOESM6_ESM.rds.gz')
newBC <- colnames(PMID35365629)
newBC <- do.call(rbind.data.frame, strsplit(newBC, '_'))
newBC <- newBC[,-1]
colnames(newBC) <- c('Patient', 'Type', 'Barcode')
newBC[,1] <- toupper(newBC[,1])
colnames(PMID35365629) <- paste0(newBC[,1], '-', newBC[,2], '_', newBC[,3])
writeMM(PMID35365629, '../Results/PreprocessedData/PMID35365629_matrix.mtx')
writeLines(rownames(PMID35365629), '../Results/PreprocessedData/PMID35365629_genes.txt')
writeLines(colnames(PMID35365629), '../Results/PreprocessedData/PMID35365629_barcodes.txt')

# GSE178318
GSE178318 <- readMM('../Data/GSE178318_matrix.mtx.gz')
rownames(GSE178318) <- read.table('../Data/GSE178318_genes.tsv.gz', sep = '\t')[,2]
colnames(GSE178318) <- readLines('../Data/GSE178318_barcodes.tsv.gz')
newBC <- colnames(GSE178318)
newBC <- do.call(rbind.data.frame, strsplit(newBC, '_'))
colnames(newBC) <- c('Barcode', 'Patient', 'Type')
newBC[,2] <- toupper(newBC[,2])
GSE178318 <- GSE178318[,newBC$Type %in% c('CRC', 'LM')]
newBC <- newBC[newBC$Type %in% c('CRC', 'LM'),]
newBC$Type <- as.factor(newBC$Type)
levels(newBC$Type) <- c('T', 'M')
colnames(GSE178318) <- paste0(newBC[,2], '-', newBC[,3], '_', newBC[,1])
writeMM(GSE178318, '../Results/PreprocessedData/GSE178318_matrix.mtx')
writeLines(rownames(GSE178318), '../Results/PreprocessedData/GSE178318_genes.txt')
writeLines(colnames(GSE178318), '../Results/PreprocessedData/GSE178318_barcodes.txt')


library(Seurat)
library(harmony)
library(Nebulosa)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')

getSeuratObjects <- function(path, X){
  W <- Matrix::readMM(paste0(path, X, '_matrix.mtx'))
  rownames(W) <- readLines(paste0(path, X, '_genes.txt'))
  colnames(W) <- readLines(paste0(path, X, '_barcodes.txt'))
  W <- Seurat::CreateSeuratObject(W, project = X)
  return(W)
}

EMTAB8107 <- getSeuratObjects(path = '../Results/PreprocessedData/', X = 'EMTAB8107')
GSE132465 <- getSeuratObjects(path = '../Results/PreprocessedData/', X = 'GSE132465')
ALL <- merge(EMTAB8107, GSE132465)
rm(EMTAB8107, GSE132465)
gc()
GSE144735 <- getSeuratObjects(path = '../Results/PreprocessedData/', X = 'GSE144735')
ALL <- merge(ALL, GSE144735)
rm(GSE144735)
gc()
PMID35365629 <- getSeuratObjects(path = '../Results/PreprocessedData/', X = 'PMID35365629')
ALL <- merge(ALL, PMID35365629)
rm(PMID35365629)
gc()
GSE178318 <- getSeuratObjects(path = '../Results/PreprocessedData/', X = 'GSE178318')
ALL <- merge(ALL, GSE178318)
rm(GSE178318)
gc()

ALL$MT <- colSums(ALL@assays$RNA@counts[grepl('^MT-', rownames(ALL@assays$RNA@counts)),])/ALL$nCount_RNA
ALL$Donor <- unlist(lapply(strsplit(colnames(ALL), '-'), function(X){X[1]}))
ALL$Region <- unlist(lapply(strsplit(colnames(ALL), '-|\\_'), function(X){X[2]}))

ALL <- scQC(ALL)
ALL <- NormalizeData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- ScaleData(ALL)
ALL <- RunPCA(ALL)
ALL <- RunHarmony(ALL, 'orig.ident')
ALL <- RunUMAP(ALL, reduction = 'harmony', dims = 1:50)
ALL <- FindNeighbors(ALL, reduction = 'harmony', dims = 1:50)
ALL <- FindClusters(ALL)

save(ALL, file = '../Results/PreprocessedData/ALL.RData')

markerList <- c('EPCAM', 'CDH1', 'COL1A1', 'COL3A1', 'MS4A1', 
                'CDH5', 'PECAM1', 'S100B', 'CDH2', 'PTPRC', 
                'CD3E', 'CD14', 'IL1RL1', 'MZB1', 'MKI67', 'PLD4')

png('atlasMarkers.png', width = 4000*.85, height = 3000*.85, res = 300)
plot_density(ALL, markerList) & 
  theme_light() & 
  theme(plot.title = element_text(face = 4), legend.key.width = unit(0.2, 'cm')) &
  scale_color_gradient(low = '#e5e5e5', high = '#d90429')
dev.off()

newID <- rep(NA, 39)
newID[c(1, 2, 3, 5, 11, 13, 18, 26, 29:39)] <- 'T'
newID[c(4, 24)] <- 'Epithelial'
newID[c(6, 25)] <- 'B'
newID[c(7, 8, 14)] <- 'Myeloid'
newID[c(9, 20, 29)] <- 'Plasma'
newID[c(10, 16, 17, 21, 28)] <- 'MSC'
newID[c(12, 27)] <- 'Endothelial'
newID[c(15)] <- 'Plasmablast'
newID[c(19)] <- 'Mast'
newID[c(22)] <- 'Glial'
newID[c(23)] <- 'Dendritic'


Idents(ALL) <- ALL$seurat_clusters
levels(Idents(ALL)) <- newID
ALL$CT <- Idents(ALL)
UMAPPlot(ALL, label = TRUE)
save(ALL, file = '../Results/PreprocessedData/ALL.RData')


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
