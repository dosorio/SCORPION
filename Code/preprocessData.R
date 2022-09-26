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
