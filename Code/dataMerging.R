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
