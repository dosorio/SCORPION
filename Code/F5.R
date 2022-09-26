library(ggplot2)
library(patchwork)
library(pbapply)
library(Rfast)
library(dplyr)
library(fgsea)
library(ggpubr)
library(Matrix)
library(SummarizedExperiment)
library(dorothea)
library(igraph)

H <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')

loadNetwork <- function(X){
  load(X)
  return(N)
}

makeComparable <- function(X, tfList, gList){
  X <- as.matrix(X)
  O <- matrix(data = 0, nrow = length(tfList), ncol = length(gList), dimnames = list(tfList, gList))
  O[rownames(X), colnames(X)] <- X
  O <- Matrix(O)
  gc()
  return(O)
}

fileList <- list.files('../Results/Networks/', full.names = TRUE)
fileList <- fileList[grepl('T-Epithelial', fileList)]
fileContent <- lapply(fileList, loadNetwork)
gList <- unique(unlist(lapply(fileContent, colnames)))
tfList <- unique(unlist(lapply(fileContent, rownames)))
fileContent <- pblapply(fileContent, function(X){makeComparable(X = X, tfList = tfList, gList = gList)})
fileContent <- sapply(fileContent, function(X){reshape2::melt(as.matrix(X))[,3]})                           
edgeList <- expand.grid(tfList, gList)

donorList <- unlist(lapply(strsplit(basename(fileList), '-'), function(X){X[1]}))
sideList <- read.csv('../Data/sideList.csv', header = FALSE, row.names = 1)
sideList <- sideList[donorList,]

fileContent <- fileContent[,!is.na(sideList)]
sideList <- sideList[!is.na(sideList)]

O <- ttests(t(fileContent), ina = (sideList == 'L')+1)
O <- data.frame(edgeList, O)
O <- na.omit(O)
O$FDR <- p.adjust(O$pvalue, method = 'fdr')
O <- O[O$FDR < 0.05,]
colnames(O) <- c('tf', 'target', 'stat', 'pvalue', 'dof', 'fdr')

NFKB2 <- O[O$tf == 'NFKB2',]

O <- O[O$tf %in% c('ZNF350', 'CEBPG', 'HOXD8', 'HOXA13', 'FOXD2', 'MEOX2', 'NFYA', 'NANOG', 'SOX21', 'IRF5',
                   'NFKB2', 'EBF1', 'NFKB1', 'ZNF768', 'FLI1', 'KLF9', 'SP1', 'ZNF520', 'KLF1', 'KLF12'),]
O <- as.data.frame.array(O)
O <- O[order(abs(O$stat), decreasing = TRUE),]
O <- lapply(split(O, O$target), function(X){na.omit(X[1:2,])})
O <- do.call(rbind.data.frame, O)
O <- O[order(abs(O$stat), decreasing = TRUE),]
# O <- split(O, O$tf)
# O <- O[unlist(lapply(O, nrow)) > 5]
# O <- do.call(rbind.data.frame, O)

X <- graph_from_data_frame(O)
set.seed(1)
oLayout <- layout_with_dh(X)
W <- O$stat * -log10(O$pvalue)
vColor <- rgb(1, 1, 1, 0.2)
eColor <- ifelse(W > 0, yes = "#d90429", no = "#1d3557")
nC <- nchar(names(degree(X)))
vShape <- 'circle' #
vFont <- ifelse(names(V(X)) %in% tfList, 4, 3)
vSize <- (4 + (nC/max(nC) * 8))
vSize[vFont %in% 3] <- nchar(names(V(X))[vFont %in% 3])

png('F5.png', width = 2000, height = 2000, res = 300)
par(mar=c(0,0,0,0))
plot(
  X,
  layout = oLayout,
  vertex.size = vSize,
  edge.arrow.size = 0.15,
  edge.width = (abs(W)/max(abs(W))) * 1, 
  edge.color = eColor,
  edge.curved = .2,
  vertex.label.cex = 0.45,
  vertex.shape = vShape,
  vertex.color = vColor,
  vertex.frame.color = vColor,
  vertex.label.family = "Arial",
  vertex.label.color = 'black',
  vertex.label.font = vFont
)
dev.off()
