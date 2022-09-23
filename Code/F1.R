library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(Matrix)
library(igraph)

AUPRC <- list.files(pattern = 'AUPRC')
AUPRC <- lapply(AUPRC, function(X){
  D <- read.csv(X)
  colnames(D)[1] <- 'METHOD'
  colnames(D)[2] <- paste0('AUPRC_', colnames(D)[2])
  return(D)
})[[1]]
rownames(AUPRC) <- AUPRC[,1]

AUROC <- list.files(pattern = 'AUROC')
AUROC <- lapply(AUROC, function(X){
  D <- read.csv(X)
  colnames(D)[1] <- 'METHOD'
  colnames(D)[2] <- paste0('AUROC_', colnames(D)[2])
  return(D)
})[[1]]
rownames(AUROC) <- AUROC[,1]

TIME <- list.files(pattern = 'Time')
TIME <- lapply(TIME, function(X){
  D <- read.csv(X)
  colnames(D)[1] <- 'METHOD'
  colnames(D)[2] <- paste0('TIME_', colnames(D)[2])
  D[,2] <- 1/D[,2]
  return(D)
})[[1]]
rownames(TIME) <- TIME[,1]

countMatrix <- read.csv('../../inputs/example/GSD/ExpressionData.csv', row.names = 1)
gMean <- rowMeans(countMatrix)

COR <- list.files('GSD', full.names = TRUE)
COR <- paste0(COR, '/rankedEdges.csv')
COR <- COR[!grepl('.csv/', COR)]
COR <- sapply(COR, function(W){
  X <- read.csv(W, sep = '\t')
  colnames(X) <- c('from','to', 'weight')
  X <- graph_from_data_frame(X)
  
  wDegree <- (rowSums(X[]) + colSums(X[]))/2
  #plot(gMean[names(wDegree)], wDegree)
  cor(gMean[names(wDegree)], wDegree, method = 'sp')
})
names(COR) <- unlist(lapply(strsplit(names(COR), '/'), function(X){X[2]}))

D <- merge(AUPRC,AUROC, by = 'METHOD')
D <- merge(D, TIME, by = 'METHOD')
D$COR <- 1/abs(COR[D$METHOD])
FBL <- read.csv('GSD-NetworkMotifs-FBL.csv', row.names = 1)
FFL <- read.csv('GSD-NetworkMotifs-FFL.csv', row.names = 1)
MI <- read.csv('GSD-NetworkMotifs-MI.csv', row.names = 1)
D$FBL <- FBL[D$METHOD,]
D$FFL <- FFL[D$METHOD,]
D$MI <- MI[D$METHOD,]
rD <- apply(abs(D[,2:ncol(D)]),2,rank)
#D <- D[order(rowMeans(rD), decreasing = FALSE),]
D <- D[c(9, 4, 2, 13, 8, 10, 1, 3, 5, 7, 6, 12, 11),]
write.csv(D, file = 'allBenchmark.csv', row.names = FALSE)

rownames(D) <- D$METHOD
D <- D[,-1]
D <- as.matrix(D)
rank <- reshape2::melt(as.matrix(apply(D, 2, rank)))[,3]
D[,3] <- round(1/D[,3],1)
D[,4] <- abs(COR[rownames(D)])
O <- reshape2::melt(D)
colnames(O) <- c('method', 'feature', 'value')
O$value <- round(O$value,2)
O$rank <- rank
O$rank[is.na(O$value)] <- NA
O$method <- factor(O$method)
levels(O$feature) <- c('AUPRC', 'AUROC', 'TIME[s]', 'BIAS', 'FBL', 'FFL', 'MI')
PA <- ggplot(O, aes(feature, method)) + 
  geom_tile(aes(fill = rank), color = 'white', lwd = 1) + 
  geom_text(aes(label = value), color = "white", size = 4) +
  theme_bw() +
  theme(legend.position = 'None') +
  ylab('') +
  xlab('') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(face = c(1,1,1,1,1,1,1,1,1,1,1,1,2))) +
  scale_x_discrete(labels = c('AUPRC', 'AUROC', parse(text = 'TIME[s]'), 'BIAS', parse(text = 'MOTIF[FBL]'), parse(text = 'MOTIF[FFL]'), parse(text = 'MOTIF[MI]'))) + 
  scale_fill_gradientn(colours = c('#1d3557', '#d90429')) + labs(tag = 'B')
  
PB <- ggplot() + theme_void() + labs(tag = 'A')
pDesign <- 'AAABB'
png('bm.png', width = 3300, height = 1650, res = 300)  
PB + PA + plot_layout(design = pDesign)
dev.off()
