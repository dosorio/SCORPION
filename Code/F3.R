library(ggplot2)
library(uwot)
library(RSpectra)
library(Seurat)
library(ggforce)
library(fgsea)
KEGG <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human')

fileList <- list.files('../Results/Networks/', full.names = TRUE)
inDegrees <- lapply(fileList, function(X){
  load(X)
  colSums(N)
})
geneList <- unique(unlist(lapply(inDegrees, names)))
inDegrees <- sapply(inDegrees, function(X){X[geneList]})
rownames(inDegrees) <- geneList
colnames(inDegrees) <- gsub('.RData', '', basename(fileList))
inDegrees[is.na(inDegrees)] <- 0

load('../Results/ALL.RData')
sideList <- read.csv('../Data/sideList.csv', row.names = 1)
ALL$Side <- sideList[ALL$Donor,1]
levels(ALL$CT)[3] <- 'B'
ctList <- sort(unique(ALL$CT))

DF <- data.frame(ALL@reductions$umap@cell.embeddings, CT = ALL$CT, R = ALL$Region)
DF$R <- as.factor(DF$R)
levels(DF$R) <- c('Border', 'Metastasic', 'Normal', 'Core')
DF$CT <- factor(DF$CT, levels = ctList)
PA <- ggplot(DF[DF$R == 'Normal',], aes(UMAP_1, UMAP_2, color = CT)) + 
  geom_point(cex = 0.01) + theme_light() +
  labs(color = 'Cell Type', shape = 'Region', subtitle = parse(text = paste0('italic(n) == ', nrow(DF[DF$R == 'Normal',])))) +
  theme(legend.position = 'None') +
  xlab('UMAP 1') + ylab('UMAP 2') +
  scale_shape_manual(values= 2) +
  labs(title = 'Normal') +
  theme(plot.title = element_text(face = 2)) + 
  scale_color_manual(values = c('#1d3557', 
                                '#e63946', 
                                '#fb8500', 
                                '#008000', 
                                '#a7c957', 
                                '#0077b6', 
                                '#80ced7', 
                                '#9e0059', 
                                '#9d4edd', 
                                '#b75d69',
                                '#480ca8'))


PB <- ggplot(DF[DF$R == 'Border',], aes(UMAP_1, UMAP_2, color = CT)) + 
  geom_point(cex = 0.01) + theme_light() +
  labs(color = 'Cell Type', shape = 'Region') +
  labs(color = 'Cell Type', shape = 'Region', subtitle = parse(text = paste0('italic(n) == ', nrow(DF[DF$R == 'Border',])))) +
  theme(legend.position = 'None') +
  xlab('UMAP 1') + ylab('UMAP 2') +
  scale_shape_manual(values= 1) +
  labs(title = 'Border') +
  theme(plot.title = element_text(face = 2)) +
  scale_color_manual(values = c('#1d3557', 
                                '#e63946', 
                                '#fb8500', 
                                '#008000', 
                                '#a7c957', 
                                '#0077b6', 
                                '#80ced7', 
                                '#9e0059', 
                                '#9d4edd', 
                                '#b75d69',
                                '#480ca8'))


PC <- ggplot(DF[DF$R == 'Core',], aes(UMAP_1, UMAP_2, color = CT)) + 
  geom_point(cex = 0.01) + theme_light() +
  labs(color = 'Cell Type', shape = 'Region', subtitle = parse(text = paste0('italic(n) == ', nrow(DF[DF$R == 'Core',])))) +
  theme(legend.position = 'None') +
  xlab('UMAP 1') + ylab('UMAP 2') + 
  scale_shape_manual(values= 3) +
  labs(title = 'Core') +
  theme(plot.title = element_text(face = 2)) +
  scale_color_manual(values = c('#1d3557', 
                                '#e63946', 
                                '#fb8500', 
                                '#008000', 
                                '#a7c957', 
                                '#0077b6', 
                                '#80ced7', 
                                '#9e0059', 
                                '#9d4edd', 
                                '#b75d69',
                                '#480ca8'))

PE <- ggplot(DF[DF$R == 'Metastasic',], aes(UMAP_1, UMAP_2, color = CT)) + 
  geom_point(cex = 0.01) + theme_light() +
  labs(color = 'Cell Type', shape = 'Region', subtitle = parse(text = paste0('italic(n) == ', nrow(DF[DF$R == 'Metastasic',])))) +
  theme(legend.position = 'None') +
  xlab('UMAP 1') + ylab('UMAP 2') + 
  scale_shape_manual(values= 3) +
  labs(title = 'Metastasic') +
  theme(plot.title = element_text(face = 2)) +
  scale_color_manual(values = c('#1d3557', 
                                '#e63946', 
                                '#fb8500', 
                                '#008000', 
                                '#a7c957', 
                                '#0077b6', 
                                '#80ced7', 
                                '#9e0059', 
                                '#9d4edd', 
                                '#b75d69',
                                '#480ca8'))


# PA + PB + PC + PE
O <- do.call(rbind.data.frame, strsplit(colnames(inDegrees), '-'))
PCA <- irlba::prcomp_irlba(t(inDegrees), 50)$x
set.seed(1)
DF <- data.frame(Rtsne::Rtsne(PCA[,1:50], check_duplicates = FALSE)$Y, O[,3], O[,2])
#DF <- data.frame(uwot::umap(PCA[,1:30]), O[,3], O[,2])

colnames(DF)[3:4] <- c('CT', 'T')
DF$T <- as.factor(DF$T)
levels(DF$T) <- c('Border', 'Metastasic', 'Normal', 'Core')
DF$T <- factor(DF$T, levels = c('Normal', 'Border', 'Core', 'Metastasic'))
DF$CT <- factor(DF$CT, levels = ctList)
DF$E <- 0
DF$E[((DF$X1 < 0) & (DF$X1 > -10) & (DF$X2 < 0))] <- 1
PD <- ggplot(DF, aes(X1, X2, color = CT, shape = T)) +
  geom_point() + theme_light() +
  xlab('t-SNE 1') + 
  ylab('t-SNE 2') +
  labs(color = 'Cell Type', shape = 'Region', subtitle = parse(text = paste0('italic(n) == ', nrow(DF)))) +
  guides(color=guide_legend(ncol=1), shape = guide_legend(ncol = 1)) +
  scale_shape_manual(values=c(1,2,3,4)) +
  labs(title = 'SCORPION Gene Regulatory Networks') +
  theme(plot.title = element_text(face = 2), legend.title = element_text(face = 2)) +
  scale_color_manual(values = c('#1d3557', 
                                '#e63946', 
                                '#fb8500', 
                                '#008000', 
                                '#a7c957', 
                                '#0077b6', 
                                '#80ced7', 
                                '#9e0059', 
                                '#9d4edd', 
                                '#b75d69',
                                '#480ca8'))

library(patchwork)
pDesign <- 
'
ABEEE
CDEEE
'

png('../Figures/F3.png', width = 3000, height = 1500, res = 300)
PA + PB + PC + PE + PD + plot_layout(design = pDesign) + plot_annotation(tag_levels = 'A')
dev.off()

