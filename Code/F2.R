library(Matrix)
library(Seurat)
library(SCORPION)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(harmony)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
library(Nebulosa)
cellTypes <- fgsea::gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=PanglaoDB_Augmented_2021')
Hallmarks <- fgsea::gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')
KEGG <- fgsea::gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
GOCC <- fgsea::gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2021')

DKO <- Read10X_h5('GSM3477500_Hnf4alphagammaDKO_ScRNAseq_filtered_gene_bc_matrices.h5')
WT <- Read10X_h5('GSM3477499_WT_ScRNAseq_filtered_gene_bc_matrices.h5')

DKO <- DKO[rowMeans(DKO != 0) > 0.05,]
WT <- WT[rowMeans(WT != 0) > 0.05,]

rownames(DKO) <- make.unique(rownames(DKO))
rownames(WT) <- make.unique(rownames(WT))

MOTIF <- read.table('mm10_TF.tsv', header = TRUE)
PPI <- read.table('mm10_PPI.tsv', header = TRUE)

# wtN <- scorpion(motif = MOTIF, expr = WT, PPI)
# save(wtN, file = 'hnf4ag_WT.RData')
# dkoN <- scorpion(motif = MOTIF, expr = DKO, PPI)
# save(dkoN, file = 'hnf4ag_DKO.RData')
load('hnf4ag_WT.RData')
load('hnf4ag_DKO.RData')

gList <- intersect(colnames(wtN$regNet), colnames(dkoN$regNet))

# t.test(wtN$regNet['Hnf4a',gList], dkoN$regNet['Hnf4a',gList])
# t.test(wtN$regNet['Hnf4g',gList], dkoN$regNet['Hnf4g',gList])

# png('HNF4AG.png', width = 2000, height = 1000, res = 300)
# par(mfrow=c(1,2), mar=c(4,4,1,1))
# plot(density(wtN$regNet['Hnf4a',gList]), main = 'Hnf4a')
# lines(density(dkoN$regNet['Hnf4a',gList]), col = 'red')
# 
# plot(density(dkoN$regNet['Hnf4a',gList] - wtN$regNet['Hnf4a',gList]))
# abline(v = 0, col = 'red', lty = 2)
# 
# plot(density(wtN$regNet['Hnf4g',gList]), main = 'Hnf4g')
# lines(density(dkoN$regNet['Hnf4g',gList]), col = 'red')
# dev.off()

D <- colSums(wtN$regNet[,gList]) - colSums(dkoN$regNet[,gList])
names(D) <- toupper(names(D))
D <- D[!grepl('RPL|RPS', names(D))]

X <- data.frame('ko-Hnf4a' = dkoN$regNet['Hnf4a',gList], 
           'wt-Hnf4a' = wtN$regNet['Hnf4a',gList], 
           'ko-Hnf4g' = dkoN$regNet['Hnf4g',gList], 
           'wt-Hnf4g' = wtN$regNet['Hnf4g',gList],
           'd-Hnf4a' = dkoN$regNet['Hnf4a',gList] - wtN$regNet['Hnf4a',gList],
           'd-Hnf4g' = dkoN$regNet['Hnf4g',gList] - wtN$regNet['Hnf4g',gList])

P1 <- reshape2::melt(X[,1:2])
P1$variable <- toupper(gsub('.Hnf4a', '', P1$variable))
P1$variable <- factor(P1$variable, levels = c('WT', 'KO'))
P1 <- ggplot(P1, aes(value, linetype = variable)) +
  geom_density(fill="grey", alpha = 0.2) + 
  theme_light() + 
  xlab(parse(text = 'italic(Hnf4a)~Outdegrees')) +
  ylab('Kernel Density') +
  labs(title = 'Hnf4a', subtitle = parse(text = paste0('italic(n) == ',nrow(X),'~genes')), linetype = '') +
  theme(legend.position = c(0.2, 0.95),legend.background = element_blank(), plot.title = element_text(face = 4), legend.direction = "vertical",
        legend.text = element_text(size = 6), legend.key.size = unit(3, 'mm'), legend.text.align = 0, plot.subtitle = element_text(size = 8)) +
  scale_linetype_manual(values = c(1, 3),
                        labels = c(parse(text = "italic(Hnf4ag)[WT]"), parse(text = "italic(Hnf4ag)[DKO]")))

tP2 <- t.test(X$d.Hnf4a)
P2 <- ggplot(X, aes(d.Hnf4a)) +
  geom_vline(xintercept = 0, lty = 2, color = 'red') + 
  geom_density(fill="grey", alpha = 0.2) + 
  theme_light() + 
  xlab(parse(text = 'KO[italic(Hnf4a)]-WT[italic(Hnf4a)]')) +
  ylab('Kernel Density') +
  labs(subtitle = parse(text = paste0('hat(mu) == ', formatC(tP2$estimate[[1]], digits = 2), '~~~P-value == ', formatC(tP2$p.value, digits = 2)))) +
  theme(plot.subtitle = element_text(size = 8))

P3 <- reshape2::melt(X[,3:4])
P3$variable <- toupper(gsub('.Hnf4g', '', P3$variable))
P3$variable <- factor(P3$variable, levels = c('WT', 'KO'))
P3 <- ggplot(P3, aes(value, linetype = variable)) +
  geom_density(fill="grey", alpha = 0.2) + 
  theme_light() + 
  xlab(parse(text = 'italic(Hnf4g)~Outdegrees')) +
  ylab('Kernel Density') +
  labs(title = 'Hnf4g', subtitle = parse(text = paste0('italic(n) == ',nrow(X),'~genes')), linetype='') +
  theme(legend.position = c(0.2, 0.95),legend.background = element_blank(), plot.title = element_text(face = 4), legend.direction = "vertical",
        legend.text = element_text(size = 6), legend.key.size = unit(3, 'mm'), legend.text.align = 0, plot.subtitle = element_text(size = 8)) + 
  scale_linetype_manual(values = c(1, 3),
                        labels = c(parse(text = "italic(Hnf4ag)[WT]"), parse(text = "italic(Hnf4ag)[DKO]")))

tP4 <- t.test(X$d.Hnf4g)
P4 <- ggplot(X, aes(d.Hnf4g)) +
  geom_vline(xintercept = 0, lty = 2, color = 'red') + 
  geom_density(fill="grey", alpha = 0.2) + 
  theme_light() + 
  xlab(parse(text = 'KO[italic(Hnf4g)]-WT[italic(Hnf4g)]')) +
  ylab('') + 
  labs(subtitle = parse(text = paste0('hat(mu) == ', formatC(tP4$estimate[[1]], digits = 2), '~~~P-value == ', formatC(tP4$p.value, digits = 2)))) +
  theme(plot.subtitle = element_text(size = 8))


# png('dko_Hnf4ag.png', width = 1500, height = 1250, res = 300)
# (P1 | P2)/(P3| P4)
# dev.off()

D <- data.frame(WT = wtN$regNet['Hnf4a', gList], DKO = dkoN$regNet['Hnf4a', gList])
O <- predict(lm(DKO~WT, D), D[,1,drop=FALSE], interval = 'prediction')
O <- O[order(O[,1]),]
O <- as.data.frame(O)
D <- D[rownames(O),]
A <- data.frame(D,O)
A$color <- 'black'
A$color[(A$DKO < A$lwr)] <- 'blue'
A$color[(A$DKO > A$upr)] <- 'red'
A$label <- rownames(A)
A$label[A$color %in% c('black', 'red')] <- NA
A$label[!toupper(A$label) %in% cellTypes$Enterocytes] <- NA
A$nudge <- 0
A$nudge[A$color == 'blue'] <- -0.05

alphaDown <- rownames(A)[A$color == 'blue']
alphaUp <- rownames(A)[A$color == 'red']

tP5 <- statsExpressions::corr_test(A, WT, DKO, type = 'nonp')
P5 <- ggplot(A, aes(WT, DKO, label = label)) +
  geom_line(aes(fit, lwr), lty = 2, col = 'gray', cex = 0.1) +
  geom_line(aes(fit, upr), lty = 2, col = 'gray', cex = 0.1) +
  geom_abline(lty = 1, col = 'gray', cex = 0.5) +
  geom_point(alpha = 0.25, pch = 16, cex = 0.6, color = A$color) +
  geom_density_2d(cex = 0.25, color = 'gray') +
  theme_light() +
  geom_text_repel(color= A$color, 
                  fontface=3, min.segment.length = 0, 
                  nudge_y = A$nudge, 
                  size = 2.5, 
                  segment.size = 0.1, max.overlaps = 75,
                  bg.color = 'white') +
  xlab(parse(text = 'italic(Hnf4ag)[WT]')) +
  ylab(parse(text = 'italic(Hnf4ag)[DKO]')) +
  labs(subtitle = parse(text = paste0('hat(rho) == ', formatC(tP5$estimate, digits = 2)))) +
  theme(plot.subtitle = element_text(size = 8))

W <- scale(D[,2] - D[,1])[,1]
names(W) <- toupper(rownames(D))
set.seed(1)
E <- fgsea::fgseaMultilevel(cellTypes, W, eps = 0)
E$leadingEdge <- unlist(lapply(E$leadingEdge,function(X){paste0(X, collapse = ';')}))
E <- E[order(E$padj),]
write.csv(E, 'Hnf4a_EnrichmentCellTypes.csv', row.names = FALSE)
EP <- formatC(E$padj[E$pathway == 'Enterocytes'], digits = 2)
EN <- formatC(E$NES[E$pathway == 'Enterocytes'], digits = 2)
E1 <- fgsea::plotEnrichment(cellTypes$Enterocytes, W) + 
  ylab('Enrichment\nScore') + 
  xlab(parse(text = 'Outdegree~Difference~Rank')) + 
  labs(title = 'Enterocytes Markers', subtitle = parse(text = paste0('NES: ', EN, '~~~P-adj: ', EP))) +
  theme_light() +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 8))
set.seed(2)
E <- fgsea::fgseaMultilevel(Hallmarks, W, eps = 0)
E$leadingEdge <- unlist(lapply(E$leadingEdge,function(X){paste0(X, collapse = ';')}))
E <- E[order(E$padj),]
write.csv(E, 'Hnf4a_EnrichmentPathways.csv', row.names = FALSE)
EP <- formatC(E$padj[E$pathway == 'Fatty Acid Metabolism'], digits = 2)
EN <- formatC(E$NES[E$pathway == 'Fatty Acid Metabolism'], digits = 2)
E2 <- fgsea::plotEnrichment(Hallmarks$`Fatty Acid Metabolism`, W) + 
  ylab('Enrichment\nScore') + 
  xlab('Gene Rank') + 
  labs(title = 'Fatty Acid\nMetabolism', subtitle = parse(text = paste0('NES: ', EN, '~~~P-adj: ', EP))) +
  theme_light() +
  theme(plot.title = element_text(face = 2, size = 10), plot.subtitle = element_text(size = 8))


D <- data.frame(WT = wtN$regNet['Hnf4g', gList], DKO = dkoN$regNet['Hnf4g', gList])
O <- predict(lm(DKO~WT, D), D[,1,drop=FALSE], interval = 'prediction')
O <- O[order(O[,1]),]
O <- as.data.frame(O)
D <- D[rownames(O),]
A <- data.frame(D,O)
A$color <- 'black'
A$color[(A$DKO < A$lwr)] <- 'blue'
A$color[(A$DKO > A$upr)] <- 'red'
A$label <- rownames(A)
A$label[A$color %in% c('black', 'red')] <- NA
A$label[!toupper(A$label) %in% cellTypes$Enterocytes] <- NA
A$nudge <- 0
A$nudge[A$color == 'blue'] <- -0.05
gammaDown <- rownames(A)[A$color == 'blue']
gammaUp <- rownames(A)[A$color == 'red']

intersect(alphaDown, gammaDown)
intersect(alphaUp, gammaUp)

tP6 <- cor(A$WT, A$DKO, method = 'sp')
P6 <- ggplot(A, aes(WT, DKO, label = label)) +
  geom_line(aes(fit, lwr), lty = 2, col = 'gray', cex = 0.1) +
  geom_line(aes(fit, upr), lty = 2, col = 'gray', cex = 0.1) +
  geom_abline(lty = 1, col = 'gray', cex = 0.5) +
  geom_point(alpha = 0.25, pch = 16, cex = 0.6, color = A$color) +
  geom_density_2d(cex = 0.25, color = 'gray') +
  theme_light() +
  geom_text_repel(color= A$color, 
                  fontface=3, min.segment.length = 0, 
                  nudge_y = A$nudge, 
                  size = 2.5, 
                  segment.size = 0.1, max.overlaps = 75,
                  bg.color = 'white') +
  xlab(parse(text = 'italic(Hnf4ag)[WT]')) +
  ylab(parse(text = 'italic(Hnf4ag)[DKO]')) +
  labs(subtitle = parse(text = paste0('hat(rho) == ', formatC(tP6, digits = 2)))) +
  theme(plot.subtitle = element_text(size = 8))

W <- scale(D[,2] - D[,1])[,1]
names(W) <- toupper(rownames(D))
set.seed(1)
E <- fgsea::fgseaMultilevel(cellTypes, W, eps = 0)
E$leadingEdge <- unlist(lapply(E$leadingEdge,function(X){paste0(X, collapse = ';')}))
E <- E[order(E$padj),]
write.csv(E, 'Hnf4g_EnrichmentCellTypes.csv', row.names = FALSE)
EP <- formatC(E$padj[E$pathway == 'Enterocytes'], digits = 2)
EN <- formatC(E$NES[E$pathway == 'Enterocytes'], digits = 2)
E3 <- fgsea::plotEnrichment(cellTypes$Enterocytes, W) + 
  ylab('Enrichment\nScore') + 
  xlab(parse(text = 'Outdegree~Difference~Rank')) + 
  labs(title = 'Enterocytes Markers', subtitle = parse(text = paste0('NES: ', EN, '~~~P-adj: ', EP))) +
  theme_light() +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 8))
set.seed(2)
E <- fgsea::fgseaMultilevel(Hallmarks, W, eps = 0)
E$leadingEdge <- unlist(lapply(E$leadingEdge,function(X){paste0(X, collapse = ';')}))
E <- E[order(E$padj),]
write.csv(E, 'Hnf4g_EnrichmentPathways.csv', row.names = FALSE)
EP <- formatC(E$padj[E$pathway == 'Fatty Acid Metabolism'], digits = 2)
EN <- formatC(E$NES[E$pathway == 'Fatty Acid Metabolism'], digits = 2)
E4 <- fgsea::plotEnrichment(Hallmarks$`Fatty Acid Metabolism`, W) + 
  ylab('Enrichment\nScore') + 
  xlab('Gene Rank') + 
  labs(title = 'Fatty Acid\nMetabolism', subtitle = parse(text = paste0('NES: ', EN, '~~~P-adj: ', EP))) +
  theme_light() +
  theme(plot.title = element_text(face = 2, size = 10), plot.subtitle = element_text(size = 8))

P2 <- P2 + ylab('Kernel Density')
P4 <- P4 + ylab('Kernel Density')

# pLayout <- '
# AAIBBC
# AAIBBD
# EEJFFG
# EEJFFH
# '
# png('dko_Hnf4ag.png', width = 3200, height = 1900, res = 300)
# P1 + P5 + E1 + E2 + P3 + P6 + E3 + E4 + P2 + P4 + patchwork::plot_layout(design = pLayout)
# #(P1 | P5 | (E1/E2))/(P3 | P6 | (E3/E4))
# dev.off()


#### DUX4 
WT <- readMM('~/MAMMOTH/DUX4/GSM5694433_HNES1_wildtype_matrix.mtx.gz')
rownames(WT) <- read.table('~/MAMMOTH/DUX4/GSM5694433_HNES1_wildtype_features.tsv.gz', sep = '\t')[,2]
colnames(WT) <- readLines('~/MAMMOTH/DUX4/GSM5694433_HNES1_wildtype_barcodes.tsv.gz')[seq_len(ncol(WT))]
WT <- WT[!grepl(':', rownames(WT)),]


OE <- readMM('~/MAMMOTH/DUX4/GSM5694434_HNES1_DUX4_overexpression_matrix.mtx.gz')
rownames(OE) <- read.table('~/MAMMOTH/DUX4/GSM5694434_HNES1_DUX4_overexpression_features.tsv.gz', sep = '\t')[,2]
colnames(OE) <- readLines('~/MAMMOTH/DUX4/GSM5694434_HNES1_DUX4_overexpression_barcodes.tsv.gz')
OE <- OE[!grepl(':', rownames(OE)),]

WT <- CreateSeuratObject(WT, project = 'WT')
OE <- CreateSeuratObject(OE, project = 'OE')
ALL <- merge(WT, OE)
ALL <- scQC(ALL, mtThreshold = 1)
ALL <- NormalizeData(ALL)
gc()
ALL <- FindVariableFeatures(ALL)
ALL <- ScaleData(ALL)
ALL <- RunPCA(ALL)
ALL <- RunHarmony(ALL, group.by.vars = 'orig.ident')
set.seed(1)
ALL <- RunUMAP(ALL, dims = 1:5, reduction = 'harmony')
ALL <- FindNeighbors(ALL, reduction = 'harmony', dims = 1:5)
ALL <- FindClusters(ALL)
table(ALL$seurat_clusters, ALL$orig.ident)
UMAPPlot(ALL, label=TRUE)

# C9 <- ALL@reductions$umap@cell.embeddings
# C9 <- C9[ALL$seurat_clusters == 9,]
# C9 <- C9[C9[,1] > 5.2,]
# C9 <- C9[C9[,2] > 1,]
# C9 <- C9[C9[,2] < 1.5,]
# C9 <- as.data.frame(C9)
UMAPP <- plot_density(ALL, c('ZSCAN4', 'DUXA', 'CCNA1', 'KDM4E'), joint = TRUE, size = 0.4)
UMAPP <- UMAPP[[5]]
UMAPP <- UMAPP +
  theme_light() +
  labs(title = '8C-Like Cells', subtitle = 'ZSCAN4+, DUXA+, CCNA1+, and KDM4E+') +
  theme(legend.position = 'none', plot.title = element_text(face = 2), plot.subtitle = element_text(size = 8, face = 3)) +
  # geom_mark_circle(data = C9, mapping = aes(UMAP_1, UMAP_2), inherit.aes = FALSE, color = 'red', lty = 2) +
  xlim(c(-6.5, 6.5))

UMAPP

markers8CL <- rownames(FindMarkers(ALL, ident.1 = 10, logfc.threshold = 1))


WT <- ALL@assays$RNA@counts[,ALL$orig.ident == 'WT']
OE <- ALL@assays$RNA@counts[,(ALL$orig.ident == 'OE') & (ALL$seurat_clusters == 9)]

WT <- WT[rowMeans(WT != 0) > 0.05,]
OE <- OE[rowMeans(OE != 0) > 0.05,]
WT <- WT[sort(intersect(rownames(WT), rownames(OE))),]
OE <- OE[sort(intersect(rownames(WT), rownames(OE))),]

MOTIF <- read.table('~/hg38_TF.tsv', header = TRUE, sep = '\t')
PPI <- read.table('~/hg38_PPI.tsv', header = TRUE, sep = '\t')

#wtN <- scorpion(motif = MOTIF, expr = WT, PPI)
#save(wtN, file = 'dux4_WT.RData')
#oeN <- scorpion(motif = MOTIF, expr = OE, PPI)
#save(oeN, file = 'dux4_OE.RData')
load('dux4_WT.RData')
load('dux4_OE.RData')

gList <- intersect(colnames(wtN$regNet), colnames(oeN$regNet))

dTF <- rowSums(oeN$regNet) - rowSums(wtN$regNet)
dTF <- cbind(dTF)
plot(scale(dTF))

X <- data.frame(WT=wtN$regNet['DUX4',][gList], 
                OE = oeN$regNet['DUX4',][gList], 
                D = oeN$regNet['DUX4',][gList] - wtN$regNet['DUX4',][gList])
rownames(X) <- gList

P7 <- reshape2::melt(X[,1:2])
P7$variable <- factor(P7$variable, levels = c('WT', 'OE'))
P7 <- ggplot(P7, aes(value, linetype = variable)) +
  geom_density(fill="grey", alpha = 0.2) + 
  theme_light() + 
  xlab(parse(text = 'italic(DUX4)~Outdegrees')) +
  ylab('Kernel Density') +
  labs(title = 'DUX4', subtitle = parse(text = paste0('italic(n) == ',nrow(X),'~genes')), linetype = '') +
  theme(legend.position = c(0.8, 0.95),legend.background = element_blank(), plot.title = element_text(face = 4), legend.direction = "vertical",
        legend.text = element_text(size = 6), legend.key.size = unit(3, 'mm'), legend.text.align = 0, plot.subtitle = element_text(size = 8)) +
  scale_linetype_manual(values = c(1, 3),
                        labels = c(parse(text = "italic(DUX4)[WT]"), parse(text = "italic(DUX4)[OE]")))

tP8 <- t.test(X$D)
P8 <- ggplot(X, aes(D)) +
  geom_vline(xintercept = 0, lty = 2, color = 'red') + 
  geom_density(fill="grey", alpha = 0.2) + 
  theme_light() + 
  xlab(parse(text = 'OE[italic(DUX4)]-WT[italic(DUX4)]')) +
  ylab('Kernel Density') +
  labs(subtitle = parse(text = paste0('hat(mu) == ', formatC(tP8$estimate[[1]], digits = 2), '~~~P-value == ', formatC(tP8$p.value, digits = 2)))) +
  theme(plot.subtitle = element_text(size = 8))

D <- data.frame(WT = wtN$regNet['DUX4',][gList], OE = oeN$regNet['DUX4',][gList])
rownames(D) <- gList
O <- predict(lm(OE~WT, D), D[,1,drop=FALSE], interval = 'prediction')
O <- O[order(O[,1]),]
O <- as.data.frame(O)
D <- D[rownames(O),]
A <- data.frame(D,O)
A$color <- 'black'
A$color[(A$OE < A$lwr)] <- 'blue'
A$color[(A$OE > A$upr)] <- 'red'
A$label <- rownames(A)
A$label[A$color == 'black'] <- NA
A$label[!toupper(A$label) %in% markers8CL] <- NA
A$nudge <- 0
A$nudge[A$color == 'red'] <- 0.5
A$nudge[A$color == 'blue'] <- -0.5

tP9 <- cor(A$WT, A$OE)
P9 <- ggplot(A, aes(WT, OE, label = label)) +
  geom_line(aes(fit, lwr), lty = 2, col = 'gray', cex = 0.1) +
  geom_line(aes(fit, upr), lty = 2, col = 'gray', cex = 0.1) +
  #geom_abline(lty = 1, col = 'gray', cex = 0.5) +
  stat_smooth(method = 'lm', se = FALSE, lty = 1, col = 'gray', cex = 0.5) + 
  geom_point(alpha = 0.25, pch = 16, cex = 0.5, color = A$color) +
  geom_density_2d(cex = 0.25, color = 'gray') +
  theme_light() +
  geom_text_repel(color= A$color, 
                  fontface=3, min.segment.length = 0, 
                  nudge_y = A$nudge, 
                  size = 2.5, 
                  segment.size = 0.1, 
                  bg.color = 'white') +
  xlab(parse(text = 'italic(DUX4)[WT]')) +
  ylab(parse(text = 'italic(DUX4)[OE]')) +
  labs(subtitle = parse(text = paste0('hat(rho) == ', formatC(tP9, digits = 2)))) +
  theme(plot.subtitle = element_text(size = 8))

W <- D$OE-D$WT
names(W) <- rownames(D)
W <- W[complete.cases(W)]
set.seed(3)
E <- fgsea::fgseaMultilevel(list(ECL = markers8CL), W, eps = 0)

EP <- formatC(E$padj[1], digits = 2)
EN <- formatC(E$NES[1], digits = 2)
E5 <- fgsea::plotEnrichment(markers8CL, W) + 
  ylab('Enrichment\nScore') + 
  labs(title = '8C-Like Markers', subtitle = parse(text = paste0('NES: ', EN, '~~~P-adj: ', EP))) +
  theme_light() +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 8)) +
  xlab(parse(text = 'Outdegree~Difference~Rank')) 

dux4Up <- rownames(A)[A$color == 'red']
dux4Down <- rownames(A)[A$color == 'blue']

pLayout <- '
AABCD
AABCD
FFGHI
FFGHI
KKLLN
KKLLN
KKMMO
'
png('F2.png', width = 3600, height = (0.7)*3600, res = 300)
P1 + P2 + P5 + E1 + P3 + P4 + P6 + E3 + UMAPP + P7 + P8 + P9 + E5 + E5 + plot_layout(design = pLayout) + plot_annotation(tag_levels = 'A')
dev.off()

