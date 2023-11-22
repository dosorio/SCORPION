library(Matrix)
library(Seurat)
library(SCORPION)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(harmony)
library(ppcor)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
library(Nebulosa)
cellTypes <- fgsea::gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=PanglaoDB_Augmented_2021')
Hallmarks <- fgsea::gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')
KEGG <- fgsea::gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
GOCC <- fgsea::gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2021')

DKO <- Read10X_h5('GSM3477500_Hnf4agDKO.h5')
WT <- Read10X_h5('GSM3477499_Hnf4agWT.h5')

DKO <- DKO[rowMeans(DKO != 0) > 0.05,]
WT <- WT[rowMeans(WT != 0) > 0.05,]

rownames(DKO) <- make.unique(rownames(DKO))
rownames(WT) <- make.unique(rownames(WT))

WT <- SCORPION:::makeSuperCells(WT[rowSums(WT) > 0,], gamma = 10, n.pc = 25, fast.pca = FALSE)
DKO <- SCORPION:::makeSuperCells(DKO[rowSums(DKO) > 0,], gamma = 10, n.pc = 25, fast.pca = FALSE)

gList <- intersect(rownames(WT), rownames(DKO))

ppcor_WT <- ppcor::pcor(as.matrix(t(WT)), method = 'pearson')
ppcor_DKO <- ppcor::pcor(as.matrix(t(DKO)), method = 'pearson')


colnames(ppcor_WT$estimate) <- rownames(ppcor_WT$estimate) <- rownames(WT)
colnames(ppcor_DKO$estimate) <- rownames(ppcor_DKO$estimate) <- rownames(DKO)

plot(ppcor_WT$estimate['Hnf4a', gList], ppcor_DKO$estimate['Hnf4a', gList])

cor(ppcor_WT$estimate['Hnf4a', gList], ppcor_DKO$estimate['Hnf4a', gList])
plot(density(ppcor_DKO$estimate['Hnf4a', gList] - ppcor_WT$estimate['Hnf4a', gList]))


X <- data.frame('ko-Hnf4a' = ppcor_DKO$estimate['Hnf4a',gList], 
                'wt-Hnf4a' = ppcor_WT$estimate['Hnf4a',gList], 
                'ko-Hnf4g' = ppcor_DKO$estimate['Hnf4g',gList], 
                'wt-Hnf4g' = ppcor_WT$estimate['Hnf4g',gList],
                'd-Hnf4a' = ppcor_DKO$estimate['Hnf4a',gList] - ppcor_WT$estimate['Hnf4a',gList],
                'd-Hnf4g' = ppcor_DKO$estimate['Hnf4g',gList] - ppcor_WT$estimate['Hnf4g',gList])

P1 <- reshape2::melt(X[,1:2])
P1$variable <- toupper(gsub('.Hnf4a', '', P1$variable))
P1$variable <- factor(P1$variable, levels = c('WT', 'KO'))
P1 <- ggplot(P1, aes(value, linetype = variable)) +
  geom_density(fill="grey", alpha = 0.2) + 
  theme_light() + 
  xlab(parse(text = 'italic(Hnf4a)~Edges~Weight')) +
  ylab('Kernel Density') +
  labs(title = 'Hnf4a', subtitle = parse(text = paste0('italic(n) == ',nrow(X),'~genes')), linetype = '') +
  theme(legend.position = c(0.1, 0.95),legend.background = element_blank(), plot.title = element_text(face = 4), legend.direction = "vertical",
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
  xlab(parse(text = 'italic(Hnf4g)~Edges~Weight')) +
  ylab('Kernel Density') +
  labs(title = 'Hnf4g', subtitle = parse(text = paste0('italic(n) == ',nrow(X),'~genes')), linetype='') +
  theme(legend.position = c(0.1, 0.95),legend.background = element_blank(), plot.title = element_text(face = 4), legend.direction = "vertical",
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

D <- data.frame(WT = ppcor_WT$estimate['Hnf4a', gList], DKO = ppcor_DKO$estimate['Hnf4a', gList])
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
EP <- formatC(E$padj[E$pathway == 'Enterocytes'], digits = 2)
EN <- formatC(E$NES[E$pathway == 'Enterocytes'], digits = 2)
E1 <- fgsea::plotEnrichment(cellTypes$Enterocytes, W) + 
  ylab('Enrichment\nScore') + 
  xlab(parse(text = 'Edges~Weight~Difference~Rank')) + 
  labs(title = 'Enterocytes Markers', subtitle = parse(text = paste0('NES: ', EN, '~~~P-adj: ', EP))) +
  theme_light() +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 8))
set.seed(2)
E <- fgsea::fgseaMultilevel(Hallmarks, W, eps = 0)
E$leadingEdge <- unlist(lapply(E$leadingEdge,function(X){paste0(X, collapse = ';')}))
E <- E[order(E$padj),]
EP <- formatC(E$padj[E$pathway == 'Fatty Acid Metabolism'], digits = 2)
EN <- formatC(E$NES[E$pathway == 'Fatty Acid Metabolism'], digits = 2)
E2 <- fgsea::plotEnrichment(Hallmarks$`Fatty Acid Metabolism`, W) + 
  ylab('Enrichment\nScore') + 
  xlab('Gene Rank') + 
  labs(title = 'Fatty Acid\nMetabolism', subtitle = parse(text = paste0('NES: ', EN, '~~~P-adj: ', EP))) +
  theme_light() +
  theme(plot.title = element_text(face = 2, size = 10), plot.subtitle = element_text(size = 8))


D <- data.frame(WT = ppcor_WT$estimate['Hnf4g', gList], DKO = ppcor_DKO$estimate['Hnf4g', gList])
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
EP <- formatC(E$padj[E$pathway == 'Enterocytes'], digits = 2)
EN <- formatC(E$NES[E$pathway == 'Enterocytes'], digits = 2)
E3 <- fgsea::plotEnrichment(cellTypes$Enterocytes, W) + 
  ylab('Enrichment\nScore') + 
  xlab(parse(text = 'Edges~Weight~Difference~Rank')) + 
  labs(title = 'Enterocytes Markers', subtitle = parse(text = paste0('NES: ', EN, '~~~P-adj: ', EP))) +
  theme_light() +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 8))
set.seed(2)
E <- fgsea::fgseaMultilevel(Hallmarks, W, eps = 0)
E$leadingEdge <- unlist(lapply(E$leadingEdge,function(X){paste0(X, collapse = ';')}))
E <- E[order(E$padj),]
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

pLayout <- '
AAIBC
AAIBC
EEJFG
EEJFG
'
png('ppcor_dko_Hnf4ag.png', width = 3500, height = 1500, res = 300)
P1 + P5 + E1 + P3 + P6 + E3 + P2 + P4 + patchwork::plot_layout(design = pLayout)
#(P1 | P5 | (E1/E2))/(P3 | P6 | (E3/E4))
dev.off()