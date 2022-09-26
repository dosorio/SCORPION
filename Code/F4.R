library(ggplot2)
library(patchwork)
library(pbapply)
library(Rfast)
library(dplyr)
library(fgsea)
library(ggrepel)
library(Matrix)
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
fileList <- fileList[grepl('Epithelial', fileList)]
groupFactorALL <- factor(unlist(lapply(strsplit(basename(fileList), '-'), function(X){X[2]})), levels = c('N', 'B', 'T', 'M'))
fileList <- fileList[!is.na(groupFactorALL)]
groupFactorALL <- factor(unlist(lapply(strsplit(basename(fileList), '-'), function(X){X[2]})), levels = c('N', 'B', 'T', 'M'))
fileList <- fileList[grepl('Epithelial', fileList)]
fileContent <- lapply(fileList, loadNetwork)
gList <- unique(unlist(lapply(fileContent, colnames)))
tfList <- unique(unlist(lapply(fileContent, rownames)))
fileContent <- pblapply(fileContent, function(X){makeComparable(X = X, tfList = tfList, gList = gList)})
fileContent <- sapply(fileContent, function(X){reshape2::melt(as.matrix(X))[,3]})                           
edgeList <- expand.grid(tfList, gList)

qLM <- function(X, y){
  X <- cbind(1,X)
  n <- nrow(X)
  p <- ncol(X)
  m <- solve(t(X) %*% X) %*% t(X) %*% t(y)
  yHat <- X %*% m
  y <- t(y)
  MSM <- rowSums((t(yHat) - colMeans(y))^2)/(p-1)
  MSE <- colSums((y - yHat)^2)/(n-p)
  P <- pf(MSM/MSE, df1 = p-1, df2 = n-p, lower.tail = FALSE)
  O <- data.frame(beta = t(m)[,2],P)
  return(O)
}

O <- data.frame(edgeList, qLM(as.numeric(groupFactorALL),fileContent))
O$FDR <- p.adjust(O$P, method = 'fdr')
O <- O[O$FDR < 0.05,]
colnames(O) <- c('tf', 'target', 'beta', 'P', 'FDR')
O <- O[order(abs(O$beta), decreasing = TRUE),]
write.csv(O, '../Results/dNet-Epithelial-Betas.csv', row.names = FALSE)


crcBetas <- read.csv('../Results/dNet-Epithelial-Betas.csv')
crcBetas <- crcBetas[crcBetas$target %in% unique(unlist(H[c('Wnt-beta Catenin Signaling', 'Hedgehog Signaling')])),]
crcBetas <- crcBetas[order(abs(crcBetas$beta), decreasing = TRUE),]

A <- rbind(
  crcBetas[crcBetas$tf %in% c('EGR2') & crcBetas$target %in% c('HDAC5'),],
  crcBetas[crcBetas$tf %in% c('HOXD8') & crcBetas$target %in% c('VEGFA'),],
  crcBetas[crcBetas$tf %in% c('TFDP1') & crcBetas$target %in% c('NCSTN'),],
  crcBetas[crcBetas$tf %in% c('NANOG') & crcBetas$target %in% c('ADGRG1'),],
  crcBetas[crcBetas$tf %in% c('PAX6') & crcBetas$target %in% c('HDAC2'),],
  crcBetas[crcBetas$tf %in% c('SP1') & crcBetas$target %in% c('CCND2'),]
)

A <- A[order(A$beta, decreasing = TRUE),]
betasPlot <- apply(A, 1, function(X){
  df <- data.frame(W = fileContent[edgeList$Var1 %in% X[1] &edgeList$Var2 %in% X[2],], G = groupFactorALL)
  levels(df$G) <- c('N', 'B', 'C', 'M')
  
  summaryDF <- df %>% group_by(G) %>% summarise(m = mean(W), lb = m - sd(W), ub = m + sd(W))
  bCoef <- coefficients(lm(df$W~as.numeric(df$G)))[2]
  bColor <- ifelse(bCoef > 0, '#d90429', '#1d3557')
  arrowSymbol <- ifelse(bCoef > 0, '\u279e', '\u21e5')
  ggplot(df, aes(G, W)) + 
    geom_boxplot(outlier.colour = NA, color = 'gray80') + 
    geom_jitter(width = 0.1, alpha = 0.1) + 
    theme_bw() +
    xlab('Tumor Region') +
    ylab('Edge Weight') +
    ggtitle(paste0(X[1],' ', arrowSymbol,' ', X[2])) +
    #labs(title = parse(text = paste0('italic(', X[1], ')', '~symbol("\u21e5")~', 'italic(', X[2], ')'))) +
    labs(subtitle = parse(text = paste0('beta == ', round(bCoef,3)))) +
    theme(panel.border = element_rect(color = bColor, fill = NA, size = 1)) +
    geom_point(data = summaryDF, mapping = aes(G, m), color = bColor, pch = 8, cex = 3) +
    geom_errorbar(data = summaryDF, mapping = aes(x = G, y = m, ymin = lb, ymax = ub), width = 0.1, color = bColor) +
    theme(plot.title = element_text(face = 3))
  })

betasPlot[[1]] <- betasPlot[[1]] + labs(tag = 'A')
# eval(parse(text = paste0(paste0('betasPlot[[', seq_along(betasPlot), ']]'), collapse = ' +')))
# 
# png('../Figures/betasPlot.png', width = 2000, height = 1800, res = 300)
# eval(parse(text = paste0(paste0('betasPlot[[', seq_along(betasPlot), ']]'), collapse = ' +')))
# dev.off()
# "betasPlot[[1]] +betasPlot[[2]] +betasPlot[[3]] +betasPlot[[4]] +betasPlot[[5]] +betasPlot[[6]]"

crcBetas <- read.csv('../Results/dNet-Epithelial-Betas.csv')
outDegree <- crcBetas %>% group_by(tf) %>% summarise(W = sum(beta))
outDegree <- outDegree[order(outDegree$W, decreasing = TRUE),]
outDegree$R <- seq_along(outDegree[[1]])
outDegree <- outDegree[complete.cases(outDegree),]
outDegree$label <- outDegree$tf
outDegree$label[11:(nrow(outDegree)-10)] <- NA
write.csv(outDegree[,1:3], row.names = FALSE, quote = FALSE, file = '../Results/CRC-tfRanks.csv')

PB <- ggplot(outDegree, aes(R, W, color = W, label = label)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_gradient2(low = '#1d3557', high = '#d90429', midpoint = 0) +
  geom_text_repel(min.segment.length = 0, size = 3, 
                  nudge_x = ifelse(outDegree$W > 0, 100,-150), 
                  point.padding = 0, 
                  segment.size = 0.1, 
                  segment.curvature = 0.01, 
                  bg.color = 'white') +
  xlab('Transcription Factor Rank') +
  ylab('Outdegrees') +
  theme(legend.position = 'None') +
  labs(tag = 'B')

inDegree <- crcBetas %>% group_by(target) %>% summarise(W = sum(beta))
I <- inDegree[[2]]
names(I) <- inDegree[[1]]
I <- I[is.finite(I)]

set.seed(1)
E <- fgseaMultilevel(H, I)
E <- E[E$padj < 0.05,]
E <- E[order(E$NES),]

UR <- E[E$NES > 0,]
DR <- E[E$NES < 0,]
UR$pathway <- factor(UR$pathway, levels = UR$pathway)
DR$pathway <- factor(DR$pathway, levels = DR$pathway)

E <- rbind(UR,DR)
E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
write.csv(E, quote = FALSE, row.names = FALSE, file = '../Results/CRC-EnrichmentPathways.csv')

PC <- ggplot(UR, aes(NES, pathway, fill = -log10(pval))) +
  geom_bar(stat = 'identity') +
  scale_fill_gradient(high = '#d90429', low = '#edf2f4', limits = c(0,max(-log10(UR$pval)))) + 
  theme_bw() +
  ylab('Hallmarks') +
  labs(fill = parse(text = '-log[10]~(P-value)')) +
  labs(tag = 'C')

PD <- ggplot(DR, aes(NES, pathway, fill = -log10(pval))) +
  geom_bar(stat = 'identity') +
  scale_fill_gradient(high = '#1d3557', low = '#edf2f4', limits = c(0,max(-log10(DR$pval)))) + 
  theme_bw() +
  ylab('Hallmarks') +
  labs(fill = parse(text = '-log[10]~(P-value)')) +
  labs(tag = 'D')

library(patchwork)
pDesign <- '
AABBCCGGGGG
AABBCCGGGGG
DDEEFFGGGGG
DDEEFFGGGGG
'

png('../Figures/F4A.png', width = 3500, height = 1600, res = 300)
betasPlot[[1]] + betasPlot[[2]] + betasPlot[[3]] + betasPlot[[4]] + betasPlot[[5]] + betasPlot[[6]] + PB + PC + PD + plot_layout(design = pDesign)
dev.off()

png('../Figures/F4B.png', width = 3500, height = 650, res = 300)
PC + PD
dev.off()
