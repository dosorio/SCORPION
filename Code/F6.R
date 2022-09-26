library(ggplot2)
library(patchwork)
library(pbapply)
library(Rfast)
library(dplyr)
library(fgsea)
library(ggpubr)
library(Matrix)
library(SummarizedExperiment)
library(enrichR)

H <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020')
BP <- gmtPathways('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2021')
MMR <- BP$`mismatch repair (GO:0006298)`

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

fileContent <- data.frame(edgeList, fileContent)
R <- fileContent %>% group_by(Var1) %>% summarise(across(colnames(fileContent)[-c(1:2)], sum))
C <- as.matrix(R[,2:ncol(R)])
C <- C[,!is.na(sideList)]
sideList <- sideList[!is.na(sideList)]
rownames(C) <- R$Var1
O <- ttests(t(C), ina = (sideList == 'L')+1)
O <- as.data.frame(O)
O$FDR <- p.adjust(O[,2], method = 'fdr')
rownames(O) <- R$Var1
O$L <- NA
O <- O[order(O$stat),]
O$L[1:10] <- rownames(O)[1:10]
O <- O[order(O$stat, decreasing = TRUE),]
O$L[1:10] <- rownames(O)[1:10]
colFun <- colorRampPalette(c('#d90429', 'gray90', '#1d3557'))
O$stat <- -1 * O$stat
O$color <- colFun(nrow(O))
write.csv(O[,1:4], quote = FALSE, file = 'LVR_TFDifferences.csv')

# rTF <- enrichr(rownames(O)[O$stat < 0 & O$FDR < 0.05], 'MSigDB_Hallmark_2020')
# rTF <- do.call(rbind.data.frame, rTF)
# lTF <- enrichr(rownames(O)[O$stat > 0 & O$FDR < 0.05], 'MSigDB_Hallmark_2020')
# lTF <- do.call(rbind.data.frame, lTF)

png('../Figures/LvR.png', width = 2000, height = 1000, res = 300)
set.seed(1)
LVR <- ggplot(O, aes(stat, -log10(pvalue),label = L)) + 
  geom_point(color = O$color, pch = ifelse(O$FDR < 0.05, 8, 16)) + theme_bw() +
  ggrepel::geom_text_repel(size = 2.5, 
                           min.segment.length = 0, 
                           #nudge_x = ifelse(O$stat > 0, -0.5, 0.5), 
                           fontface = 3, 
                           segment.size = 0.1, 
                           bg.color = 'white', max.overlaps = 100) +
  xlim(c(1.1*min(O$stat),1.1*max(O$stat))) +
  ylim(0,3.6) +
  xlab(parse(text = 'italic(t)')) +
  ylab(parse(text = '-log[10]~(P-value)')) +
  labs(tag = 'A')
print(LVR)
dev.off()

bpList <- lapply(rev(na.omit(O$L)), function(i){
  DF <- data.frame(W = scale(C[i,])[,1], S = sideList)
  DF$S <- factor(DF$S, levels = c('R', 'L'))
  bColor <- O[i,6]
  P <- ggplot(DF, aes(S, W)) + 
    geom_boxplot(outlier.colour = NA) + 
    labs(title = i) + 
    theme_bw() + 
    theme(panel.border = element_rect(color = bColor, fill = NA, size = 1)) +
    theme(plot.title = element_text(face = 4)) +
    geom_jitter(alpha = 0.25, width = 0.1) +
    stat_compare_means(comparisons = list(c(1,2)), method = 't.test', label = "p.signif") +
    ylim(c(-3,4)) +
    ylab(parse(text = 'Z~(Outdegree)')) +
    xlab('Side')
  # if(i %in% c('ZNF350', 'MEOX2', 'KLF12', 'FLI1')){
  #   P <- P + ylab(parse(text = 'Z~(Outdegree)'))
  # }
  # if(i %in% c('NANOG','NFKB1')){
  #   P <- P + xlab(parse(text = 'Side'))
  # }
  return(P)
})

png('../Figures/LvR2.png', width = 2000, height = 2000, res = 300)
eval(parse(text = paste0(paste0('bpList[[',seq_along(bpList),']]'), collapse = ' + ')))
dev.off()

pDesign <- '
AAAAA####
AAAAA####
AAAAA####
AAAAA####
BCDEF####
GHIJK####
LMNOP####
QRSTU####'
bpList[[1]] <- bpList[[1]] + labs(tag = 'B')
bpList[[11]] <- bpList[[11]] + labs(tag = 'C')

png('../Figures/LvR3.png', width = 2000, height = 3000, res = 300)
LVR + bpList[[1]] + bpList[[2]] + bpList[[3]] + bpList[[4]] + bpList[[5]] + 
  bpList[[6]] + bpList[[7]] + bpList[[8]] + bpList[[9]] + bpList[[10]] + 
  bpList[[11]] + bpList[[12]] + bpList[[13]] + bpList[[14]] + bpList[[15]] + 
  bpList[[16]] + bpList[[17]] + bpList[[18]] + bpList[[19]] + bpList[[20]] + plot_layout(design = pDesign)
dev.off()

# TCGA
library(TCGAbiolinks)
library(SummarizedExperiment)
coadQuery <- GDCquery(
  project = c("TCGA-COAD", 'TCGA-READ'),
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)
GDCdownload(coadQuery)
TCGA <- GDCprepare(coadQuery)

sideTCGA <- as.factor(TCGA@colData$tissue_or_organ_of_origin)
levels(sideTCGA) <- c('R', 'R', NA, NA, 'L', 'R', 'L', 'L', 'L', 'L', 'R')
GEX <- assay(TCGA)
GEX <- GEX[,TCGA@colData$definition %in% c('Primary solid Tumor')]
#GEX <- log1p(t(t(GEX)/colSums(GEX))* 1e6)
sideTCGA <- sideTCGA[TCGA@colData$definition %in% c('Primary solid Tumor')]

DF <- data.frame(W = GEX['NFKB2',], S = sideTCGA)
DF$S <- factor(DF$S, levels = c('R', 'L'))
DF <- DF[complete.cases(DF),]
png('../Figures/sideNFKB2.png', width = 750, height = 1500, res = 300)
PD <- ggplot(na.omit(DF), aes(S,W)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.1, alpha = 0.1) +
  theme_bw() +
  ylab(parse(text = 'FPKM')) +
  xlab('Side') +
  stat_compare_means(comparisons = list(c(1,2)), label = 'p.signif') +
  ylim(c(min(DF$W), 1.2*max(DF$W))) +
  labs(title = 'NFKB2') +
  theme(plot.title = element_text(face = 4)) +
  theme(panel.border = element_rect(color = '#d90429', fill = NA, size = 1))
print(PD)
dev.off()

DF <- data.frame(W = GEX['ZNF350',], S = sideTCGA)
DF$S <- factor(DF$S, levels = c('R', 'L'))
DF <- DF[complete.cases(DF),]
PF <- ggplot(na.omit(DF), aes(S,W)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.1, alpha = 0.1) +
  theme_bw() +
  ylab(parse(text = 'FPKM')) +
  xlab('Side') +
  stat_compare_means(comparisons = list(c(1,2)), label = 'p.signif') +
  ylim(c(min(DF$W), 1.25*max(DF$W))) +
  labs(title = 'ZNF350') +
  theme(plot.title = element_text(face = 4)) +
  theme(panel.border = element_rect(color = '#1d3557', fill = NA, size = 1))

A <- as.data.frame(TCGA@colData)
A <- A[A$barcode %in% colnames(GEX),]
library(survival)
plot(survfit(Surv(A$days_to_last_follow_up/365)~GEX['NFKB2',]>=median(GEX['NFKB2',])))
survdiff(Surv(A$days_to_last_follow_up/365)~GEX['NFKB2',]>=median(GEX['NFKB2',]))

SD <- summary(survfit(Surv(A$days_to_last_follow_up/365)~GEX['NFKB2',]>=median(GEX['NFKB2',])))
SD <- data.frame(T = SD$time, S = SD$surv, G = SD$strata)
levels(SD$G) <- c('Low', 'High')

png('../Figures/SPlot.png', width = 1000, height = 1000, res = 300)
PE <- ggplot(SD, aes(T, S, lty = G)) + 
  geom_step() + 
  theme_bw() + 
  labs(lty = parse(text = 'italic(NFKB2)')) +
  ylab('Survival Probability') +
  xlab('Time (Years)') +
  annotate(geom="text", x=7, y=0.85, label="Log Rank Test\nP-value = 0.042\nHR = 1.771") +
  theme(legend.position = c(0.15,0.12), legend.background = element_rect(fill = NA), legend.key.height = unit(.3, units = 'cm'))
print(PE)
dev.off()

pDesign <- '
XXAAAAAAA
XXAAAAAAA
XXAAAAAAA
XXAAAAAAA
BCDEFVVVV
GHIJKVVVV
LMNOPWWWW
QRSTUYYYY'

PA <- ggplot() + theme_minimal() + labs(tag = 'A')
LVR <- LVR + labs(tag = 'B')
bpList[[1]] <- bpList[[1]] + labs(tag = 'C')
bpList[[11]] <- bpList[[11]] + labs(tag = 'D')
PE <- PE + labs(tag = 'D')
PF <- PF + labs(tag = 'E')
PD <- PD + labs(tag = 'F')

library(ggplot2)
library(ggpubr)

fileList <- list.files(pattern = '.tab', path = '~/Desktop/LVR/', full.names = TRUE)

fileContent <- lapply(fileList, function(X){read.table(X, row.names = 1)})
geneList <- rownames(fileContent[[1]])
fileContent <- sapply(fileContent, function(X){X[geneList,1]})
rownames(fileContent) <- geneList
colnames(fileContent) <- gsub('ReadsPerGene.out.tab|_','',basename(fileList))
fileContent <- fileContent[-c(1:4),]
fileContent <- log1p(t(t(fileContent)/colSums(fileContent)) * 1e6)

df <- data.frame(E = fileContent[grepl('ENSG00000256683', rownames(fileContent)),], S = ifelse(grepl('CRC', colnames(fileContent)), 'L', 'R'))
df$S <- factor(df$S, levels = c('R', 'L'))
P1 <- ggplot(df, aes(S,E)) + 
  geom_boxplot(outlier.color = NA) + 
  stat_compare_means(comparisons = list(c(1,2)), 
                     method.args = list(alternative = "less"),
                     label = 'p.signif') +
  geom_jitter(width = 0.1, alpha = 0.1) +
  labs(title = 'ZNF350') +
  ylab(parse(text = 'log(CPM + 1)')) +
  theme_light() +
  ylim(c(min(df$E), 1.2*max(df$E))) + 
  theme(plot.title = element_text(face = 4)) +
  xlab('Side') +
  theme(panel.border = element_rect(color = '#1d3557', fill = NA, size = 1))


df <- data.frame(E = fileContent[grepl('ENSG00000077150', rownames(fileContent)),], S = ifelse(grepl('CRC', colnames(fileContent)), 'L', 'R'))
df$S <- factor(df$S, levels = c('R', 'L'))
P2 <- ggplot(df, aes(S,E)) + 
  geom_boxplot(outlier.color = NA) + 
  stat_compare_means(comparisons = list(c(1,2)), 
                     method.args = list(alternative = "greater"),
                     label = 'p.signif') +
  geom_jitter(width = 0.1, alpha = 0.1) +
  labs(title = 'NFKB2') +
  ylab(parse(text = 'log(CPM + 1)')) +
  theme_light() +
  ylim(c(min(df$E), 1.2*max(df$E))) + 
  theme(plot.title = element_text(face = 4)) +
  xlab('Side') +
  theme(panel.border = element_rect(color = '#d90429', fill = NA, size = 1))


pDesign <- '
@@AAAAA
@@AAAAA
@@AAAAA
@@AAAAA
@@AAAAA
BCDEFVV
GHIJKVV
LMNOPWY
QRSTUXZ
'
PE <- PE + labs(tag = 'E')
PF <- PF + labs(tag = 'F')
PD <- PD + labs(tag = 'G')

png('../Figures/F6.png', width = 3200, height = 3600, res = 300)
PA + LVR + bpList[[1]] + bpList[[2]] + bpList[[3]] + bpList[[4]] + bpList[[5]] + 
  bpList[[6]] + bpList[[7]] + bpList[[8]] + bpList[[9]] + bpList[[10]] + 
  bpList[[11]] + bpList[[12]] + bpList[[13]] + bpList[[14]] + bpList[[15]] + 
  bpList[[16]] + bpList[[17]] + bpList[[18]] + bpList[[19]] + bpList[[20]] + 
  PE + PF + PD + P1 + P2 + plot_layout(design = pDesign)
dev.off()

