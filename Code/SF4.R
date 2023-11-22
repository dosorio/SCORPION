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

doPlots <- function(donorID, tag = FALSE){
  C <- loadNetwork(paste0('~/scCRC/Results/Networks/',donorID,'-T-Epithelial.RData'))
  B <- loadNetwork(paste0('~/scCRC/Results/Networks/',donorID,'-B-Epithelial.RData'))
  N <- loadNetwork(paste0('~/scCRC/Results/Networks/',donorID,'-N-Epithelial.RData'))
  
  gList <- intersect(intersect(colnames(C), colnames(B)), colnames(N))
  df <- data.frame(C = c(as.matrix(C[,gList])), B = c(as.matrix(B[,gList])), N = c(as.matrix(N[,gList])))
  cor(df, method = 'sp')
  
  corValue <- cor(df[,1], df[,2], method = 'sp')
  P <- densCols(df[,1], df[,2], colramp = colorRampPalette(viridis::cividis(5)))
  PA <- ggplot(df, aes(B, C)) + 
    geom_point(pch = 16, alpha = 0.01, color = P) +
    theme_light() +
    geom_abline(intercept = 0, slope = 1, color = 'red', lty = 2) +
    geom_density_2d(color = 'gray90', alpha = 0.5) +
    xlab(parse(text = 'Tumor~Border[Epithelial~Cells]')) +
    ylab(parse(text = 'Tumor~Core[Epithelial~Cells]')) +
    labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
    labs(title = donorID) +
    theme(plot.title = element_text(face = 2)) +
    theme(legend.position = 'None')
  if(tag){
    PA <- PA + labs(tag = 'A')
  }
  
  corValue <- cor(df[,2], df[,3], method = 'sp')
  P <- densCols(df[,2], df[,3], colramp = colorRampPalette(viridis::cividis(5)))
  PB <- ggplot(df, aes(B, N)) + 
    geom_point(color = P, alpha = 0.01, pch = 16) +
    theme_light() +
    geom_abline(intercept = 0, slope = 1, color = 'red', lty = 2) +
    geom_density_2d(color = 'gray90', alpha = 0.5) +
    xlab(parse(text = 'Tumor~Border[Epithelial~Cells]')) +
    ylab(parse(text = 'Normal~Tissue[Epithelial~Cells]')) +
    labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
    scale_color_viridis_c(end = 1) +
    theme(legend.position = 'None')
  if(tag){
    PB <- PB + labs(tag = 'B')
  }
  
  corValue <- cor(df[,1], df[,3], method = 'sp')
  P <- densCols(df[,1], df[,3], colramp = colorRampPalette(viridis::cividis(5)))
  PC <- ggplot(df, aes(C, N, color = P)) + 
    geom_point(color = P, alpha = 0.01, pch = 16) +
    theme_light() +
    geom_abline(intercept = 0, slope = 1, color = 'red', lty = 2) +
    geom_density_2d(color = 'gray90', alpha = 0.5) +
    xlab(parse(text = 'Tumor~Core[Epithelial~Cells]')) +
    ylab(parse(text = 'Normal~Tissue[Epithelial~Cells]')) + 
    labs(subtitle = parse(text = paste0('rho == ', round(corValue,3)))) +
    scale_color_viridis_c(end = 1) +
    theme(legend.position = 'None')
  
  O <- list(PA, PB, PC)  
  return(O)
}


P31 <- doPlots('P31', TRUE)
KUL01 <- doPlots('KUL01', FALSE)
P33 <- doPlots('P33', FALSE)
KUL21 <- doPlots('KUL21', FALSE)

png('../Figures/SF2.png', width = 3000, height = 3000, res = 300)
plotLayout <-
  'ABC
DEF
GHI
JKL'
P31[[1]] + P31[[2]] + P31[[3]] +
  KUL01[[1]] + KUL01[[2]] + KUL01[[3]] +
  P33[[1]] + P33[[2]] + P33[[3]] +
  KUL21[[1]] + KUL21[[2]] + KUL21[[3]] + plot_layout(design = plotLayout)
dev.off()
