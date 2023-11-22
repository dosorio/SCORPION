library(Matrix)
library(Seurat)
library(SCORPION)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(harmony)
library(ggridges)
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
library(Nebulosa)


DKO <- Read10X_h5('GSM3477500_Hnf4agDKO.h5')
WT <- Read10X_h5('GSM3477499_Hnf4agWT.h5')

DKO <- DKO[rowMeans(DKO != 0) > 0.05,]
WT <- WT[rowMeans(WT != 0) > 0.05,]

rownames(DKO) <- make.unique(rownames(DKO))
rownames(WT) <- make.unique(rownames(WT))

load('mm10_TF.RData')
load('mm_PPI.RData')

pbapply::pbsapply(1:50, function(i){
  if(!file.exists(paste0('random_hnf4ag_DKO_',i,'.RData'))){
    random_mmTF <- mmTF
    set.seed(i)
    random_mmTF[,1] <- sample(random_mmTF[,1])
    random_mmTF[,2] <- sample(random_mmTF[,2])
    random_mmTF[,3] <- sample(random_mmTF[,3])
    
    random_mmPPI <- mmPPI
    set.seed(i)
    random_mmPPI[,1] <- sample(random_mmPPI[,1])
    random_mmPPI[,2] <- sample(random_mmPPI[,2])
    random_mmPPI[,3] <- sample(random_mmPPI[,3])
    
    random_wtN <- scorpion(tfMotifs = random_mmTF, gexMatrix = WT, ppiNet = random_mmPPI)
    save(random_wtN, file = paste0('random_hnf4ag_WT_',i,'.RData'))
    random_dkoN <- scorpion(tfMotifs = random_mmTF, gexMatrix = DKO, ppiNet = random_mmPPI)
    save(random_dkoN, file = paste0('random_hnf4ag_DKO_',i,'.RData'))
  }
})

random_results <- t(sapply(1:50, function(i){
  load(paste0('random_hnf4ag_DKO_',i,'.RData'))
  load(paste0('random_hnf4ag_WT_',i,'.RData'))
  gList <- intersect(colnames(random_dkoN$regNet), colnames(random_wtN$regNet))
  random_dkoN$regNet <- random_dkoN$regNet[,gList]
  random_wtN$regNet <- random_wtN$regNet[,gList]
  mean_hnf4a <- mean(random_dkoN$regNet['Hnf4a',] - random_wtN$regNet['Hnf4a',])
  mean_hnf4g <- mean(random_dkoN$regNet['Hnf4g',] - random_wtN$regNet['Hnf4g',])
  network_similarity <- cor(c(as.matrix(random_dkoN$regNet)), c(as.matrix(random_wtN$regNet)))
  data.frame(mean_hnf4a = mean_hnf4a, mean_hnf4g = mean_hnf4g, network_similarity = network_similarity)
}))

random_results <- data.frame(mean_hnf4a = unlist(random_results[,1]), 
           mean_hnf4g = unlist(random_results[,2]),
           network_similarity = unlist(random_results[,3]))
random_results <- reshape2::melt(random_results)

png('random_prioris.png', width = 2750, height = 500, res = 300)
ggplot(random_results, aes(value, variable)) +
  ggridges::geom_density_ridges(scale=0.8,quantile_lines = TRUE, quantiles = 2, fill="grey", alpha = 0.2) +
  theme_light() +
  geom_segment(aes(x = -0.24, xend = -0.24, y = 1,yend = 1.9), color = 'red', lty = 2) +
  geom_segment(aes(x = -0.21, xend = -0.21, y = 2,yend = 2.9), color = 'red', lty = 2) + 
  geom_segment(aes(x = 0.88, xend = 0.88, y = 3,yend = 3.9), color = 'red', lty = 2) +
  scale_y_discrete(NULL, labels = c(parse(text = 'hat(mu)(KO[italic(Hnf4a)]-WT[italic(Hnf4a)])'), 
                                   parse(text = 'hat(mu)(KO[italic(Hnf4g)]-WT[italic(Hnf4g)])'), 
                                   parse(text = 'hat(rho)'))) +
  xlab('Estimate')
dev.off()

t.test(random_results[random_results[,1] == 'mean_hnf4a',2], mu = -0.24, alternative = 'greater')
t.test(random_results[random_results[,1] == 'mean_hnf4g',2], mu = -0.21, alternative = 'greater')
t.test(random_results[random_results[,1] == 'network_similarity',2], mu = 0.88, alternative = 'greater')
