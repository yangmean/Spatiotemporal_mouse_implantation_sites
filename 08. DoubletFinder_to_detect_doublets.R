#Author: Min Yang
#We used this package to detect doublets in each sample.

library("DoubletFinder")

#start from one.sample.seurat object
#求最有pK
one.sample.seurat.doublet.list <- paramSweep_v3(one.sample.seurat, PCs = 1:30, sct = F)
one.sample.seurat.doublet.stats <- summarizeSweep(one.sample.seurat.doublet.list, GT = FALSE)
one.sample.seurat.doublet.bcmvn <- find.pK(one.sample.seurat.doublet.stats)
one.sample.seurat.doublet.bcmvn <- one.sample.seurat.doublet.bcmvn$pK[which.max(one.sample.seurat.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#one.sample.seurat.doublet.celltype <- modelHomotypic(one.sample.seurat$subname) 
nExp_poi <- round(doubletrate*ncol(one.sample.seurat)) #predicte the number of doublets
one.sample.seurat <- doubletFinder_v3(one.sample.seurat, PCs = 1:30, pN = 0.25, pK = one.sample.seurat.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
