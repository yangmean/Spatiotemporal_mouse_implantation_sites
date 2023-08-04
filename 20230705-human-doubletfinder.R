#Nature
human.S.meta <- read.table("../../20220508-download-data/Nature-2018/meta_10x.txt", sep = "\t", header = T, check.names = F)
human.S.data <- read.table("../../20220508-download-data/Nature-2018/raw_data_10x.txt", sep = "\t", header = T, check.names = F)
human.S.data.c <- human.S.data
human.S.data.c$Gene <- gsub("_ENSG.*", "", human.S.data.c$Gene)
human.S.data.c <- data.frame(human.S.data.c %>% group_by(Gene) %>% filter(!duplicated(Gene)), check.names = F)
rownames(human.S.data.c) <- human.S.data.c$Gene
human.S.data.c <- human.S.data.c[, -1]
human.all.seurat <- CreateSeuratObject(counts = human.S.data.c, min.cells = 3, min.features = 200, meta.data = human.S.meta)
save(human.all.seurat, file = "../../20220508-download-data/Nature-2018/Human-nature-2018-All-20230705.RData")

#FCA7167219
FCA7167219.doublet <- subset(human.all.seurat, orig.ident == "FCA7167219")
FCA7167219.doublet <- NormalizeData(FCA7167219.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7167219.doublet <- FindVariableFeatures(FCA7167219.doublet, selection.method = "vst", nfeatures = 2000)
FCA7167219.doublet <- ScaleData(FCA7167219.doublet, features = VariableFeatures(FCA7167219.doublet))
FCA7167219.doublet <- RunPCA(FCA7167219.doublet, features = VariableFeatures(FCA7167219.doublet), npcs = 50)
#求最有pK
FCA7167219.doublet.list <- paramSweep_v3(FCA7167219.doublet, PCs = 1:30, sct = F)
FCA7167219.doublet.stats <- summarizeSweep(FCA7167219.doublet.list, GT = FALSE)
FCA7167219.doublet.bcmvn <- find.pK(FCA7167219.doublet.stats)
FCA7167219.doublet.bcmvn <- FCA7167219.doublet.bcmvn$pK[which.max(FCA7167219.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7167219.doublet.celltype <- modelHomotypic(FCA7167219.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7167219.doublet)) #预估doublet数量
FCA7167219.doublet <- doubletFinder_v3(FCA7167219.doublet, PCs = 1:30, pN = 0.25, pK = FCA7167219.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#FCA7167221
FCA7167221.doublet <- subset(human.all.seurat, orig.ident == "FCA7167221")
FCA7167221.doublet <- NormalizeData(FCA7167221.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7167221.doublet <- FindVariableFeatures(FCA7167221.doublet, selection.method = "vst", nfeatures = 2000)
FCA7167221.doublet <- ScaleData(FCA7167221.doublet, features = VariableFeatures(FCA7167221.doublet))
FCA7167221.doublet <- RunPCA(FCA7167221.doublet, features = VariableFeatures(FCA7167221.doublet), npcs = 50)
#求最有pK
FCA7167221.doublet.list <- paramSweep_v3(FCA7167221.doublet, PCs = 1:30, sct = F)
FCA7167221.doublet.stats <- summarizeSweep(FCA7167221.doublet.list, GT = FALSE)
FCA7167221.doublet.bcmvn <- find.pK(FCA7167221.doublet.stats)
FCA7167221.doublet.bcmvn <- FCA7167221.doublet.bcmvn$pK[which.max(FCA7167221.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7167221.doublet.celltype <- modelHomotypic(FCA7167221.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7167221.doublet)) #预估doublet数量
FCA7167221.doublet <- doubletFinder_v3(FCA7167221.doublet, PCs = 1:30, pN = 0.25, pK = FCA7167221.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#FCA7167222
FCA7167222.doublet <- subset(human.all.seurat, orig.ident == "FCA7167222")
FCA7167222.doublet <- NormalizeData(FCA7167222.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7167222.doublet <- FindVariableFeatures(FCA7167222.doublet, selection.method = "vst", nfeatures = 2000)
FCA7167222.doublet <- ScaleData(FCA7167222.doublet, features = VariableFeatures(FCA7167222.doublet))
FCA7167222.doublet <- RunPCA(FCA7167222.doublet, features = VariableFeatures(FCA7167222.doublet), npcs = 50)
#求最有pK
FCA7167222.doublet.list <- paramSweep_v3(FCA7167222.doublet, PCs = 1:30, sct = F)
FCA7167222.doublet.stats <- summarizeSweep(FCA7167222.doublet.list, GT = FALSE)
FCA7167222.doublet.bcmvn <- find.pK(FCA7167222.doublet.stats)
FCA7167222.doublet.bcmvn <- FCA7167222.doublet.bcmvn$pK[which.max(FCA7167222.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7167222.doublet.celltype <- modelHomotypic(FCA7167222.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7167222.doublet)) #预估doublet数量
FCA7167222.doublet <- doubletFinder_v3(FCA7167222.doublet, PCs = 1:30, pN = 0.25, pK = FCA7167222.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#FCA7167223
FCA7167223.doublet <- subset(human.all.seurat, orig.ident == "FCA7167223")
FCA7167223.doublet <- NormalizeData(FCA7167223.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7167223.doublet <- FindVariableFeatures(FCA7167223.doublet, selection.method = "vst", nfeatures = 2000)
FCA7167223.doublet <- ScaleData(FCA7167223.doublet, features = VariableFeatures(FCA7167223.doublet))
FCA7167223.doublet <- RunPCA(FCA7167223.doublet, features = VariableFeatures(FCA7167223.doublet), npcs = 50)
#求最有pK
FCA7167223.doublet.list <- paramSweep_v3(FCA7167223.doublet, PCs = 1:30, sct = F)
FCA7167223.doublet.stats <- summarizeSweep(FCA7167223.doublet.list, GT = FALSE)
FCA7167223.doublet.bcmvn <- find.pK(FCA7167223.doublet.stats)
FCA7167223.doublet.bcmvn <- FCA7167223.doublet.bcmvn$pK[which.max(FCA7167223.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7167223.doublet.celltype <- modelHomotypic(FCA7167223.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7167223.doublet)) #预估doublet数量
FCA7167223.doublet <- doubletFinder_v3(FCA7167223.doublet, PCs = 1:30, pN = 0.25, pK = FCA7167223.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#FCA7167224
FCA7167224.doublet <- subset(human.all.seurat, orig.ident == "FCA7167224")
FCA7167224.doublet <- NormalizeData(FCA7167224.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7167224.doublet <- FindVariableFeatures(FCA7167224.doublet, selection.method = "vst", nfeatures = 2000)
FCA7167224.doublet <- ScaleData(FCA7167224.doublet, features = VariableFeatures(FCA7167224.doublet))
FCA7167224.doublet <- RunPCA(FCA7167224.doublet, features = VariableFeatures(FCA7167224.doublet), npcs = 50)
#求最有pK
FCA7167224.doublet.list <- paramSweep_v3(FCA7167224.doublet, PCs = 1:30, sct = F)
FCA7167224.doublet.stats <- summarizeSweep(FCA7167224.doublet.list, GT = FALSE)
FCA7167224.doublet.bcmvn <- find.pK(FCA7167224.doublet.stats)
FCA7167224.doublet.bcmvn <- FCA7167224.doublet.bcmvn$pK[which.max(FCA7167224.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7167224.doublet.celltype <- modelHomotypic(FCA7167224.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7167224.doublet)) #预估doublet数量
FCA7167224.doublet <- doubletFinder_v3(FCA7167224.doublet, PCs = 1:30, pN = 0.25, pK = FCA7167224.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#FCA7167226
FCA7167226.doublet <- subset(human.all.seurat, orig.ident == "FCA7167226")
FCA7167226.doublet <- NormalizeData(FCA7167226.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7167226.doublet <- FindVariableFeatures(FCA7167226.doublet, selection.method = "vst", nfeatures = 2000)
FCA7167226.doublet <- ScaleData(FCA7167226.doublet, features = VariableFeatures(FCA7167226.doublet))
FCA7167226.doublet <- RunPCA(FCA7167226.doublet, features = VariableFeatures(FCA7167226.doublet), npcs = 50)
#求最有pK
FCA7167226.doublet.list <- paramSweep_v3(FCA7167226.doublet, PCs = 1:30, sct = F)
FCA7167226.doublet.stats <- summarizeSweep(FCA7167226.doublet.list, GT = FALSE)
FCA7167226.doublet.bcmvn <- find.pK(FCA7167226.doublet.stats)
FCA7167226.doublet.bcmvn <- FCA7167226.doublet.bcmvn$pK[which.max(FCA7167226.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7167226.doublet.celltype <- modelHomotypic(FCA7167226.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7167226.doublet)) #预估doublet数量
FCA7167226.doublet <- doubletFinder_v3(FCA7167226.doublet, PCs = 1:30, pN = 0.25, pK = FCA7167226.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#FCA7196218
FCA7196218.doublet <- subset(human.all.seurat, orig.ident == "FCA7196218")
FCA7196218.doublet <- NormalizeData(FCA7196218.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7196218.doublet <- FindVariableFeatures(FCA7196218.doublet, selection.method = "vst", nfeatures = 2000)
FCA7196218.doublet <- ScaleData(FCA7196218.doublet, features = VariableFeatures(FCA7196218.doublet))
FCA7196218.doublet <- RunPCA(FCA7196218.doublet, features = VariableFeatures(FCA7196218.doublet), npcs = 50)
#求最有pK
FCA7196218.doublet.list <- paramSweep_v3(FCA7196218.doublet, PCs = 1:30, sct = F)
FCA7196218.doublet.stats <- summarizeSweep(FCA7196218.doublet.list, GT = FALSE)
FCA7196218.doublet.bcmvn <- find.pK(FCA7196218.doublet.stats)
FCA7196218.doublet.bcmvn <- FCA7196218.doublet.bcmvn$pK[which.max(FCA7196218.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7196218.doublet.celltype <- modelHomotypic(FCA7196218.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7196218.doublet)) #预估doublet数量
FCA7196218.doublet <- doubletFinder_v3(FCA7196218.doublet, PCs = 1:30, pN = 0.25, pK = FCA7196218.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#FCA7196219
FCA7196219.doublet <- subset(human.all.seurat, orig.ident == "FCA7196219")
FCA7196219.doublet <- NormalizeData(FCA7196219.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7196219.doublet <- FindVariableFeatures(FCA7196219.doublet, selection.method = "vst", nfeatures = 2000)
FCA7196219.doublet <- ScaleData(FCA7196219.doublet, features = VariableFeatures(FCA7196219.doublet))
FCA7196219.doublet <- RunPCA(FCA7196219.doublet, features = VariableFeatures(FCA7196219.doublet), npcs = 50)
#求最有pK
FCA7196219.doublet.list <- paramSweep_v3(FCA7196219.doublet, PCs = 1:30, sct = F)
FCA7196219.doublet.stats <- summarizeSweep(FCA7196219.doublet.list, GT = FALSE)
FCA7196219.doublet.bcmvn <- find.pK(FCA7196219.doublet.stats)
FCA7196219.doublet.bcmvn <- FCA7196219.doublet.bcmvn$pK[which.max(FCA7196219.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7196219.doublet.celltype <- modelHomotypic(FCA7196219.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7196219.doublet)) #预估doublet数量
FCA7196219.doublet <- doubletFinder_v3(FCA7196219.doublet, PCs = 1:30, pN = 0.25, pK = FCA7196219.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#FCA7196224
FCA7196224.doublet <- subset(human.all.seurat, orig.ident == "FCA7196224")
FCA7196224.doublet <- NormalizeData(FCA7196224.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7196224.doublet <- FindVariableFeatures(FCA7196224.doublet, selection.method = "vst", nfeatures = 2000)
FCA7196224.doublet <- ScaleData(FCA7196224.doublet, features = VariableFeatures(FCA7196224.doublet))
FCA7196224.doublet <- RunPCA(FCA7196224.doublet, features = VariableFeatures(FCA7196224.doublet), npcs = 50)
#求最有pK
FCA7196224.doublet.list <- paramSweep_v3(FCA7196224.doublet, PCs = 1:30, sct = F)
FCA7196224.doublet.stats <- summarizeSweep(FCA7196224.doublet.list, GT = FALSE)
FCA7196224.doublet.bcmvn <- find.pK(FCA7196224.doublet.stats)
FCA7196224.doublet.bcmvn <- FCA7196224.doublet.bcmvn$pK[which.max(FCA7196224.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7196224.doublet.celltype <- modelHomotypic(FCA7196224.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7196224.doublet)) #预估doublet数量
FCA7196224.doublet <- doubletFinder_v3(FCA7196224.doublet, PCs = 1:30, pN = 0.25, pK = FCA7196224.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#FCA7196225
FCA7196225.doublet <- subset(human.all.seurat, orig.ident == "FCA7196225")
FCA7196225.doublet <- NormalizeData(FCA7196225.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7196225.doublet <- FindVariableFeatures(FCA7196225.doublet, selection.method = "vst", nfeatures = 2000)
FCA7196225.doublet <- ScaleData(FCA7196225.doublet, features = VariableFeatures(FCA7196225.doublet))
FCA7196225.doublet <- RunPCA(FCA7196225.doublet, features = VariableFeatures(FCA7196225.doublet), npcs = 50)
#求最有pK
FCA7196225.doublet.list <- paramSweep_v3(FCA7196225.doublet, PCs = 1:30, sct = F)
FCA7196225.doublet.stats <- summarizeSweep(FCA7196225.doublet.list, GT = FALSE)
FCA7196225.doublet.bcmvn <- find.pK(FCA7196225.doublet.stats)
FCA7196225.doublet.bcmvn <- FCA7196225.doublet.bcmvn$pK[which.max(FCA7196225.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7196225.doublet.celltype <- modelHomotypic(FCA7196225.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7196225.doublet)) #预估doublet数量
FCA7196225.doublet <- doubletFinder_v3(FCA7196225.doublet, PCs = 1:30, pN = 0.25, pK = FCA7196225.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#FCA7196229
FCA7196229.doublet <- subset(human.all.seurat, orig.ident == "FCA7196229")
FCA7196229.doublet <- NormalizeData(FCA7196229.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7196229.doublet <- FindVariableFeatures(FCA7196229.doublet, selection.method = "vst", nfeatures = 2000)
FCA7196229.doublet <- ScaleData(FCA7196229.doublet, features = VariableFeatures(FCA7196229.doublet))
FCA7196229.doublet <- RunPCA(FCA7196229.doublet, features = VariableFeatures(FCA7196229.doublet), npcs = 50)
#求最有pK
FCA7196229.doublet.list <- paramSweep_v3(FCA7196229.doublet, PCs = 1:30, sct = F)
FCA7196229.doublet.stats <- summarizeSweep(FCA7196229.doublet.list, GT = FALSE)
FCA7196229.doublet.bcmvn <- find.pK(FCA7196229.doublet.stats)
FCA7196229.doublet.bcmvn <- FCA7196229.doublet.bcmvn$pK[which.max(FCA7196229.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7196229.doublet.celltype <- modelHomotypic(FCA7196229.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7196229.doublet)) #预估doublet数量
FCA7196229.doublet <- doubletFinder_v3(FCA7196229.doublet, PCs = 1:30, pN = 0.25, pK = FCA7196229.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#FCA7474062
FCA7474062.doublet <- subset(human.all.seurat, orig.ident == "FCA7474062")
FCA7474062.doublet <- NormalizeData(FCA7474062.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7474062.doublet <- FindVariableFeatures(FCA7474062.doublet, selection.method = "vst", nfeatures = 2000)
FCA7474062.doublet <- ScaleData(FCA7474062.doublet, features = VariableFeatures(FCA7474062.doublet))
FCA7474062.doublet <- RunPCA(FCA7474062.doublet, features = VariableFeatures(FCA7474062.doublet), npcs = 50)
#求最有pK
FCA7474062.doublet.list <- paramSweep_v3(FCA7474062.doublet, PCs = 1:30, sct = F)
FCA7474062.doublet.stats <- summarizeSweep(FCA7474062.doublet.list, GT = FALSE)
FCA7474062.doublet.bcmvn <- find.pK(FCA7474062.doublet.stats)
FCA7474062.doublet.bcmvn <- FCA7474062.doublet.bcmvn$pK[which.max(FCA7474062.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7474062.doublet.celltype <- modelHomotypic(FCA7474062.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7474062.doublet)) #预估doublet数量
FCA7474062.doublet <- doubletFinder_v3(FCA7474062.doublet, PCs = 1:30, pN = 0.25, pK = FCA7474062.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#FCA7474063
FCA7474063.doublet <- subset(human.all.seurat, orig.ident == "FCA7474063")
FCA7474063.doublet <- NormalizeData(FCA7474063.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7474063.doublet <- FindVariableFeatures(FCA7474063.doublet, selection.method = "vst", nfeatures = 2000)
FCA7474063.doublet <- ScaleData(FCA7474063.doublet, features = VariableFeatures(FCA7474063.doublet))
FCA7474063.doublet <- RunPCA(FCA7474063.doublet, features = VariableFeatures(FCA7474063.doublet), npcs = 50)
#求最有pK
FCA7474063.doublet.list <- paramSweep_v3(FCA7474063.doublet, PCs = 1:30, sct = F)
FCA7474063.doublet.stats <- summarizeSweep(FCA7474063.doublet.list, GT = FALSE)
FCA7474063.doublet.bcmvn <- find.pK(FCA7474063.doublet.stats)
FCA7474063.doublet.bcmvn <- FCA7474063.doublet.bcmvn$pK[which.max(FCA7474063.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7474063.doublet.celltype <- modelHomotypic(FCA7474063.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7474063.doublet)) #预估doublet数量
FCA7474063.doublet <- doubletFinder_v3(FCA7474063.doublet, PCs = 1:30, pN = 0.25, pK = FCA7474063.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#FCA7511881
FCA7511881.doublet <- subset(human.all.seurat, orig.ident == "FCA7511881")
FCA7511881.doublet <- NormalizeData(FCA7511881.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7511881.doublet <- FindVariableFeatures(FCA7511881.doublet, selection.method = "vst", nfeatures = 2000)
FCA7511881.doublet <- ScaleData(FCA7511881.doublet, features = VariableFeatures(FCA7511881.doublet))
FCA7511881.doublet <- RunPCA(FCA7511881.doublet, features = VariableFeatures(FCA7511881.doublet), npcs = 50)
#求最有pK
FCA7511881.doublet.list <- paramSweep_v3(FCA7511881.doublet, PCs = 1:30, sct = F)
FCA7511881.doublet.stats <- summarizeSweep(FCA7511881.doublet.list, GT = FALSE)
FCA7511881.doublet.bcmvn <- find.pK(FCA7511881.doublet.stats)
FCA7511881.doublet.bcmvn <- FCA7511881.doublet.bcmvn$pK[which.max(FCA7511881.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7511881.doublet.celltype <- modelHomotypic(FCA7511881.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7511881.doublet)) #预估doublet数量
FCA7511881.doublet <- doubletFinder_v3(FCA7511881.doublet, PCs = 1:30, pN = 0.25, pK = FCA7511881.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#FCA7511882
FCA7511882.doublet <- subset(human.all.seurat, orig.ident == "FCA7511882")
FCA7511882.doublet <- NormalizeData(FCA7511882.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7511882.doublet <- FindVariableFeatures(FCA7511882.doublet, selection.method = "vst", nfeatures = 2000)
FCA7511882.doublet <- ScaleData(FCA7511882.doublet, features = VariableFeatures(FCA7511882.doublet))
FCA7511882.doublet <- RunPCA(FCA7511882.doublet, features = VariableFeatures(FCA7511882.doublet), npcs = 50)
#求最有pK
FCA7511882.doublet.list <- paramSweep_v3(FCA7511882.doublet, PCs = 1:30, sct = F)
FCA7511882.doublet.stats <- summarizeSweep(FCA7511882.doublet.list, GT = FALSE)
FCA7511882.doublet.bcmvn <- find.pK(FCA7511882.doublet.stats)
FCA7511882.doublet.bcmvn <- FCA7511882.doublet.bcmvn$pK[which.max(FCA7511882.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7511882.doublet.celltype <- modelHomotypic(FCA7511882.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7511882.doublet)) #预估doublet数量
FCA7511882.doublet <- doubletFinder_v3(FCA7511882.doublet, PCs = 1:30, pN = 0.25, pK = FCA7511882.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#FCA7511884
FCA7511884.doublet <- subset(human.all.seurat, orig.ident == "FCA7511884")
FCA7511884.doublet <- NormalizeData(FCA7511884.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
FCA7511884.doublet <- FindVariableFeatures(FCA7511884.doublet, selection.method = "vst", nfeatures = 2000)
FCA7511884.doublet <- ScaleData(FCA7511884.doublet, features = VariableFeatures(FCA7511884.doublet))
FCA7511884.doublet <- RunPCA(FCA7511884.doublet, features = VariableFeatures(FCA7511884.doublet), npcs = 50)
#求最有pK
FCA7511884.doublet.list <- paramSweep_v3(FCA7511884.doublet, PCs = 1:30, sct = F)
FCA7511884.doublet.stats <- summarizeSweep(FCA7511884.doublet.list, GT = FALSE)
FCA7511884.doublet.bcmvn <- find.pK(FCA7511884.doublet.stats)
FCA7511884.doublet.bcmvn <- FCA7511884.doublet.bcmvn$pK[which.max(FCA7511884.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#FCA7511884.doublet.celltype <- modelHomotypic(FCA7511884.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(FCA7511884.doublet)) #预估doublet数量
FCA7511884.doublet <- doubletFinder_v3(FCA7511884.doublet, PCs = 1:30, pN = 0.25, pK = FCA7511884.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#extract the doublet information
human.all.doublet.meta <- c(FCA7167219.doublet$DF.classifications_0.25_0.12_23, FCA7167221.doublet$DF.classifications_0.25_0.03_30, 
                            FCA7167222.doublet$DF.classifications_0.25_0.26_59, FCA7167223.doublet$DF.classifications_0.25_0.14_229, 
                            FCA7167224.doublet$DF.classifications_0.25_0.01_86, FCA7167226.doublet$DF.classifications_0.25_0.23_270, 
                            FCA7196218.doublet$DF.classifications_0.25_0.24_290, FCA7196219.doublet$DF.classifications_0.25_0.12_172, 
                            FCA7196224.doublet$DF.classifications_0.25_0.25_345, FCA7196225.doublet$DF.classifications_0.25_0.3_131, 
                            FCA7196229.doublet$DF.classifications_0.25_0.01_165, FCA7474062.doublet$DF.classifications_0.25_0.26_96, 
                            FCA7474063.doublet$DF.classifications_0.25_0.08_67, FCA7511881.doublet$DF.classifications_0.25_0.01_88, 
                            FCA7511882.doublet$DF.classifications_0.25_0.13_103, FCA7511884.doublet$DF.classifications_0.25_0.24_124)
save(human.all.doublet.meta, file = "../../20220508-download-data/Nature-2018/Human-all-doubletfinder-info.RData")

load("20220622-Comments/02.human-data/nature/Nature-human-DSC-recluster-CCA.RData")
DimPlot(human.S.combined)
human.all.doublet.meta.sub <- human.all.doublet.meta[colnames(human.S.combined)]
human.S.combined <- AddMetaData(human.S.combined, metadata = human.all.doublet.meta.sub, col.name = "Doublet")
DimPlot(human.S.combined, group.by = "Doublet", pt.size = 0.5)
pdf("20220622-Comments/02.human-data/00.figures/21.Human-nature-doubletfinder.pdf", width = 10, height = 10)
DimPlot(human.S.combined, group.by = "Doublet", pt.size = 0.5)
dev.off()

###################SCIENCE ADVANCE#########################################
load("20220622-Comments/02.human-data/science_advance/decidua_commbined.RData")
DimPlot(decidua.combined)

#P1D.DS
P1D.DS.doublet <- subset(decidua.combined, orig.ident == "P1D.DS")
P1D.DS.doublet <- NormalizeData(P1D.DS.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
P1D.DS.doublet <- FindVariableFeatures(P1D.DS.doublet, selection.method = "vst", nfeatures = 2000)
P1D.DS.doublet <- ScaleData(P1D.DS.doublet, features = VariableFeatures(P1D.DS.doublet))
P1D.DS.doublet <- RunPCA(P1D.DS.doublet, features = VariableFeatures(P1D.DS.doublet), npcs = 50)
#求最有pK
P1D.DS.doublet.list <- paramSweep_v3(P1D.DS.doublet, PCs = 1:30, sct = F)
P1D.DS.doublet.stats <- summarizeSweep(P1D.DS.doublet.list, GT = FALSE)
P1D.DS.doublet.bcmvn <- find.pK(P1D.DS.doublet.stats)
P1D.DS.doublet.bcmvn <- P1D.DS.doublet.bcmvn$pK[which.max(P1D.DS.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#P1D.DS.doublet.celltype <- modelHomotypic(P1D.DS.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(P1D.DS.doublet)) #预估doublet数量
P1D.DS.doublet <- doubletFinder_v3(P1D.DS.doublet, PCs = 1:30, pN = 0.25, pK = P1D.DS.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#P2D.DS
P2D.DS.doublet <- subset(decidua.combined, orig.ident == "P2D.DS")
P2D.DS.doublet <- NormalizeData(P2D.DS.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
P2D.DS.doublet <- FindVariableFeatures(P2D.DS.doublet, selection.method = "vst", nfeatures = 2000)
P2D.DS.doublet <- ScaleData(P2D.DS.doublet, features = VariableFeatures(P2D.DS.doublet))
P2D.DS.doublet <- RunPCA(P2D.DS.doublet, features = VariableFeatures(P2D.DS.doublet), npcs = 50)
#求最有pK
P2D.DS.doublet.list <- paramSweep_v3(P2D.DS.doublet, PCs = 1:30, sct = F)
P2D.DS.doublet.stats <- summarizeSweep(P2D.DS.doublet.list, GT = FALSE)
P2D.DS.doublet.bcmvn <- find.pK(P2D.DS.doublet.stats)
P2D.DS.doublet.bcmvn <- P2D.DS.doublet.bcmvn$pK[which.max(P2D.DS.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#P2D.DS.doublet.celltype <- modelHomotypic(P2D.DS.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(P2D.DS.doublet)) #预估doublet数量
P2D.DS.doublet <- doubletFinder_v3(P2D.DS.doublet, PCs = 1:30, pN = 0.25, pK = P2D.DS.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#P3D.DS
P3D.DS.doublet <- subset(decidua.combined, orig.ident == "P3D.DS")
P3D.DS.doublet <- NormalizeData(P3D.DS.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
P3D.DS.doublet <- FindVariableFeatures(P3D.DS.doublet, selection.method = "vst", nfeatures = 2000)
P3D.DS.doublet <- ScaleData(P3D.DS.doublet, features = VariableFeatures(P3D.DS.doublet))
P3D.DS.doublet <- RunPCA(P3D.DS.doublet, features = VariableFeatures(P3D.DS.doublet), npcs = 50)
#求最有pK
P3D.DS.doublet.list <- paramSweep_v3(P3D.DS.doublet, PCs = 1:30, sct = F)
P3D.DS.doublet.stats <- summarizeSweep(P3D.DS.doublet.list, GT = FALSE)
P3D.DS.doublet.bcmvn <- find.pK(P3D.DS.doublet.stats)
P3D.DS.doublet.bcmvn <- P3D.DS.doublet.bcmvn$pK[which.max(P3D.DS.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#P3D.DS.doublet.celltype <- modelHomotypic(P3D.DS.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(P3D.DS.doublet)) #预估doublet数量
P3D.DS.doublet <- doubletFinder_v3(P3D.DS.doublet, PCs = 1:30, pN = 0.25, pK = P3D.DS.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#P4D.DS
P4D.DS.doublet <- subset(decidua.combined, orig.ident == "P4D.DS")
P4D.DS.doublet <- NormalizeData(P4D.DS.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
P4D.DS.doublet <- FindVariableFeatures(P4D.DS.doublet, selection.method = "vst", nfeatures = 2000)
P4D.DS.doublet <- ScaleData(P4D.DS.doublet, features = VariableFeatures(P4D.DS.doublet))
P4D.DS.doublet <- RunPCA(P4D.DS.doublet, features = VariableFeatures(P4D.DS.doublet), npcs = 50)
#求最有pK
P4D.DS.doublet.list <- paramSweep_v3(P4D.DS.doublet, PCs = 1:30, sct = F)
P4D.DS.doublet.stats <- summarizeSweep(P4D.DS.doublet.list, GT = FALSE)
P4D.DS.doublet.bcmvn <- find.pK(P4D.DS.doublet.stats)
P4D.DS.doublet.bcmvn <- P4D.DS.doublet.bcmvn$pK[which.max(P4D.DS.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#P4D.DS.doublet.celltype <- modelHomotypic(P4D.DS.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(P4D.DS.doublet)) #预估doublet数量
P4D.DS.doublet <- doubletFinder_v3(P4D.DS.doublet, PCs = 1:30, pN = 0.25, pK = P4D.DS.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#P5D.DS
P5D.DS.doublet <- subset(decidua.combined, orig.ident == "P5D.DS")
P5D.DS.doublet <- NormalizeData(P5D.DS.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
P5D.DS.doublet <- FindVariableFeatures(P5D.DS.doublet, selection.method = "vst", nfeatures = 2000)
P5D.DS.doublet <- ScaleData(P5D.DS.doublet, features = VariableFeatures(P5D.DS.doublet))
P5D.DS.doublet <- RunPCA(P5D.DS.doublet, features = VariableFeatures(P5D.DS.doublet), npcs = 50)
#求最有pK
P5D.DS.doublet.list <- paramSweep_v3(P5D.DS.doublet, PCs = 1:30, sct = F)
P5D.DS.doublet.stats <- summarizeSweep(P5D.DS.doublet.list, GT = FALSE)
P5D.DS.doublet.bcmvn <- find.pK(P5D.DS.doublet.stats)
P5D.DS.doublet.bcmvn <- P5D.DS.doublet.bcmvn$pK[which.max(P5D.DS.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#P5D.DS.doublet.celltype <- modelHomotypic(P5D.DS.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(P5D.DS.doublet)) #预估doublet数量
P5D.DS.doublet <- doubletFinder_v3(P5D.DS.doublet, PCs = 1:30, pN = 0.25, pK = P5D.DS.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#P6D.10x
P6D.10x.doublet <- subset(decidua.combined, orig.ident == "P6D.10x")
P6D.10x.doublet <- NormalizeData(P6D.10x.doublet, normalization.method = "LogNormalize", scale.factor = 10000)
P6D.10x.doublet <- FindVariableFeatures(P6D.10x.doublet, selection.method = "vst", nfeatures = 2000)
P6D.10x.doublet <- ScaleData(P6D.10x.doublet, features = VariableFeatures(P6D.10x.doublet))
P6D.10x.doublet <- RunPCA(P6D.10x.doublet, features = VariableFeatures(P6D.10x.doublet), npcs = 50)
#求最有pK
P6D.10x.doublet.list <- paramSweep_v3(P6D.10x.doublet, PCs = 1:30, sct = F)
P6D.10x.doublet.stats <- summarizeSweep(P6D.10x.doublet.list, GT = FALSE)
P6D.10x.doublet.bcmvn <- find.pK(P6D.10x.doublet.stats)
P6D.10x.doublet.bcmvn <- P6D.10x.doublet.bcmvn$pK[which.max(P6D.10x.doublet.bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
doubletrate <- 0.055 ## Assuming 5.5% doublet formation rate when loading 12,000 cells - tailor for your dataset
#P6D.10x.doublet.celltype <- modelHomotypic(P6D.10x.doublet$cellname) #不提供细胞类型
nExp_poi <- round(doubletrate*ncol(P6D.10x.doublet)) #预估doublet数量
P6D.10x.doublet <- doubletFinder_v3(P6D.10x.doublet, PCs = 1:30, pN = 0.25, pK = P6D.10x.doublet.bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


#extract the doublet information
human.SA.all.doublet.meta <- c(P1D.DS.doublet$DF.classifications_0.25_0.05_19, P2D.DS.doublet$DF.classifications_0.25_0.16_24, 
                               P3D.DS.doublet$DF.classifications_0.25_0.3_45, P4D.DS.doublet$DF.classifications_0.25_0.3_57, 
                               P5D.DS.doublet$DF.classifications_0.25_0.04_58, P6D.10x.doublet$DF.classifications_0.25_0.005_244)
save(human.SA.all.doublet.meta, file = "../../20220508-download-data/Nature-2018/human-SA-all-doubletfinder-info.RData")

load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-CCA.RData")
DimPlot(FB.DSC.decidua.combined)
human.SA.all.doublet.meta.sub <- human.SA.all.doublet.meta[colnames(FB.DSC.decidua.combined)]
FB.DSC.decidua.combined <- AddMetaData(FB.DSC.decidua.combined, metadata = human.SA.all.doublet.meta.sub, col.name = "Doublet")
pdf("20220622-Comments/02.human-data/00.figures/21.human-SA-nature-doubletfinder.pdf", width = 10, height = 10)
DimPlot(FB.DSC.decidua.combined, group.by = "Doublet", pt.size = 2)
dev.off()

#merged sample
load("20220622-Comments/02.human-data/Human-Nature-SD-DSCs-downsample-500-integrate.RData")
DimPlot(human.decidua.combined.sub)
human.decidua.combined.sub.doublet.meta <- c(human.SA.all.doublet.meta, human.all.doublet.meta)
human.decidua.combined.sub.doublet <- human.decidua.combined.sub.doublet.meta[colnames(human.decidua.combined.sub)]
human.decidua.combined.sub <- AddMetaData(human.decidua.combined.sub, metadata = human.decidua.combined.sub.doublet, col.name = "Doublet")
pdf("20220622-Comments/02.human-data/00.figures/22.Human-Nature-SA-integration-doublet.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub, group.by = "Doublet", pt.size = 1)
dev.off()
#remove doublet
human.decidua.combined.sub.rm <- subset(human.decidua.combined.sub, Doublet == "Singlet")
DefaultAssay(human.decidua.combined.sub.rm) <- "integrated"
human.decidua.combined.sub.rm <- ScaleData(human.decidua.combined.sub.rm, features = rownames(human.decidua.combined.sub.rm))
human.decidua.combined.sub.rm <- RunPCA(human.decidua.combined.sub.rm, features = VariableFeatures(human.decidua.combined.sub.rm))
ElbowPlot(human.decidua.combined.sub.rm, ndims = 50)
human.decidua.combined.sub.rm <- RunUMAP(human.decidua.combined.sub.rm, reduction = "pca", dims = 1:30)
human.decidua.combined.sub.rm <- FindNeighbors(human.decidua.combined.sub.rm, reduction = "pca", dims = 1:30)
human.decidua.combined.sub.rm <- FindClusters(human.decidua.combined.sub.rm, resolution = 0.7)
DimPlot(human.decidua.combined.sub.rm, label = T)
DimPlot(human.decidua.combined.sub.rm, group.by = "literature")
DimPlot(human.decidua.combined.sub.rm, group.by = "orig.cellname")
human.decidua.combined.sub.rm <- RenameIdents(human.decidua.combined.sub.rm, "0" = "hFB1", "1" = "hFB2", "2" = "hFB1", "3" = "hDSC", "5" = "hFB2", 
                                              "4" = "hDSC", "6" = "hDSC", "7" = "hFB1", "8" = "hiDSC")
pdf("20220622-Comments/02.human-data/00.figures/23.Human-DSC-integrate-UMAP-sub1000-rm-doublet-color-cluster.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub.rm, label = F, pt.size = 1, cols = c("#57BBC2", "#C0A335", "#006400", "#6A3D9A"))
dev.off()
pdf("20220622-Comments/02.human-data/00.figures/24.Human-DSC-integrate-UMAP-sub1000-rm-doublet-color-literature.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub.rm, label = F, pt.size = 1, group.by = "literature", cols = c("#9370DB", "#3CB371"))
dev.off()
pdf("20220622-Comments/02.human-data/00.figures/25.Human-DSC-integrate-UMAP-sub1000-rm-doublet-color-orig-cellname", width = 10, height = 10)
DimPlot(human.decidua.combined.sub.rm, label = F, pt.size = 1, group.by = "orig.cellname")
dev.off()
DefaultAssay(human.decidua.combined.sub.rm) <- "RNA"
human.decidua.combined.sub.rm.markers <- FindAllMarkers(human.decidua.combined.sub.rm, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
human.decidua.combined.sub.rm.markers <- human.decidua.combined.sub.rm.markers[human.decidua.combined.sub.rm.markers$gene != "AC090498.1", ]
human.decidua.combined.sub.rm.markers <- human.decidua.combined.sub.rm.markers[human.decidua.combined.sub.rm.markers$gene != "MALAT1", ]
human.decidua.combined.sub.rm.markers.top10 <- human.decidua.combined.sub.rm.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
write.table(human.decidua.combined.sub.rm.markers, file = "20220622-Comments/02.human-data/Human-DSC-subclusters-rm-doublet-FB1-FB2-DSC-DEGs-fc0.25.csv", sep = ",", col.names = T, row.names = F, quote = F)

pdf("20220622-Comments/02.human-data/00.figures/25.Human-DSC-sub-rm-doublet-DEGs-top10.pdf", width = 10, height = 10)
DoHeatmap(human.decidua.combined.sub.rm, features = human.decidua.combined.sub.rm.markers.top10$gene)+ 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()

#GOs
library("clusterProfiler")
enrich.GO.function <- function(geneset){
  enrich.GO <- enrichGO(gene = geneset,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "SYMBOL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
  enrich.GO <- simplify(enrich.GO, cutoff = 0.7, by = "p.adjust", select_fun = min)
  return(enrich.GO)
}
human.iDSC.GO.df <- data.frame()
for (c in unique(human.decidua.combined.sub.rm.markers$cluster)) {
  print(c)
  geneset <- human.decidua.combined.sub.rm.markers[human.decidua.combined.sub.rm.markers$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  human.iDSC.GO.df <- rbind(human.iDSC.GO.df, tmp)
}
human.iDSC.GO.df$GO_Term <- factor(human.iDSC.GO.df$GO_Term, levels = rev(unique(human.iDSC.GO.df$GO_Term)))
write.table(human.iDSC.GO.df, file = "20220622-Comments/02.human-data/Human-DSC-subclusters-rm-doublet-FB1-FB2-DSC-GOs.csv", sep = ",", col.names = T, row.names = F, quote = F)

#chemokine
human.DSC.iDSC.avg <- AverageExpression(human.decidua.combined.sub.rm, slot = "data", assays = "RNA", use.counts = T) #use count matrix
human.DSC.iDSC.avg <- human.DSC.iDSC.avg$RNA
human.DSC.iDSC.avg.chemokine <- human.DSC.iDSC.avg[toupper(c(chemokine.ligands[1:7], "CCL27", chemokine.ligands[9:15], "CCL21", chemokine.ligands[17:18], "Mdk", "Lgals9")),]
color <- rev(brewer.pal(11, "Spectral"))
pdf("20220622-Comments/02.human-data/00.figures/25.human-DSC-iDSC-rm-doublet-chemokines-MK-GALECTIN-L-heatmap.pdf", width = 10, height = 10)
pheatmap(log2(human.DSC.iDSC.avg.chemokine+1), scale = "none", cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(colors = color)(100))
dev.off()

vascular.gene <- c("CALCA", "ADM", "ANGPT2", "ANGPTL2", "ANGPTL4", "RARRES2", "SEMA3A", "SEMA3B", "SEMA3C", "SEMA3E", "VEGFA", "VEGFB", "VEGFD")
pdf("20220622-Comments/02.human-data/00.figures/26.Human-DSC-iDSC-rm-doublet-Vascularization-expression-vlnplot.pdf", width = 10, height = 10)
VlnPlot(human.decidua.combined.sub.rm, features = vascular.gene, pt.size = 0, cols = c("#57BBC2", "#C0A335", "#006400", "#6A3D9A"))
dev.off()

#subset iDSC
human.decidua.combined.sub.iDSC <- subset(human.decidua.combined.sub, seurat_clusters == 4)
DefaultAssay(human.decidua.combined.sub.iDSC) <- "RNA"
human.decidua.combined.sub.iDSC <- NormalizeData(human.decidua.combined.sub.iDSC, normalization.method = "LogNormalize", scale.factor = 10000)
human.decidua.combined.sub.iDSC <- FindVariableFeatures(human.decidua.combined.sub.iDSC, selection.method = "vst", nfeatures = 2000)
human.decidua.combined.sub.iDSC <- ScaleData(human.decidua.combined.sub.iDSC, features = VariableFeatures(human.decidua.combined.sub.iDSC))
human.decidua.combined.sub.iDSC <- RunPCA(human.decidua.combined.sub.iDSC, features = VariableFeatures(human.decidua.combined.sub.iDSC), npcs = 50)

human.decidua.combined.sub.iDSC <- FindNeighbors(object = human.decidua.combined.sub.iDSC, dims = 1:30)
human.decidua.combined.sub.iDSC <- FindClusters(human.decidua.combined.sub.iDSC, resolution = 0.4)
human.decidua.combined.sub.iDSC <- RunUMAP(human.decidua.combined.sub.iDSC, dims = 1:30, seed.use = 417)
DimPlot(human.decidua.combined.sub.iDSC, label = T, group.by = "literature")

human.decidua.combined.sub.iDSC.list <- SplitObject(human.decidua.combined.sub.iDSC, split.by = "literature")
features <- SelectIntegrationFeatures(object.list = human.decidua.combined.sub.iDSC.list)
human.decidua.anchors <- FindIntegrationAnchors(object.list = human.decidua.combined.sub.iDSC.list, anchor.features = features, k.filter = 120)
human.decidua.combined.sub.iDSC <- IntegrateData(anchorset = human.decidua.anchors)
DefaultAssay(human.decidua.combined.sub.iDSC) <- "integrated"
human.decidua.combined.sub.iDSC <- ScaleData(human.decidua.combined.sub.iDSC, features = rownames(human.decidua.combined.sub.iDSC))
human.decidua.combined.sub.iDSC <- RunPCA(human.decidua.combined.sub.iDSC, features = VariableFeatures(human.decidua.combined.sub.iDSC))
ElbowPlot(human.decidua.combined.sub.iDSC, ndims = 50)
human.decidua.combined.sub.iDSC <- RunUMAP(human.decidua.combined.sub.iDSC, reduction = "pca", dims = 1:10)
human.decidua.combined.sub.iDSC <- FindNeighbors(human.decidua.combined.sub.iDSC, reduction = "pca", dims = 1:10)
human.decidua.combined.sub.iDSC <- FindClusters(human.decidua.combined.sub.iDSC, resolution = 0.3)
DimPlot(human.decidua.combined.sub.iDSC, label = T, group.by = "literature")
pdf("20220622-Comments/02.human-data/00.figures/27.human-DSC-recluster-color-doublet.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub.iDSC, label = F, group.by = "Doublet", pt.size = 4)
dev.off()

pdf("20220622-Comments/02.human-data/00.figures/27.human-DSC-recluster-color-cluster.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub.iDSC, cols = c("#3CB371", "#F4A460"), pt.size = 4)
dev.off()

DefaultAssay(human.decidua.combined.sub.iDSC) <- "RNA"
human.decidua.combined.sub.iDSC.markers <- FindAllMarkers(human.decidua.combined.sub.iDSC, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)

enrich.GO.function <- function(geneset){
  enrich.GO <- enrichGO(gene = geneset,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "SYMBOL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
  enrich.GO <- simplify(enrich.GO, cutoff = 0.7, by = "p.adjust", select_fun = min)
  enrich.GO <- enrich.GO@result[enrich.GO@result$pvalue<0.01, ]
  return(enrich.GO)
}
human.decidua.combined.markers.GO.df <- data.frame()
for (c in unique(human.decidua.combined.markers$cluster)) {
  print(c)
  geneset <- human.decidua.combined.markers[human.decidua.combined.markers$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  human.decidua.combined.markers.GO.df <- rbind(human.decidua.combined.markers.GO.df, tmp)
}
human.decidua.combined.markers.GO.df$GO_Term <- factor(human.decidua.combined.markers.GO.df$GO_Term, levels = rev(unique(human.decidua.combined.markers.GO.df$GO_Term)))

