#Author: Min Yang
#Download the human data first:
#Human data published:Normal samples
#1). Vento-Tormo, R., et al. Single-cell reconstruction of the early maternal–fetal interface in humans. 2018
#2). Suryawanshi, H., et al. A single-cell survey of the human first-trimester placenta and decidua. 2018

#Human data published:Disease samples
#1). Du, L., et al. Single-cell transcriptome analysis reveals defective decidua stronal niche attributes to recurrent spontaneous abortion. 2021

library(dplyr)
library(Seurat)
library(harmony)
Human.N.meta <- read.table("./Nature-2018/meta_10x.txt", sep = "\t", header = T, check.names = F)
Human.N.meta.dsc <- Human.N.meta[Human.N.meta$final_cluster %in% c(0, 3, 21), ]
Human.N.data <- read.table("./Nature-2018/raw_data_10x.txt", sep = "\t", header = T, check.names = F)
Human.N.data.c <- Human.N.data
Human.N.data.c$Gene <- gsub("_ENSG.*", "", Human.N.data.c$Gene)
Human.N.data.c <- data.frame(Human.N.data.c %>% group_by(Gene) %>% filter(!duplicated(Gene)), check.names = F)
rownames(Human.N.data.c) <- Human.N.data.c$Gene
Human.N.data.c <- Human.N.data.c[, -1]
Human.N.data.dsc <- Human.N.data.c[, colnames(Human.N.data.c) %in% rownames(Human.N.meta.dsc)]
Human.N.seurat <- CreateSeuratObject(counts = Human.N.data.dsc, min.cells = 3, min.features = 200, meta.data = Human.N.meta.dsc)
Human.N.seurat <- NormalizeData(Human.N.seurat, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = VariableFeatures(Human.N.seurat))
Human.N.seurat <- RunPCA(Human.N.seurat, features = VariableFeatures(Human.N.seurat), npcs = 50)
Human.N.seurat <- FindNeighbors(object = Human.N.seurat, dims = 1:30)
Human.N.seurat <- FindClusters(Human.N.seurat, resolution = 0.4)
Human.N.seurat <- RunUMAP(Human.N.seurat, dims = 1:30, seed.use = 417)

#seurat integration according to Fetus time
Human.N.seurat.list <- SplitObject(Human.N.seurat, split.by = "Fetus")
Human.N.seurat.list <- lapply(X = Human.N.seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = Human.N.seurat.list)
Human.N.anchors <- FindIntegrationAnchors(object.list = Human.N.seurat.list, anchor.features = features)
Human.N.combined <- IntegrateData(anchorset = Human.N.anchors, features.to.integrate = rownames(Human.N.seurat))
DefaultAssay(Human.N.combined) <- "integrated"

Human.N.combined <- ScaleData(Human.N.combined, verbose = FALSE)
Human.N.combined <- RunPCA(Human.N.combined, npcs = 50, verbose = FALSE)
ElbowPlot(Human.N.combined, ndims = 50)
Human.N.combined <- RunUMAP(Human.N.combined, reduction = "pca", dims = 1:30)
Human.N.combined <- FindNeighbors(Human.N.combined, reduction = "pca", dims = 1:30)
Human.N.combined <- FindClusters(Human.N.combined, resolution = 0.2)

#dataset2
library(Seurat)
P1D.DS <- read.table("./20220510-SienceAdvance-human/00.CountMatrix/P1D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P2D.DS <- read.table("./20220510-SienceAdvance-human/00.CountMatrix/P2D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P3D.DS <- read.table("./20220510-SienceAdvance-human/00.CountMatrix/P3D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P4D.DS <- read.table("./20220510-SienceAdvance-human/00.CountMatrix/P4D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P5D.DS <- read.table("./20220510-SienceAdvance-human/00.CountMatrix/P5D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P6D.10x <- Read10X("./20220510-SienceAdvance-human/00.CountMatrix/P6D-10X_filtered_feature_bc_matrix")

decidua.list <- list(P1D.DS, P2D.DS, P3D.DS, P4D.DS, P5D.DS, P6D.10x)
names(decidua.list) <- c("P1D.DS", "P2D.DS", "P3D.DS", "P4D.DS", "P5D.DS", "P6D.10x")
for (i in names(decidua.list)) {
  print(i)
  decidua.list[[i]] <- decidua.list[[i]][-grep("^RP[S,L]", rownames(decidua.list[[i]])), ]
  decidua.list[[i]] <- decidua.list[[i]][-grep("MALAT1", rownames(decidua.list[[i]])), ]
  decidua.list[[i]] <- CreateSeuratObject(decidua.list[[i]], min.cells = 3, min.features = 100, project = i)
  decidua.list[[i]]$percent.mt <- PercentageFeatureSet(decidua.list[[i]], pattern = "^MT-")
  decidua.list[[i]] <- subset(decidua.list[[i]], subset = percent.mt<25)
  decidua.list[[i]] <- NormalizeData(decidua.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  decidua.list[[i]] <- FindVariableFeatures(decidua.list[[i]], selection.method = "vst", nfeatures = 2000)
  decidua.list[[i]] <- ScaleData(decidua.list[[i]], features = rownames(decidua.list[[i]]))
}
decidua.DS <- merge(decidua.list[[1]], y = decidua.list[2:5], project = "decidua.DS")
decidua.DS$project <- "decidua.DS"
decidua.DS$sample <- "decidua"
decidua.DS$tech <- "DS"
decidua.DS <- NormalizeData(decidua.DS)
decidua.DS <- FindVariableFeatures(decidua.DS, selection.method = "vst", nfeatures = 2000)
decidua.list$P6D.10x$project <- "decidua.10X"
decidua.list$P6D.10x$sample <- "decidua"
decidua.list$P6D.10x$tech <- "10X"

features <- SelectIntegrationFeatures(object.list = list(decidua.DS, decidua.list$P6D.10x))
decidua.anchors <- FindIntegrationAnchors(object.list = list(decidua.DS, decidua.list$P6D.10x), anchor.features = features)
decidua.combined <- IntegrateData(anchorset = decidua.anchors)
DefaultAssay(decidua.combined) <- "integrated"
decidua.combined <- ScaleData(decidua.combined, verbose = FALSE)
decidua.combined <- RunPCA(decidua.combined, npcs = 30, verbose = FALSE)
decidua.combined <- RunTSNE(decidua.combined, reduction = "pca", dims = 1:16)
decidua.combined <- FindNeighbors(decidua.combined, reduction = "pca", dims = 1:16)
decidua.combined <- FindClusters(decidua.combined, resolution = 0.5)
decidua.combined <- ScaleData(decidua.combined, features = rownames(decidua.combined))
DefaultAssay(decidua.combined) <- "RNA"

FB.DSC.decidua <- subset(decidua.combined, subset = seurat_clusters %in% c(5,7,1,2,6))
FB.DSC.decidua <- RenameIdents(FB.DSC.decidua, "5" = "DSC", "7" = "DSC", "1" = "FB1", "2" = "FB1", "6" = "FB2")

#recluster, seurat CCA integration
FB.DSC.decidua <- subset(FB.DSC.decidua, subset = nFeature_RNA >1000 & percent.mt<25)
FB.DSC.decidua <- NormalizeData(FB.DSC.decidua)
FB.DSC.decidua <- FindVariableFeatures(FB.DSC.decidua, selection.method = "vst", nfeatures = 2000)
FB.DSC.decidua.list <- SplitObject(FB.DSC.decidua, split.by = "project")
features <- SelectIntegrationFeatures(object.list = FB.DSC.decidua.list)
FB.DSC.decidua.anchors <- FindIntegrationAnchors(object.list = FB.DSC.decidua.list, anchor.features = features)
FB.DSC.decidua.combined <- IntegrateData(anchorset = FB.DSC.decidua.anchors, features.to.integrate = rownames(FB.DSC.decidua))
DefaultAssay(FB.DSC.decidua.combined) <- "integrated"
FB.DSC.decidua.combined <- ScaleData(FB.DSC.decidua.combined, features = rownames(FB.DSC.decidua.combined))
FB.DSC.decidua.combined <- RunPCA(FB.DSC.decidua.combined, features = VariableFeatures(FB.DSC.decidua.combined))
ElbowPlot(FB.DSC.decidua.combined, ndims = 50)
FB.DSC.decidua.combined <- RunTSNE(FB.DSC.decidua.combined, reduction = "pca", dims = 1:20)
FB.DSC.decidua.combined <- RunUMAP(FB.DSC.decidua.combined, reduction = "pca", dims = 1:20)
FB.DSC.decidua.combined <- FindNeighbors(FB.DSC.decidua.combined, reduction = "pca", dims = 1:20)
FB.DSC.decidua.combined <- FindClusters(FB.DSC.decidua.combined, resolution = 0.4)

#downsample and integration two datasets
FB.DSC.decidua.combined$literature <- "SD"
FB.DSC.decidua.combined$primary_cluster <- paste("SD_", FB.DSC.decidua.combined$seurat_clusters)
set.seed(0)
FB.DSC.decidua.combined.sub <- subset(FB.DSC.decidua.combined, downsample = 500)
Human.N.combined$literature <- "Nature"
Human.N.combined$primary_cluster <- paste("Nature_", Human.N.combined$seurat_clusters)
set.seed(0)
Human.N.combined.sub <- subset(Human.N.combined, downsample = 500)

features <- SelectIntegrationFeatures(object.list = c(FB.DSC.decidua.combined.sub, Human.N.combined))
human.decidua.anchors <- FindIntegrationAnchors(object.list = c(FB.DSC.decidua.combined.sub, Human.N.combined), anchor.features = features)
human.decidua.combined.sub <- IntegrateData(anchorset = human.decidua.anchors, features.to.integrate = intersect(rownames(FB.DSC.decidua.combined.sub), rownames(Human.N.combined)))
DefaultAssay(human.decidua.combined.sub) <- "integrated"
human.decidua.combined.sub <- ScaleData(human.decidua.combined.sub, features = rownames(human.decidua.combined.sub))
human.decidua.combined.sub <- RunPCA(human.decidua.combined.sub, features = VariableFeatures(human.decidua.combined.sub))
ElbowPlot(human.decidua.combined.sub, ndims = 50)
human.decidua.combined.sub <- RunUMAP(human.decidua.combined.sub, reduction = "pca", dims = 1:50)
human.decidua.combined.sub <- FindNeighbors(human.decidua.combined.sub, reduction = "pca", dims = 1:50)
human.decidua.combined.sub <- FindClusters(human.decidua.combined.sub, resolution = 0.3)
DimPlot(human.decidua.combined.sub, group.by = "literature")

#Normal vs RSA
library(dplyr)
library(Seurat)
library(harmony)
library(clusterProfiler)
human.RSA1 <- Read10X("./20220930-human-disease/01.cellranger/RSA1/outs/filtered_feature_bc_matrix/")
human.RSA1 <- CreateSeuratObject(human.RSA1, min.cells = 3, min.features = 200, project = "human.RSA1")
human.RSA2 <- Read10X("./20220930-human-disease/01.cellranger/RSA2/outs/filtered_feature_bc_matrix/")
human.RSA2 <- CreateSeuratObject(human.RSA2, min.cells = 3, min.features = 200, project = "human.RSA2")
human.RSA3 <- Read10X("./20220930-human-disease/01.cellranger/RSA3/outs/filtered_feature_bc_matrix/")
human.RSA3 <- CreateSeuratObject(human.RSA3, min.cells = 3, min.features = 200, project = "human.RSA3")
human.RSA4 <- Read10X("./20220930-human-disease/01.cellranger/RSA4/outs/filtered_feature_bc_matrix/")
human.RSA4 <- CreateSeuratObject(human.RSA4, min.cells = 3, min.features = 200, project = "human.RSA4")
human.RSA5 <- Read10X("./20220930-human-disease/01.cellranger/RSA5/outs/filtered_feature_bc_matrix/")
human.RSA5 <- CreateSeuratObject(human.RSA5, min.cells = 3, min.features = 200, project = "human.RSA5")
human.RSA6 <- Read10X("./20220930-human-disease/01.cellranger/RSA6/outs/filtered_feature_bc_matrix/")
human.RSA6 <- CreateSeuratObject(human.RSA6, min.cells = 3, min.features = 200, project = "human.RSA6")
human.RSA <- merge(human.RSA1, y = list(human.RSA2, human.RSA3, human.RSA4, human.RSA5, human.RSA6), 
                   add.cell.ids = c("human.RSA1", "human.RSA2", "human.RSA3", "human.RSA4", "human.RSA5", "human.RSA6"))
human.RSA$percent.mt <- PercentageFeatureSet(human.RSA, pattern = "^MT-")
VlnPlot(human.RSA, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), pt.size = 0)
human.RSA <- subset(human.RSA, subset = percent.mt < 5 & 
                      nFeature_RNA > 500)
human.RSA <- NormalizeData(human.RSA, normalization.method = "LogNormalize", scale.factor = 10000)
human.RSA <- FindVariableFeatures(human.RSA, selection.method = "vst", nfeatures = 2000)
human.RSA <- ScaleData(human.RSA, features = rownames(human.RSA))
human.RSA <- RunPCA(human.RSA, features = VariableFeatures(human.RSA), npcs = 70)
human.RSA <- RunHarmony(human.RSA, group.by.vars = "orig.ident", plot_convergence = TRUE, max.iter.harmony = 50)
ElbowPlot(human.RSA, ndims = 70, reduction = "harmony")
human.RSA <- FindNeighbors(object = human.RSA, dims = 1:30, reduction = "harmony")
human.RSA <- RunUMAP(human.RSA, dims = 1:30, seed.use = 0912, reduction = "harmony")
human.RSA <- FindClusters(human.RSA, resolution = 0.5)
DimPlot(human.RSA, label = T)
save(human.RSA, file = "./08.human_RSA/human-RSA-seurat-harmony.RData")
FeaturePlot(human.RSA, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"))
FeaturePlot(human.RSA, features = c("HAND2", "PGR", "PTPRC", "RGS5", "PECAM1", "CDH5", "EPCAM"))

#PTPRC positive clusters
human.RSA.PTPRC <- subset(human.RSA, subset = seurat_clusters %in% c(0,1,6,11,15,17,19))
human.RSA.PTPRC <- NormalizeData(human.RSA.PTPRC, normalization.method = "LogNormalize", scale.factor = 10000)
human.RSA.PTPRC <- FindVariableFeatures(human.RSA.PTPRC, selection.method = "vst", nfeatures = 2000)
human.RSA.PTPRC <- ScaleData(human.RSA.PTPRC, features = rownames(human.RSA.PTPRC))
human.RSA.PTPRC <- RunPCA(human.RSA.PTPRC, features = VariableFeatures(human.RSA.PTPRC), npcs = 70)
human.RSA.PTPRC <- RunHarmony(human.RSA.PTPRC, group.by.vars = "orig.ident", plot_convergence = TRUE, max.iter.harmony = 50)
ElbowPlot(human.RSA.PTPRC, ndims = 70, reduction = "harmony")
human.RSA.PTPRC <- FindNeighbors(object = human.RSA.PTPRC, dims = 1:30, reduction = "harmony")
human.RSA.PTPRC <- RunUMAP(human.RSA.PTPRC, dims = 1:30, seed.use = 0912, reduction = "harmony")
human.RSA.PTPRC <- FindClusters(human.RSA.PTPRC, resolution = 0.5)
DimPlot(human.RSA.PTPRC, label = T)
FeaturePlot(human.RSA.PTPRC, features = c("PTPRC", "HAND2"))

#HAND2 positive clusters
human.RSA.HAND2 <- subset(human.RSA, subset = seurat_clusters %in% c(2,3,4,8,9,12))
human.RSA.HAND2 <- NormalizeData(human.RSA.HAND2, normalization.method = "LogNormalize", scale.factor = 10000)
human.RSA.HAND2 <- FindVariableFeatures(human.RSA.HAND2, selection.method = "vst", nfeatures = 2000)
human.RSA.HAND2 <- ScaleData(human.RSA.HAND2, features = rownames(human.RSA.HAND2))
human.RSA.HAND2 <- RunPCA(human.RSA.HAND2, features = VariableFeatures(human.RSA.HAND2), npcs = 70)
human.RSA.HAND2 <- RunHarmony(human.RSA.HAND2, group.by.vars = "orig.ident", plot_convergence = TRUE, max.iter.harmony = 50)
ElbowPlot(human.RSA.HAND2, ndims = 70, reduction = "harmony")
human.RSA.HAND2 <- FindNeighbors(object = human.RSA.HAND2, dims = 1:30, reduction = "harmony")
human.RSA.HAND2 <- RunUMAP(human.RSA.HAND2, dims = 1:30, seed.use = 0912, reduction = "harmony")
human.RSA.HAND2 <- FindClusters(human.RSA.HAND2, resolution = 0.5)
DimPlot(human.RSA.HAND2, label = T)
FeaturePlot(human.RSA.HAND2, features = c("PTPRC", "HAND2"))
VlnPlot(human.RSA.HAND2, features = c("PTPRC", "HAND2"), pt.size = 0)
save(human.RSA.HAND2, file = "./08.human_RSA/human-RSA-HAND2-20221228.RData")
human.RSA.HAND2.markers <- FindAllMarkers(human.RSA.HAND2, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
human.RSA.HAND2.markers.top5 <- human.RSA.HAND2.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(human.RSA.HAND2, features = human.RSA.HAND2.markers.top5$gene)


#Healthy samples
human.N1 <- Read10X("./20220930-human-disease/01.cellranger/Normal1/outs/filtered_feature_bc_matrix/")
human.N1 <- CreateSeuratObject(human.N1, min.cells = 3, min.features = 200, project = "human.N1")
human.N2 <- Read10X("./20220930-human-disease/01.cellranger/Normal2/outs/filtered_feature_bc_matrix/")
human.N2 <- CreateSeuratObject(human.N2, min.cells = 3, min.features = 200, project = "human.N2")
human.N3 <- Read10X("./20220930-human-disease/01.cellranger/Normal3/outs/filtered_feature_bc_matrix/")
human.N3 <- CreateSeuratObject(human.N3, min.cells = 3, min.features = 200, project = "human.N3")
human.N4 <- Read10X("./20220930-human-disease/01.cellranger/Normal4/outs/filtered_feature_bc_matrix/")
human.N4 <- CreateSeuratObject(human.N4, min.cells = 3, min.features = 200, project = "human.N4")
human.N5 <- Read10X("./20220930-human-disease/01.cellranger/Normal5/outs/filtered_feature_bc_matrix/")
human.N5 <- CreateSeuratObject(human.N5, min.cells = 3, min.features = 200, project = "human.N5")
human.N <- merge(human.N1, y = list(human.N2, human.N3, human.N4, human.N5), 
                 add.cell.ids = c("human.N1", "human.N2", "human.N3", "human.N4", "human.N5"))
human.N$percent.mt <- PercentageFeatureSet(human.N, pattern = "^MT-")
VlnPlot(human.N, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), pt.size = 0)
human.N <- subset(human.N, subset = percent.mt < 5 & 
                    nFeature_RNA > 500)
human.N <- merge(human.N1, y = list(human.N2, human.N3, human.N4, human.N5), 
                 add.cell.ids = c("human.N1", "human.N2", "human.N3", "human.N4", "human.N5"))
human.N$percent.mt <- PercentageFeatureSet(human.N, pattern = "^MT-")
VlnPlot(human.N, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), pt.size = 0)
human.N <- subset(human.N, subset = percent.mt < 5 & 
                    nFeature_RNA > 500)
human.N <- NormalizeData(human.N, normalization.method = "LogNormalize", scale.factor = 10000)
human.N <- FindVariableFeatures(human.N, selection.method = "vst", nfeatures = 2000)
human.N <- ScaleData(human.N, features = rownames(human.N))
human.N <- RunPCA(human.N, features = VariableFeatures(human.N), npcs = 70)
human.N <- RunHarmony(human.N, group.by.vars = "orig.ident", plot_convergence = TRUE, max.iter.harmony = 50)
ElbowPlot(human.N, ndims = 70, reduction = "harmony")
human.N <- FindNeighbors(object = human.N, dims = 1:30, reduction = "harmony")
human.N <- RunUMAP(human.N, dims = 1:30, seed.use = 0912, reduction = "harmony")
human.N <- FindClusters(human.N, resolution = 0.5)
DimPlot(human.N, label = T)
FeaturePlot(human.N, features = c("HAND2", "PGR", "PTPRC", "RGS5", "PECAM1", "CDH5", "EPCAM"))

#HAND2
human.N.HAND2 <- subset(human.N, subset = seurat_clusters %in% c(0,1,8,16,17,13))
human.N.HAND2 <- NormalizeData(human.N.HAND2, normalization.method = "LogNormalize", scale.factor = 10000)
human.N.HAND2 <- FindVariableFeatures(human.N.HAND2, selection.method = "vst", nfeatures = 2000)
human.N.HAND2 <- ScaleData(human.N.HAND2, features = rownames(human.N.HAND2))
human.N.HAND2 <- RunPCA(human.N.HAND2, features = VariableFeatures(human.N.HAND2), npcs = 70)
human.N.HAND2 <- RunHarmony(human.N.HAND2, group.by.vars = "orig.ident", plot_convergence = TRUE, max.iter.harmony = 50)
ElbowPlot(human.N.HAND2, ndims = 70, reduction = "harmony")
human.N.HAND2 <- FindNeighbors(object = human.N.HAND2, dims = 1:30, reduction = "harmony")
human.N.HAND2 <- RunUMAP(human.N.HAND2, dims = 1:30, seed.use = 0912, reduction = "harmony")
human.N.HAND2 <- FindClusters(human.N.HAND2, resolution = 0.5)
DimPlot(human.N.HAND2, label = T, group.by = "orig.ident")
save(human.N.HAND2, file = "./08.human_RSA/human-Normal-HAND2-20221228.RData")
FeaturePlot(human.N.HAND2, features = c("PTPRC", "HAND2"))
human.N.HAND2.markers <- FindAllMarkers(human.N.HAND2, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
human.N.HAND2.markers.top5 <- human.N.HAND2.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(human.N.HAND2, features = human.N.HAND2.markers.top5$gene)

#PTPRC
human.N.PTPRC <- subset(human.N, subset = seurat_clusters %in% c(2,4,9,3,7,14,12))
human.N.PTPRC <- NormalizeData(human.N.PTPRC, normalization.method = "LogNormalize", scale.factor = 10000)
human.N.PTPRC <- FindVariableFeatures(human.N.PTPRC, selection.method = "vst", nfeatures = 2000)
human.N.PTPRC <- ScaleData(human.N.PTPRC, features = rownames(human.N.PTPRC))
human.N.PTPRC <- RunPCA(human.N.PTPRC, features = VariableFeatures(human.N.PTPRC), npcs = 70)
human.N.PTPRC <- RunHarmony(human.N.PTPRC, group.by.vars = "orig.ident", plot_convergence = TRUE, max.iter.harmony = 50)
ElbowPlot(human.N.PTPRC, ndims = 70, reduction = "harmony")
human.N.PTPRC <- FindNeighbors(object = human.N.PTPRC, dims = 1:30, reduction = "harmony")
human.N.PTPRC <- RunUMAP(human.N.PTPRC, dims = 1:30, seed.use = 0912, reduction = "harmony")
human.N.PTPRC <- FindClusters(human.N.PTPRC, resolution = 0.5)
DimPlot(human.N.PTPRC, label = T)
FeaturePlot(human.N.PTPRC, features = c("PTPRC", "HAND2"))


#merge RSA and Normal
human.RSA$sample <- "RSA"
human.N$sample <- "Normal"
human.RSA.N <- merge(human.RSA, human.N)
human.RSA.N <- NormalizeData(human.RSA.N, normalization.method = "LogNormalize", scale.factor = 10000)
human.RSA.N <- FindVariableFeatures(human.RSA.N, selection.method = "vst", nfeatures = 2000)
human.RSA.N <- ScaleData(human.RSA.N, features = rownames(human.RSA.N))
human.RSA.N <- RunPCA(human.RSA.N, features = VariableFeatures(human.RSA.N), npcs = 70)
human.RSA.N <- RunHarmony(human.RSA.N, group.by.vars = "orig.ident", plot_convergence = TRUE, max.iter.harmony = 50)
ElbowPlot(human.RSA.N, ndims = 70, reduction = "harmony")
human.RSA.N <- FindNeighbors(object = human.RSA.N, dims = 1:30, reduction = "harmony")
human.RSA.N <- RunUMAP(human.RSA.N, dims = 1:30, seed.use = 0912, reduction = "harmony")
human.RSA.N <- FindClusters(human.RSA.N, resolution = 0.5)
DimPlot(human.RSA.N, label = T, group.by = "orig.ident")
FeaturePlot(human.RSA.N, features = c("PTPRC", "HAND2"))

#subset iDSC
human.N.HAND2.iDSC <- subset(human.N.HAND2, subset = seurat_clusters %in% c(8,9))
human.RSA.HAND2.iDSC <- subset(human.RSA.HAND2, subset = seurat_clusters %in% c(9,10))
human.N.HAND2.iDSC$cluster <- "Normal_iDSC"
human.RSA.HAND2.iDSC$cluster <- "RSA_iDSC"
human.Normal.RSA.iDSC <- merge(human.N.HAND2.iDSC, human.RSA.HAND2.iDSC)
human.Normal.RSA.iDSC <- NormalizeData(human.Normal.RSA.iDSC)
human.Normal.RSA.iDSC <- ScaleData(human.Normal.RSA.iDSC, features = rownames(human.Normal.RSA.iDSC))
human.Normal.RSA.iDSC <- FindVariableFeatures(human.Normal.RSA.iDSC)
human.Normal.RSA.iDSC <- RunPCA(human.Normal.RSA.iDSC, features = VariableFeatures(human.Normal.RSA.iDSC), npcs = 70)
human.Normal.RSA.iDSC <- RunHarmony(human.Normal.RSA.iDSC, group.by.vars = "orig.ident", plot_convergence = TRUE, max.iter.harmony = 50)
ElbowPlot(human.Normal.RSA.iDSC, ndims = 70, reduction = "harmony")
human.Normal.RSA.iDSC <- FindNeighbors(object = human.Normal.RSA.iDSC, dims = 1:20, reduction = "harmony")
human.Normal.RSA.iDSC <- RunUMAP(human.Normal.RSA.iDSC, dims = 1:20, seed.use = 0912, reduction = "harmony")
human.Normal.RSA.iDSC <- FindClusters(human.Normal.RSA.iDSC, res = 0.4)
DimPlot(human.Normal.RSA.iDSC, label = T)
DimPlot(human.Normal.RSA.iDSC, label = T, group.by = "cluster")
human.Normal.RSA.iDSC@active.ident <- factor(human.Normal.RSA.iDSC$cluster)
human.Normal.RSA.iDSC.markers <- FindAllMarkers(human.Normal.RSA.iDSC, only.pos = T, min.pct = 0.75, logfc.threshold = 0.25)
human.Normal.RSA.iDSC.markers.top10 <- human.Normal.RSA.iDSC.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(human.Normal.RSA.iDSC, features = human.Normal.RSA.iDSC.markers.top10$gene)
human.Normal.RSA.iDSC.avg <- AverageExpression(human.Normal.RSA.iDSC, slot = "data", assays = "RNA")
human.Normal.RSA.iDSC.avg <- human.Normal.RSA.iDSC.avg$RNA
human.Normal.RSA.iDSC.avg$change <- "NO"
human.Normal.RSA.iDSC.avg[rownames(human.Normal.RSA.iDSC.avg) %in% human.Normal.RSA.iDSC.markers[human.Normal.RSA.iDSC.markers$cluster == "Normal_iDSC", ]$gene, ]$change <- "Normal_iDSC"
human.Normal.RSA.iDSC.avg[rownames(human.Normal.RSA.iDSC.avg) %in% human.Normal.RSA.iDSC.markers[human.Normal.RSA.iDSC.markers$cluster == "RSA_iDSC", ]$gene, ]$change <- "RSA_iDSC"
human.Normal.RSA.iDSC.avg$change <- factor(human.Normal.RSA.iDSC.avg$change, levels = c("RSA_iDSC", "NO", "Normal_iDSC"))
human.Normal.RSA.iDSC.avg.anno <- human.Normal.RSA.iDSC.avg[human.Normal.RSA.iDSC.avg$change != "NO", ]
human.Normal.RSA.iDSC.avg.anno$gene <- rownames(human.Normal.RSA.iDSC.avg.anno)












