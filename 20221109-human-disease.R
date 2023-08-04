library(dplyr)
library(Seurat)
library(harmony)
library(clusterProfiler)
human.RSA1 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/RSA1/outs/filtered_feature_bc_matrix/")
human.RSA1 <- CreateSeuratObject(human.RSA1, min.cells = 3, min.features = 200, project = "human.RSA1")
human.RSA2 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/RSA2/outs/filtered_feature_bc_matrix/")
human.RSA2 <- CreateSeuratObject(human.RSA2, min.cells = 3, min.features = 200, project = "human.RSA2")
human.RSA3 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/RSA3/outs/filtered_feature_bc_matrix/")
human.RSA3 <- CreateSeuratObject(human.RSA3, min.cells = 3, min.features = 200, project = "human.RSA3")
human.RSA4 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/RSA4/outs/filtered_feature_bc_matrix/")
human.RSA4 <- CreateSeuratObject(human.RSA4, min.cells = 3, min.features = 200, project = "human.RSA4")
human.RSA5 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/RSA5/outs/filtered_feature_bc_matrix/")
human.RSA5 <- CreateSeuratObject(human.RSA5, min.cells = 3, min.features = 200, project = "human.RSA5")
human.RSA6 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/RSA6/outs/filtered_feature_bc_matrix/")
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
save(human.RSA, file = "20220622-Comments/08.human_RSA/human-RSA-seurat-harmony.RData")
FeaturePlot(human.RSA, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"))
FeaturePlot(human.RSA, features = c("HAND2", "PGR", "PTPRC", "RGS5", "PECAM1", "CDH5", "EPCAM"))

#PTPRC
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

#HAND2
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
save(human.RSA.HAND2, file = "20220622-Comments/08.human_RSA/human-RSA-HAND2-20221228.RData")
human.RSA.HAND2.markers <- FindAllMarkers(human.RSA.HAND2, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
human.RSA.HAND2.markers.top5 <- human.RSA.HAND2.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(human.RSA.HAND2, features = human.RSA.HAND2.markers.top5$gene)


#Healthy samples
human.N1 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/Normal1/outs/filtered_feature_bc_matrix/")
human.N1 <- CreateSeuratObject(human.N1, min.cells = 3, min.features = 200, project = "human.N1")
human.N2 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/Normal2/outs/filtered_feature_bc_matrix/")
human.N2 <- CreateSeuratObject(human.N2, min.cells = 3, min.features = 200, project = "human.N2")
human.N3 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/Normal3/outs/filtered_feature_bc_matrix/")
human.N3 <- CreateSeuratObject(human.N3, min.cells = 3, min.features = 200, project = "human.N3")
human.N4 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/Normal4/outs/filtered_feature_bc_matrix/")
human.N4 <- CreateSeuratObject(human.N4, min.cells = 3, min.features = 200, project = "human.N4")
human.N5 <- Read10X("/Data3/yangmin/Jennie/20220930-human-disease/01.cellranger/Normal5/outs/filtered_feature_bc_matrix/")
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
save(human.N.HAND2, file = "20220622-Comments/08.human_RSA/human-Normal-HAND2-20221228.RData")
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
pdf("20220622-Comments/05.DBA/00.figures/46.Human-Normal-RSA-iDSC-all-compare-dotplot.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(human.Normal.RSA.iDSC.avg, mapping = aes(x = log2(Normal_iDSC + 1), y = log2(RSA_iDSC + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  theme_classic()
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/46.Human-Normal-RSA-iDSC-all-compare-dotplot-annotation.pdf", width = 10, height = 9.5)
ggplot() + 
  geom_point(human.Normal.RSA.iDSC.avg, mapping = aes(x = log2(Normal_iDSC + 1), y = log2(RSA_iDSC + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  annotate(geom = "text", x = log2(human.Normal.RSA.iDSC.avg.anno$Normal_iDSC+1), y = log2(human.Normal.RSA.iDSC.avg.anno$RSA_iDSC+1), 
           label = human.Normal.RSA.iDSC.avg.anno$gene, size = 2) + 
  theme_classic()
dev.off()



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
human.iDSC.GO.df <- data.frame()
for (c in unique(human.Normal.RSA.iDSC.markers$cluster)) {
  print(c)
  geneset <- human.Normal.RSA.iDSC.markers[human.Normal.RSA.iDSC.markers$cluster == c, ]$gene
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
human.iDSC.GO.df.top5 <- human.iDSC.GO.df %>% group_by(cluster) %>% top_n(n = 10, -log(pvalue))
pdf("20220622-Comments/05.DBA/00.figures/47.Human-Normal-RSA-iDSC-all-compare-DEGs-GO.pdf")
ggplot() + 
  geom_bar(human.iDSC.GO.df.top5, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()


