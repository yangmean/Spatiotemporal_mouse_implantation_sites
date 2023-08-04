library(dplyr)
library(Seurat)
library(harmony)
human.S.meta <- read.table("../../20220508-download-data/Nature-2018/meta_10x.txt", sep = "\t", header = T, check.names = F)
human.S.meta.dsc <- human.S.meta[human.S.meta$final_cluster %in% c(0, 3, 21), ]
human.S.data <- read.table("../../20220508-download-data/Nature-2018/raw_data_10x.txt", sep = "\t", header = T, check.names = F)
human.S.data.c <- human.S.data
human.S.data.c$Gene <- gsub("_ENSG.*", "", human.S.data.c$Gene)
human.S.data.c <- data.frame(human.S.data.c %>% group_by(Gene) %>% filter(!duplicated(Gene)), check.names = F)
rownames(human.S.data.c) <- human.S.data.c$Gene
human.S.data.c <- human.S.data.c[, -1]
human.S.data.dsc <- human.S.data.c[, colnames(human.S.data.c) %in% rownames(human.S.meta.dsc)]
human.S.seurat <- CreateSeuratObject(counts = human.S.data.dsc, min.cells = 3, min.features = 200, meta.data = human.S.meta.dsc)
save(human.S.seurat, file = "../../20220508-download-data/Nature-2018/Human-nature-2018-DSC-20220509.RData")
human.S.seurat <- NormalizeData(human.S.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
human.S.seurat <- FindVariableFeatures(human.S.seurat, selection.method = "vst", nfeatures = 2000)
#human.S.seurat$CC.Score <- human.S.seurat$S.Score - human.S.seurat$G2M.Score
human.S.seurat <- ScaleData(human.S.seurat, features = VariableFeatures(human.S.seurat))
human.S.seurat <- RunPCA(human.S.seurat, features = VariableFeatures(human.S.seurat), npcs = 50)
ElbowPlot(human.S.seurat, ndims = 50)
human.S.seurat <- FindNeighbors(object = human.S.seurat, dims = 1:30)
human.S.seurat <- FindClusters(human.S.seurat, resolution = 0.4)
human.S.seurat <- RunUMAP(human.S.seurat, dims = 1:30, seed.use = 417)
human.S.seurat.markers <- FindAllMarkers(human.S.seurat, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
human.S.seurat.markers.top10 <- human.S.seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(human.S.seurat, features = human.S.seurat.markers.top10$gene)
FeaturePlot(human.S.seurat, features = "PTPRC")
DimPlot(human.S.seurat, group.by = "annotation")
DimPlot(human.S.seurat, group.by = "location")
DimPlot(human.S.seurat, group.by = "Fetus")
DimPlot(human.S.seurat, label = T)

#seurat integration according to Fetus time
human.S.seurat.list <- SplitObject(human.S.seurat, split.by = "Fetus")
human.S.seurat.list <- lapply(X = human.S.seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = human.S.seurat.list)
human.S.anchors <- FindIntegrationAnchors(object.list = human.S.seurat.list, anchor.features = features)
human.S.combined <- IntegrateData(anchorset = human.S.anchors, features.to.integrate = rownames(human.S.seurat))
DefaultAssay(human.S.combined) <- "integrated"

human.S.combined <- ScaleData(human.S.combined, verbose = FALSE)
human.S.combined <- RunPCA(human.S.combined, npcs = 50, verbose = FALSE)
ElbowPlot(human.S.combined, ndims = 50)
human.S.combined <- RunUMAP(human.S.combined, reduction = "pca", dims = 1:30)
human.S.combined <- FindNeighbors(human.S.combined, reduction = "pca", dims = 1:30)
human.S.combined <- FindClusters(human.S.combined, resolution = 0.2)
clustree::clustree(human.S.combined, prefix = "integrated_snn_res.")
DimPlot(human.S.combined, label = T, group.by = "integrated_snn_res.0.2")
FeaturePlot(human.S.combined, features = c("nFeature_RNA", "nCount_RNA"))
human.S.combined.markers <- FindAllMarkers(human.S.combined, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
human.S.combined.markers.top10 <- human.S.combined.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(human.S.combined, features = human.S.combined.markers.top10$gene)
DimPlot(human.S.combined, label = T)
DimPlot(human.S.combined, group.by = "Fetus")
DimPlot(human.S.combined, group.by = "annotation")
FeaturePlot(human.S.combined, features = c("CCL3", "PTPRC", "ACTA2"))
human.S.combined$cluster1 <- paste("N_CCA_", human.S.combined@active.ident, sep = "")
human.S.combined$orig.cellname <- human.S.combined$annotation
human.S.combined$literature <- "Nature"
save(human.S.combined, file = "20220622-Comments/02.human-data/nature/Nature-human-DSC-recluster-CCA.RData")
write.table(human.S.combined.markers, file = "20220622-Comments/02.human-data/nature/Nature-human-DSC-CCA-res0.2-markers.csv", sep = ",", quote = F)
#GO
library("clusterProfiler")
enrich.GO.function <- function(geneset){
  enrich.GO <- enrichGO(gene = geneset,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "SYMBOL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
  #enrich.GO <- simplify(enrich.GO, cutoff = 0.7, by = "p.adjust", select_fun = min)
  enrich.GO <- enrich.GO@result[enrich.GO@result$pvalue<0.01, ]
  return(enrich.GO)
}
human.DSC.GO.df <- data.frame()
for (c in levels(human.S.combined.markers$cluster)) {
  print(c)
  geneset <- human.S.combined.markers[human.S.combined.markers$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    gene = enrich.GO$geneID)
  human.DSC.GO.df <- rbind(human.DSC.GO.df, tmp)
}

#harmony integration according to Fetus time
human.S.seurat <- CreateSeuratObject(counts = human.S.data.dsc, min.cells = 3, min.features = 200, meta.data = human.S.meta.dsc)
human.S.seurat <- NormalizeData(human.S.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
human.S.seurat <- FindVariableFeatures(human.S.seurat, selection.method = "vst", nfeatures = 2000)
human.S.seurat <- ScaleData(human.S.seurat, features = VariableFeatures(human.S.seurat))
human.S.seurat <- RunPCA(human.S.seurat, npcs = 50, verbose = FALSE)
human.S.harmony <- subset(human.S.seurat, subset = orig.ident %in% c("FCA7167221", "FCA7167222", "FCA7167224", "FCA7167226", "FCA7196218", "FCA7196219", 
                                                                     "FCA7196224", "FCA7196225", "FCA7474062", "FCA7474063", "FCA7511881", "FCA7511882"))
human.S.harmony <- RunHarmony(human.S.harmony, group.by.vars = "Fetus")
ElbowPlot(human.S.harmony, reduction = "harmony", ndims = 50)
human.S.harmony <- RunUMAP(human.S.harmony, reduction = "harmony", dims = 1:30)
human.S.harmony <- FindNeighbors(human.S.harmony, reduction = "harmony", dims = 1:30)
human.S.harmony <- FindClusters(human.S.harmony, resolution = 0.5)
clustree::clustree(human.S.harmony, prefix = "RNA_snn_res.")

DimPlot(human.S.harmony, label = T)
DimPlot(human.S.harmony, label = T, group.by = "Fetus")
DimPlot(human.S.harmony, label = T, group.by = "annotation")
FeaturePlot(human.S.harmony, features = c("CCL3", "PTPRC", "ACTA2"))
FeaturePlot(human.S.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
human.S.harmony.markers <- FindAllMarkers(human.S.harmony, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
human.S.harmony.markers.top10 <- human.S.harmony.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(human.S.harmony, features = human.S.harmony.markers.top10$gene)

pdf("../../20220508-download-data/Nature-2018/dS-recluster-UMAP.pdf", width = 10, height = 10)
DimPlot(human.S.harmony, label = T, pt.size = 0.5)
dev.off()
DimPlot(human.S.harmony, label = T, group.by = "Fetus")
# pdf("../../20220508-download-data/Nature-2018/dS-3clusters-UMAP.pdf", width = 10, height = 10)
DimPlot(human.S.harmony, label = F, group.by = "annotation", pt.size = 1)
dev.off()
pdf("../../20220508-download-data/Nature-2018/dS-featureplot-dS-immune-markers.pdf", width = 15, height = 15)
FeaturePlot(human.S.harmony, features = c("ACTA2", "IGFBP1", "CCL3", "CCL4", "NKG7", "CD74", "GZMB"), ncol = 3)
dev.off()

human.S.harmony.markers <- FindAllMarkers(human.S.harmony, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
human.S.harmony.markers.top10 <- human.S.harmony.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
pdf("../../20220508-download-data/Nature-2018/dS-clusters-heatmap-top10.pdf", width = 10, height = 10)
DoHeatmap(human.S.harmony, features = human.S.harmony.markers.top10$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()
save(human.S.harmony, file = "../../20220508-download-data/Nature-2018/Human-nature-2018-DSC-harmony-20220510.RData")


#raw data analysis from cellranger result
#read tables
decidua.CD45neg.D6 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45-_D6/")
decidua.CD45neg.D7.1 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45-_D7_1/")
decidua.CD45neg.D7.2 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45-_D7_2/")
decidua.CD45neg.D8 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45-_D8/")
decidua.CD45neg.D9 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45-_D9/")
decidua.CD45neg.D10 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45-_D10/")
decidua.CD45neg.D12 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45-_D12/")

decidua.CD45pos.D6.1 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45pos_D6_1/")
decidua.CD45pos.D6.2 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45pos_D6_2/")
decidua.CD45pos.D7 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45pos_D7/")
decidua.CD45pos.D8 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45pos_D8/")
decidua.CD45pos.D9 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45pos_D9/")
decidua.CD45pos.D10 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45pos_D10/")
decidua.CD45pos.D12 <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-Nature-human/decidua_CD45pos_D12/")

N.decidua.list <- list(decidua.CD45neg.D6, decidua.CD45neg.D7.1, decidua.CD45neg.D7.2, decidua.CD45neg.D8, decidua.CD45neg.D9, decidua.CD45neg.D10, decidua.CD45neg.D12, 
                       decidua.CD45pos.D6.1, decidua.CD45pos.D6.2, decidua.CD45pos.D7, decidua.CD45pos.D8, decidua.CD45pos.D9, decidua.CD45pos.D10, decidua.CD45pos.D12)
names(N.decidua.list) <- c("decidua.CD45neg.D6", "decidua.CD45neg.D7.1", "decidua.CD45neg.D7.2", "decidua.CD45neg.D8", "decidua.CD45neg.D9", "decidua.CD45neg.D10", "decidua.CD45neg.D12", 
                         "decidua.CD45pos.D6.1", "decidua.CD45pos.D6.2", "decidua.CD45pos.D7", "decidua.CD45pos.D8", "decidua.CD45pos.D9", "decidua.CD45pos.D10", "decidua.CD45pos.D12")
for (i in names(N.decidua.list)) {
  print(i)
  N.decidua.list[[i]] <- CreateSeuratObject(N.decidua.list[[i]], min.cells = 3, min.features = 500, project = i)
  N.decidua.list[[i]]$percent.mt <- PercentageFeatureSet(N.decidua.list[[i]], pattern = "^MT-")
  N.decidua.list[[i]] <- subset(N.decidua.list[[i]], subset = percent.mt<20)
  N.decidua.list[[i]] <- NormalizeData(N.decidua.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  N.decidua.list[[i]] <- FindVariableFeatures(N.decidua.list[[i]], selection.method = "vst", nfeatures = 2000)
  N.decidua.list[[i]] <- ScaleData(N.decidua.list[[i]], features = rownames(N.decidua.list[[i]]))
}
N.decidua.merge <- merge(N.decidua.list[[1]], y = N.decidua.list[2:14], project = "N.decidua.merge")
N.decidua.merge <- NormalizeData(N.decidua.merge, normalization.method = "LogNormalize", scale.factor = 10000)
N.decidua.merge <- FindVariableFeatures(N.decidua.merge, selection.method = "vst", nfeatures = 2000)
N.decidua.merge <- ScaleData(N.decidua.merge, features = rownames(N.decidua.merge))
N.decidua.merge <- RunPCA(N.decidua.merge, npcs = 50)
ElbowPlot(N.decidua.merge, ndims = 50)
N.decidua.merge <- FindNeighbors(N.decidua.merge, reduction = "pca", dims = 1:30)
N.decidua.merge <- RunUMAP(N.decidua.merge, reduction = "pca", dims = 1:30)
N.decidua.merge <- FindClusters(N.decidua.merge, resolution = 0.5)
DimPlot(N.decidua.merge, label = T)
FeaturePlot(N.decidua.merge, features = c("ACTA2", "PTPRC", "RGS5"))
N.decidua.DSC <- subset(N.decidua.merge, subset = seurat_clusters %in% c(0, 5, 6, 10, 13))
N.decidua.DSC <- NormalizeData(N.decidua.DSC, normalization.method = "LogNormalize", scale.factor = 10000)
N.decidua.DSC <- FindVariableFeatures(N.decidua.DSC, selection.method = "vst", nfeatures = 2000)
N.decidua.DSC <- ScaleData(N.decidua.DSC, features = rownames(N.decidua.DSC))
N.decidua.DSC <- RunPCA(N.decidua.DSC, npcs = 50)
ElbowPlot(N.decidua.DSC, ndims = 50)
N.decidua.DSC <- FindNeighbors(N.decidua.DSC, reduction = "pca", dims = 1:30)
N.decidua.DSC <- RunUMAP(N.decidua.DSC, reduction = "pca", dims = 1:30)
N.decidua.DSC <- FindClusters(N.decidua.DSC, resolution = 0.5)
DimPlot(N.decidua.DSC, label = T)
FeaturePlot(N.decidua.DSC, features = c("ACTA2", "CCL4", "NKG7", "PTPRC"))
DimPlot(N.decidua.DSC, group.by = "orig.ident", label = T)


