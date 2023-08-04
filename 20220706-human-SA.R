library(Seurat)
P1D.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P1D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P1V.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P1V-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P2D.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P2D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P2V.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P2V-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P3D.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P3D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P3V.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P3V-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P4D.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P4D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P4V.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P4V-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P5D.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P5D-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P5V.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P5V-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)
P6V.DS <- read.table("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P6V-DS_out_gene_exon_tagged_dge.txt", sep = "\t", header = T, check.names = F, row.names = 1)

P6D.10x <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P6D-10X_filtered_feature_bc_matrix")
P7V.10x <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P7V-10X_filtered_feature_bc_matrix/")
P8V.10x <- Read10X("/Data1/Users/yangmin/Projects/Jennie/20220510-SienceAdvance-human/00.CountMatrix/P8V-10X_filtered_feature_bc_matrix/")

#analyze decidua samples
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
#integration
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
decidua.combined.markers <- FindAllMarkers(decidua.combined, min.pct = 0.5, logfc.threshold = 0.75, only.pos = T)
decidua.combined.markers.top10 <- decidua.combined.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DimPlot(decidua.combined, label = T)
DoHeatmap(decidua.combined, features = decidua.combined.markers.top10$gene)
DSC.markers <- c("PRL", "IGFBP1", "APOA1", "CHI3L2", "SERPINA3", "IL1B", "PROK1", "DCN")
SMC.markers <- c("MYH11", "ACTA2", "RGS5", "PI15", "NDUFA4L2")
FB1.markers <- c("COL1A1", "COL1A2", "COL3A1", "DCN")
FB2.markers <- c("ARC", "GEM", "BDKRB1", "PLIN2")
DefaultAssay(decidua.combined) <- "RNA"
FeaturePlot(decidua.combined, features = c(DSC.markers, SMC.markers, FB1.markers, FB2.markers), ncol = 6)
save(decidua.combined, file = "20220622-Comments/02.human-data/science_advance/decidua_commbined.RData")

#analyze placenta samples
placenta.list <- list(P1V.DS, P2V.DS, P3V.DS, P4V.DS, P5V.DS, P6V.DS, P7V.10x, P8V.10x)
names(placenta.list) <- c("P1V.DS", "P2V.DS", "P3V.DS", "P4V.DS", "P5V.DS", "P6V.DS", "P7V.10x", "P8V.10x")
for (i in names(placenta.list)) {
  print(i)
  placenta.list[[i]] <- placenta.list[[i]][-grep("^RP[S,L]", rownames(placenta.list[[i]])), ]
  placenta.list[[i]] <- placenta.list[[i]][-grep("MALAT1", rownames(placenta.list[[i]])), ]
  placenta.list[[i]] <- CreateSeuratObject(placenta.list[[i]], min.cells = 3, min.features = 100, project = i)
  placenta.list[[i]]$percent.mt <- PercentageFeatureSet(placenta.list[[i]], pattern = "^MT-")
  placenta.list[[i]] <- subset(placenta.list[[i]], subset = percent.mt<25)
  placenta.list[[i]] <- NormalizeData(placenta.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  placenta.list[[i]] <- FindVariableFeatures(placenta.list[[i]], selection.method = "vst", nfeatures = 2000)
  placenta.list[[i]] <- ScaleData(placenta.list[[i]], features = rownames(placenta.list[[i]]))
}
placenta.DS <- merge(placenta.list[[1]], y = placenta.list[2:6], project = "placenta.DS")
placenta.DS$project <- "placenta.DS"
placenta.DS$sample <- "placenta"
placenta.DS$tech <- "DS"
placenta.DS <- NormalizeData(placenta.DS)
placenta.DS <- FindVariableFeatures(placenta.DS, selection.method = "vst", nfeatures = 2000)
placenta.10X <- merge(placenta.list[[1]], y = placenta.list[[2]], project = "placenta.10X")
placenta.10X$project <- "placenta.10X"
placenta.10X$sample <- "placenta"
placenta.10X$tech <- "10X"
placenta.10X <- NormalizeData(placenta.10X)
placenta.10X <- FindVariableFeatures(placenta.10X, selection.method = "vst", nfeatures = 2000)

#integration
features <- SelectIntegrationFeatures(object.list = list(placenta.DS, placenta.10X))
placenta.anchors <- FindIntegrationAnchors(object.list = list(placenta.DS, placenta.10X), anchor.features = features)
placenta.combined <- IntegrateData(anchorset = placenta.anchors)
DefaultAssay(placenta.combined) <- "integrated"
placenta.combined <- ScaleData(placenta.combined, verbose = FALSE)
placenta.combined <- RunPCA(placenta.combined, npcs = 30, verbose = FALSE)
placenta.combined <- RunTSNE(placenta.combined, reduction = "pca", dims = 1:16)
placenta.combined <- FindNeighbors(placenta.combined, reduction = "pca", dims = 1:16)
placenta.combined <- FindClusters(placenta.combined, resolution = 1)
placenta.combined <- ScaleData(placenta.combined, features = rownames(placenta.combined))
placenta.combined.markers <- FindAllMarkers(placenta.combined, min.pct = 0.5, logfc.threshold = 0.75, only.pos = T)
placenta.combined.markers.top10 <- placenta.combined.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DimPlot(placenta.combined, label = T)
DoHeatmap(placenta.combined, features = placenta.combined.markers.top10$gene)
placenta.combined <- subset(placenta.combined, subset = seurat_clusters != 0)#再重新使用上面的将维分析,start with DefaultAssay function

FB13.markers <- c("COL1A1", "COL1A2", "COL3A1", "DLK1", "EGFL6", "ACTA2", "VIM", "DES", "MFAP4", "OGN", "S100A4")
FB2.markers <- c("REN", "AGTR1", "IGFBP7", "AREG")
DefaultAssay(placenta.combined) <- "RNA"
FeaturePlot(placenta.combined, features = c(FB13.markers, FB2.markers), ncol = 4)
save(placenta.combined, file = "20220622-Comments/02.human-data/science_advance/placenta_commbined.RData")

#subset FB/DSC from decidua
FB.DSC.decidua <- subset(decidua.combined, subset = seurat_clusters %in% c(5,7,1,2,6))
FB.DSC.decidua <- RenameIdents(FB.DSC.decidua, "5" = "DSC", "7" = "DSC", "1" = "FB1", "2" = "FB1", "6" = "FB2")
FB.DSC.decidua$orig.cellname <- FB.DSC.decidua@active.ident
FeaturePlot(FB.DSC.decidua, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
DimPlot(FB.DSC.decidua, label = T)
p1 <- DimPlot(FB.DSC.decidua, group.by = c("project"))
p2 <- DimPlot(FB.DSC.decidua, group.by = c("sample"))
p3 <- DimPlot(FB.DSC.decidua, group.by = c("tech"))
p1 + p2 + p3
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
DimPlot(FB.DSC.decidua.combined, label = T)
DimPlot(FB.DSC.decidua.combined, group.by = c("orig.ident"))
clustree::clustree(FB.DSC.decidua.combined, prefix = "integrated_snn_res.")
FB.DSC.decidua.combined.markers <- FindAllMarkers(FB.DSC.decidua.combined, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
FB.DSC.decidua.combined.markers.top10 <- FB.DSC.decidua.combined.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(FB.DSC.decidua.combined, features = FB.DSC.decidua.combined.markers.top10$gene)

p1 <- DimPlot(FB.DSC.decidua.combined, group.by = c("project"))
p2 <- DimPlot(FB.DSC.decidua.combined, group.by = c("sample"))
p3 <- DimPlot(FB.DSC.decidua.combined, group.by = c("tech"))
p1 + p2 + p3
FeaturePlot(FB.DSC.decidua.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FB.DSC.decidua.combined <- ScaleData(FB.DSC.decidua.combined, features = rownames(FB.DSC.decidua.combined))
FB.DSC.decidua.combined.markers <- FindAllMarkers(FB.DSC.decidua.combined, min.pct = 0.5, logfc.threshold = 0.25, only.pos = T)
FB.DSC.decidua.combined.markers.top10 <- FB.DSC.decidua.combined.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(FB.DSC.decidua.combined, features = FB.DSC.decidua.combined.markers.top10$gene)
DefaultAssay(FB.DSC.decidua.combined) <- "RNA"
FeaturePlot(FB.DSC.decidua.combined, features = c("CCL4", "NKG7", "PTPRC"), ncol = 3)
FB.DSC.decidua.combined$cluster1 <- paste("SD_CCA_", FB.DSC.decidua.combined@active.ident, sep = "")
save(FB.DSC.decidua.combined, file = "20220622-Comments/02.human-data/science_advance/decidua-subset-FB-CCA.RData")
write.table(FB.DSC.decidua.combined.markers, file = "20220622-Comments/02.human-data/science_advance/decidua-SD-CCA-res0.3-markers.csv", sep = ",", quote = F)
#harmony integration
FB.DSC.decidua <- subset(decidua.combined, subset = seurat_clusters %in% c(5,7,1,2,6))
FB.DSC.decidua <- RenameIdents(FB.DSC.decidua, "5" = "DSC", "7" = "DSC", "1" = "FB1", "2" = "FB1", "6" = "FB2")
FB.DSC.decidua$orig.cellname <- FB.DSC.decidua@active.ident
DefaultAssay(FB.DSC.decidua) <- "RNA"
FB.DSC.decidua <- NormalizeData(FB.DSC.decidua)
FB.DSC.decidua <- FindVariableFeatures(FB.DSC.decidua, selection.method = "vst", nfeatures = 2000)
FB.DSC.decidua <- ScaleData(FB.DSC.decidua, features = rownames(FB.DSC.decidua))
FB.DSC.decidua <- RunPCA(FB.DSC.decidua, features = VariableFeatures(FB.DSC.decidua))
FB.DSC.decidua.harmony <- RunHarmony(FB.DSC.decidua, group.by.vars = "project")
ElbowPlot(FB.DSC.decidua.harmony, ndims = 50, reduction = "harmony")
FB.DSC.decidua.harmony <- RunUMAP(FB.DSC.decidua.harmony, reduction = "harmony", dims = 1:30)
FB.DSC.decidua.harmony <- FindNeighbors(FB.DSC.decidua.harmony, reduction = "harmony", dims = 1:30)
FB.DSC.decidua.harmony <- FindClusters(FB.DSC.decidua.harmony, resolution = 0.5)
DimPlot(FB.DSC.decidua.harmony, label = T)
p1 <- DimPlot(FB.DSC.decidua.harmony, group.by = c("project"))
p2 <- DimPlot(FB.DSC.decidua.harmony, group.by = c("sample"))
p3 <- DimPlot(FB.DSC.decidua.harmony, group.by = c("tech"))
p4 <- DimPlot(FB.DSC.decidua.harmony, group.by = "orig.ident")
DimPlot(FB.DSC.decidua.harmony, group.by = "orig.ident", split.by = "orig.ident")
p1 + p2 + p3 + p4
FeaturePlot(FB.DSC.decidua.harmony, features = c("CCL4", "NKG7", "PTPRC", "ACTA2"), ncol = 2)
FB.DSC.decidua.harmony.markers <- FindAllMarkers(FB.DSC.decidua.harmony, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
FB.DSC.decidua.harmony.markers.top10 <- FB.DSC.decidua.harmony.markers %>% group_by(cluster) %>% top_n(n = 9, avg_logFC)
DoHeatmap(FB.DSC.decidua.harmony, features = FB.DSC.decidua.harmony.markers.top10$gene)
save(FB.DSC.decidua.harmony, file = "20220622-Comments/02.human-data/science_advance/decidua-subset-FB-harmony.RData")

#subset FB/DSC from placenta
FB.DSC.placenta <- subset(placenta.combined, subset = seurat_clusters %in% c(3,9,17,5,14))
FB.DSC.placenta <- RenameIdents(FB.DSC.placenta, "3" = "Pla_FB1", "9" = "Pla_FB1", "17" = "Pla_FB1", "5" = "Pla_FB2", "14" = "Pla_FB3")
FB.DSC.placenta$orig.cellname <- FB.DSC.placenta@active.ident
FeaturePlot(FB.DSC.placenta, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
DimPlot(FB.DSC.placenta, label = T)
#recluster, seurat CCA integration
FB.DSC.placenta <- NormalizeData(FB.DSC.placenta)
FB.DSC.placenta <- FindVariableFeatures(FB.DSC.placenta, selection.method = "vst", nfeatures = 2000)
FB.DSC.placenta.list <- SplitObject(FB.DSC.placenta, split.by = "project")
features <- SelectIntegrationFeatures(object.list = FB.DSC.placenta.list)
FB.DSC.placenta.anchors <- FindIntegrationAnchors(object.list = FB.DSC.placenta.list, anchor.features = features)
FB.DSC.placenta.combined <- IntegrateData(anchorset = FB.DSC.placenta.anchors)
DefaultAssay(FB.DSC.placenta.combined) <- "integrated"
FB.DSC.placenta.combined <- ScaleData(FB.DSC.placenta.combined, features = rownames(FB.DSC.placenta.combined))
FB.DSC.placenta.combined <- RunPCA(FB.DSC.placenta.combined, features = VariableFeatures(FB.DSC.placenta.combined))
ElbowPlot(FB.DSC.placenta.combined, ndims = 50)
FB.DSC.placenta.combined <- RunTSNE(FB.DSC.placenta.combined, reduction = "pca", dims = 1:30)
FB.DSC.placenta.combined <- RunUMAP(FB.DSC.placenta.combined, reduction = "pca", dims = 1:30)
FB.DSC.placenta.combined <- FindNeighbors(FB.DSC.placenta.combined, reduction = "pca", dims = 1:30)
FB.DSC.placenta.combined <- FindClusters(FB.DSC.placenta.combined, resolution = 0.5)
DimPlot(FB.DSC.placenta.combined, label = T)
p1 <- DimPlot(FB.DSC.placenta.combined, group.by = c("project"))
p2 <- DimPlot(FB.DSC.placenta.combined, group.by = c("sample"))
p3 <- DimPlot(FB.DSC.placenta.combined, group.by = c("tech"))
p1 + p2 + p3
FeaturePlot(FB.DSC.placenta.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FB.DSC.placenta.combined <- ScaleData(FB.DSC.placenta.combined, features = rownames(FB.DSC.placenta.combined))
FB.DSC.placenta.combined.markers <- FindAllMarkers(FB.DSC.placenta.combined, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
FB.DSC.placenta.combined.markers.top10 <- FB.DSC.placenta.combined.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(FB.DSC.placenta.combined, features = FB.DSC.placenta.combined.markers.top10$gene)
DefaultAssay(FB.DSC.placenta.combined) <- "RNA"
FeaturePlot(FB.DSC.placenta.combined, features = c("CCL4", "NKG7", "PTPRC"), ncol = 3)
#harmony integration
FB.DSC.placenta <- subset(placenta.combined, subset = seurat_clusters %in% c(3,9,17,5,14))
FB.DSC.placenta <- RenameIdents(FB.DSC.placenta, "3" = "Pla_FB1", "9" = "Pla_FB1", "17" = "Pla_FB1", "5" = "Pla_FB2", "14" = "Pla_FB3")
FB.DSC.placenta$orig.cellname <- FB.DSC.placenta@active.ident
DefaultAssay(FB.DSC.placenta) <- "RNA"
FB.DSC.placenta <- NormalizeData(FB.DSC.placenta)
FB.DSC.placenta <- FindVariableFeatures(FB.DSC.placenta, selection.method = "vst", nfeatures = 2000)
FB.DSC.placenta <- ScaleData(FB.DSC.placenta, features = rownames(FB.DSC.placenta))
FB.DSC.placenta <- RunPCA(FB.DSC.placenta, features = VariableFeatures(FB.DSC.placenta))
FB.DSC.placenta.harmony <- RunHarmony(FB.DSC.placenta, group.by.vars = "project")
ElbowPlot(FB.DSC.placenta.harmony, ndims = 50, reduction = "harmony")
FB.DSC.placenta.harmony <- RunUMAP(FB.DSC.placenta.harmony, reduction = "harmony", dims = 1:30)
FB.DSC.placenta.harmony <- FindNeighbors(FB.DSC.placenta.harmony, reduction = "harmony", dims = 1:30)
FB.DSC.placenta.harmony <- FindClusters(FB.DSC.placenta.harmony, resolution = 0.5)
DimPlot(FB.DSC.placenta.harmony, label = T)
p1 <- DimPlot(FB.DSC.placenta.harmony, group.by = c("project"))
p2 <- DimPlot(FB.DSC.placenta.harmony, group.by = c("sample"))
p3 <- DimPlot(FB.DSC.placenta.harmony, group.by = c("tech"))
p4 <- DimPlot(FB.DSC.placenta.harmony, group.by = "orig.ident")
p1 + p2 + p3 + p4
FeaturePlot(FB.DSC.placenta.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeaturePlot(FB.DSC.placenta.harmony, features = c("CCL4", "NKG7", "PTPRC", "ACTA2"), ncol = 2)
FB.DSC.placenta.harmony.markers <- FindAllMarkers(FB.DSC.placenta.harmony, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
FB.DSC.placenta.harmony.markers.top10 <- FB.DSC.placenta.harmony.markers %>% group_by(cluster) %>% top_n(n = 9, avg_logFC)
DoHeatmap(FB.DSC.placenta.harmony, features = FB.DSC.placenta.harmony.markers.top10$gene)
save(FB.DSC.placenta.harmony, file = "20220622-Comments/02.human-data/science_advance/Placenta-subset-FB-harmony.RData")

#integrate FB in decidua and placenta
decidua.placenta.FB <- merge(FB.DSC.decidua.harmony, FB.DSC.placenta.harmony)
decidua.placenta.FB <- NormalizeData(decidua.placenta.FB, normalization.method = "LogNormalize", scale.factor = 10000)
decidua.placenta.FB <- FindVariableFeatures(decidua.placenta.FB, selection.method = "vst", nfeatures = 2000)
decidua.placenta.FB <- ScaleData(decidua.placenta.FB, features = rownames(decidua.placenta.FB))
decidua.placenta.FB <- RunPCA(decidua.placenta.FB, features = VariableFeatures(decidua.placenta.FB))
ElbowPlot(decidua.placenta.FB, ndims = 50)
decidua.placenta.FB <- FindNeighbors(decidua.placenta.FB, dims = 1:30)
decidua.placenta.FB <- FindClusters(decidua.placenta.FB, resolution = 0.5)
decidua.placenta.FB <- RunUMAP(decidua.placenta.FB, dims = 1:30)
DimPlot(decidua.placenta.FB, label = T)
DimPlot(decidua.placenta.FB, label = T, group.by = "project")
#harmony integration
decidua.placenta.FB <- merge(FB.DSC.decidua.harmony, FB.DSC.placenta.harmony)
DefaultAssay(decidua.placenta.FB) <- "RNA"
decidua.placenta.FB <- NormalizeData(decidua.placenta.FB)
decidua.placenta.FB <- FindVariableFeatures(decidua.placenta.FB, selection.method = "vst", nfeatures = 2000)
decidua.placenta.FB <- ScaleData(decidua.placenta.FB, features = rownames(decidua.placenta.FB))
decidua.placenta.FB <- RunPCA(decidua.placenta.FB, features = VariableFeatures(decidua.placenta.FB))
decidua.placenta.FB.harmony <- RunHarmony(decidua.placenta.FB, group.by.vars = "project")
ElbowPlot(decidua.placenta.FB.harmony, ndims = 50, reduction = "harmony")
decidua.placenta.FB.harmony <- RunUMAP(decidua.placenta.FB.harmony, reduction = "harmony", dims = 1:30)
decidua.placenta.FB.harmony <- FindNeighbors(decidua.placenta.FB.harmony, reduction = "harmony", dims = 1:30)
decidua.placenta.FB.harmony <- FindClusters(decidua.placenta.FB.harmony, resolution = 0.5)
DimPlot(decidua.placenta.FB.harmony, label = T)
decidua.placenta.FB.harmony.markers <- FindAllMarkers(decidua.placenta.FB.harmony, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
decidua.placenta.FB.harmony.markers.top10 <- decidua.placenta.FB.harmony.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(decidua.placenta.FB.harmony, features = decidua.placenta.FB.harmony.markers.top10$gene)
#CCA integration
decidua.placenta.FB <- merge(FB.DSC.decidua.harmony, FB.DSC.placenta.harmony)
decidua.placenta.FB <- NormalizeData(decidua.placenta.FB, normalization.method = "LogNormalize", scale.factor = 10000)
decidua.placenta.FB <- FindVariableFeatures(decidua.placenta.FB, selection.method = "vst", nfeatures = 2000)
decidua.placenta.FB.list <- SplitObject(decidua.placenta.FB, split.by = "project")
features <- SelectIntegrationFeatures(object.list = decidua.placenta.FB.list)
decidua.placenta.FB.anchors <- FindIntegrationAnchors(object.list = decidua.placenta.FB.list, anchor.features = features)
decidua.placenta.FB.combined <- IntegrateData(anchorset = decidua.placenta.FB.anchors)
DefaultAssay(decidua.placenta.FB.combined) <- "integrated"
decidua.placenta.FB.combined <- ScaleData(decidua.placenta.FB.combined, features = rownames(decidua.placenta.FB.combined))
decidua.placenta.FB.combined <- RunPCA(decidua.placenta.FB.combined, features = VariableFeatures(decidua.placenta.FB.combined))
ElbowPlot(decidua.placenta.FB.combined, ndims = 50)
decidua.placenta.FB.combined <- RunTSNE(decidua.placenta.FB.combined, reduction = "pca", dims = 1:30)
decidua.placenta.FB.combined <- RunUMAP(decidua.placenta.FB.combined, reduction = "pca", dims = 1:30)
decidua.placenta.FB.combined <- FindNeighbors(decidua.placenta.FB.combined, dims = 1:30)
decidua.placenta.FB.combined <- FindClusters(decidua.placenta.FB.combined, resolution = 0.5)
DimPlot(decidua.placenta.FB.combined, group.by = "orig.cellname")
DimPlot(decidua.placenta.FB.combined, label = T)
decidua.placenta.FB.combined.markers <- FindAllMarkers(decidua.placenta.FB.combined, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
decidua.placenta.FB.combined.markers.top10 <- decidua.placenta.FB.combined.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(decidua.placenta.FB.combined, features = decidua.placenta.FB.combined.markers.top10$gene)












