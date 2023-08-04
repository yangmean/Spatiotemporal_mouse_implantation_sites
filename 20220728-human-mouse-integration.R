#convert human gene names to mouse gene names
load("FB-immune-DSC-harmony-20211214.RData")
load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-harmony.RData")#Science Advance data
load("../../20220508-download-data/Nature-2018/Human-nature-2018-DSC-harmony-20220510.RData")#Nature data
DSC.human.gene <- union(rownames(FB.DSC.decidua.harmony), (rownames(human.S.harmony)))
library("biomaRt")
#human.gene <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#mouse.gene <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human.gene <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse.gene <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                  values = DSC.human.gene , mart = human.gene, 
                  attributesL = c("mgi_symbol"), martL = mouse.gene, uniqueRows=T)
save(genesV2, file = "20220622-Comments/02.human-data/DSC-human-mouse-gene-transfer.RData")
colnames(genesV2) <- c("human.gene", "mouse.gene")
genesV2 <- genesV2[!duplicated(genesV2$human.gene), ]

#read mouse data
load("FB-immune-DSC-harmony-20211214.RData")
unique(UE.FB@active.ident)
UE.FB$cellname <- UE.FB@active.ident
UE.FB$literature <- "My.mouse"

#extract human DSC data
load("../../20220508-download-data/Nature-2018/Human-nature-2018-DSC-harmony-20220510.RData")#Nature data
human.S.harmony$cellname <- human.S.harmony$annotation
human.S.harmony$time <- human.S.harmony$Fetus
human.S.harmony$literature <- "Nature.human"
human.S.harmony.meta <- human.S.harmony@meta.data
human.S.harmony.data <- data.frame(as.matrix(GetAssayData(human.S.harmony, slot = "counts")), check.names = F)
human.S.harmony.data$human.gene <- rownames(human.S.harmony.data)
human.S.harmony.data <- merge(human.S.harmony.data, genesV2, by = "human.gene")
human.S.harmony.data <- human.S.harmony.data[!duplicated(human.S.harmony.data$mouse.gene), ]
rownames(human.S.harmony.data) <- human.S.harmony.data$mouse.gene
human.S.harmony.data$human.gene <- NULL
human.S.harmony.data$mouse.gene <- NULL
DSC.Nature <- CreateSeuratObject(human.S.harmony.data, meta.data = human.S.harmony.meta, project = "human.nature")

load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-harmony.RData")#Science Advance data
DimPlot(FB.DSC.decidua.harmony)
FB.DSC.decidua.harmony$literature <- "SD.human"
FB.DSC.decidua.harmony$cellname <- FB.DSC.decidua.harmony$orig.cellname
FB.DSC.decidua.harmony$orig.cellname <- NULL
FB.DSC.decidua.harmony$time <- FB.DSC.decidua.harmony$orig.ident
FB.DSC.decidua.harmony.meta <- FB.DSC.decidua.harmony@meta.data
FB.DSC.decidua.harmony.data <- data.frame(as.matrix(GetAssayData(FB.DSC.decidua.harmony, slot = "counts")), check.names = F)
FB.DSC.decidua.harmony.data$human.gene <- rownames(FB.DSC.decidua.harmony.data)
FB.DSC.decidua.harmony.data <- merge(FB.DSC.decidua.harmony.data, genesV2, by = "human.gene")
FB.DSC.decidua.harmony.data <- FB.DSC.decidua.harmony.data[!duplicated(FB.DSC.decidua.harmony.data$mouse.gene), ]
rownames(FB.DSC.decidua.harmony.data) <- FB.DSC.decidua.harmony.data$mouse.gene
FB.DSC.decidua.harmony.data$human.gene <- NULL
FB.DSC.decidua.harmony.data$mouse.gene <- NULL
DSC.SD <- CreateSeuratObject(FB.DSC.decidua.harmony.data, meta.data = FB.DSC.decidua.harmony.meta, project = "human.SD")

#harmony integration
DSC.human.mouse <- merge(UE.FB.sub, y = list(DSC.Nature, DSC.SD))
DSC.human.mouse <- NormalizeData(DSC.human.mouse)
DSC.human.mouse <- FindVariableFeatures(DSC.human.mouse, selection.method = "vst", nfeatures = 2000)
DSC.human.mouse <- ScaleData(DSC.human.mouse, features = rownames(DSC.human.mouse))
DSC.human.mouse <- RunPCA(DSC.human.mouse, features = VariableFeatures(DSC.human.mouse))
DSC.human.mouse.harmony <- RunHarmony(DSC.human.mouse, group.by.vars = "literature")
ElbowPlot(DSC.human.mouse.harmony, ndims = 50, reduction = "harmony")
DSC.human.mouse.harmony <- RunUMAP(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindNeighbors(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindClusters(DSC.human.mouse.harmony, resolution = 0.5)
DimPlot(DSC.human.mouse.harmony, label = T)
DimPlot(DSC.human.mouse.harmony, label = T, group.by = "literature")
DimPlot(DSC.human.mouse.harmony, label = T, group.by = "cellname")
save(DSC.human.mouse.harmony, file = "20220622-Comments/DSC-human-mouse-harmony.RData")

#Seurat3 integration
DSC.human.mouse.list <- SplitObject(DSC.human.mouse, split.by = "literature")
DSC.human.mouse.list <- lapply(X = DSC.human.mouse.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = DSC.human.mouse.list)
DSC.human.mouse.anchors <- FindIntegrationAnchors(object.list = DSC.human.mouse.list, anchor.features = features)
DSC.human.mouse.combined <- IntegrateData(anchorset = DSC.human.mouse.anchors)
DefaultAssay(DSC.human.mouse.combined) <- "integrated"
DSC.human.mouse.combined <- ScaleData(DSC.human.mouse.combined, features = rownames(DSC.human.mouse.combined))
DSC.human.mouse.combined <- RunPCA(DSC.human.mouse.combined, features = VariableFeatures(DSC.human.mouse.combined))
ElbowPlot(DSC.human.mouse.combined, ndims = 50)
DSC.human.mouse.combined <- RunUMAP(DSC.human.mouse.combined, reduction = "pca", dims = 1:30)
DSC.human.mouse.combined <- FindNeighbors(DSC.human.mouse.combined, reduction = "pca", dims = 1:30)
DSC.human.mouse.combined <- FindClusters(DSC.human.mouse.combined, resolution = 0.5)
DimPlot(DSC.human.mouse.combined, label = T)
DimPlot(DSC.human.mouse.combined, label = T, group.by = "literature")
DimPlot(DSC.human.mouse.combined, label = T, group.by = "cellname")
save(DSC.human.mouse.combined, file = "20220622-Comments/DSC-human-mouse-SeuratV3.RData")

#SCTransform
DSC.human.mouse.list <- lapply(X = DSC.human.mouse.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = DSC.human.mouse.list, nfeatures = 3000)
DSC.human.mouse.list <- PrepSCTIntegration(object.list = DSC.human.mouse.list, anchor.features = features)
DSC.human.mouse.anchors <- FindIntegrationAnchors(object.list = DSC.human.mouse.list, normalization.method = "SCT",
                                                  anchor.features = features)
DSC.human.mouse.sct <- IntegrateData(anchorset = DSC.human.mouse.anchors, normalization.method = "SCT")
DSC.human.mouse.sct <- RunPCA(DSC.human.mouse.sct, verbose = FALSE)
DSC.human.mouse.sct <- RunUMAP(DSC.human.mouse.sct, reduction = "pca", dims = 1:30)
DSC.human.mouse.sct <- FindNeighbors(DSC.human.mouse.sct, reduction = "pca", dims = 1:30)
DSC.human.mouse.sct <- FindClusters(DSC.human.mouse.sct, resolution = 0.5)
save(DSC.human.mouse.sct, file = "20220622-Comments/DSC-human-mouse-SCTransform.RData")
DimPlot(DSC.human.mouse.sct, label = T)
DimPlot(DSC.human.mouse.sct, label = T, group.by = "literature")
DimPlot(DSC.human.mouse.sct, label = T, group.by = "cellname")

#subset 1000 cells per cluster of mouse DSCs
UE.FB.sub <- subset(UE.FB, downsample = 2000)
DSC.human.mouse <- merge(UE.FB.sub, y = list(DSC.Nature, DSC.SD))
variable.genes <- intersect(rownames(DSC.human.mouse), VariableFeatures(UE.FB))
DSC.human.mouse.list <- lapply(X = list(UE.FB, DSC.Nature, DSC.SD), FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = DSC.human.mouse.list, nfeatures = 2000)
DSC.human.mouse.list <- PrepSCTIntegration(object.list = DSC.human.mouse.list, anchor.features = features)
DSC.human.mouse.anchors <- FindIntegrationAnchors(object.list = DSC.human.mouse.list, normalization.method = "SCT",
                                                  anchor.features = features)
DSC.human.mouse.sct <- IntegrateData(anchorset = DSC.human.mouse.anchors, normalization.method = "SCT")
DSC.human.mouse.sct <- RunPCA(DSC.human.mouse.sct, verbose = FALSE, features = variable.genes)
ElbowPlot(DSC.human.mouse.sct, ndims = 50)
DSC.human.mouse.sct <- RunUMAP(DSC.human.mouse.sct, reduction = "pca", dims = 1:30)
DSC.human.mouse.sct <- FindNeighbors(DSC.human.mouse.sct, reduction = "pca", dims = 1:30)
DSC.human.mouse.sct <- FindClusters(DSC.human.mouse.sct, resolution = 0.5)
save(DSC.human.mouse.sct, file = "20220622-Comments/DSC-human-mouse-SCTransform-subset.RData")
DimPlot(DSC.human.mouse.sct, label = T)
DimPlot(DSC.human.mouse.sct, label = T, group.by = "literature")
DimPlot(DSC.human.mouse.sct, label = T, group.by = "cellname", split.by = "literature")
DSC.human.mouse.sct <- ScaleData(DSC.human.mouse.sct, features = rownames(DSC.human.mouse.sct))
DSC.human.mouse.sct.markers <- FindAllMarkers(DSC.human.mouse.sct, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T, assay = "SCT")
DSC.human.mouse.sct.markers.top8 <- DSC.human.mouse.sct.markers %>% group_by(cluster) %>% top_n(n = 8, avg_logFC)
DoHeatmap(DSC.human.mouse.sct, features = DSC.human.mouse.sct.markers.top8$gene)

DSC.human.mouse.list <- lapply(X = list(UE.FB.sub, DSC.Nature, DSC.SD), FUN = SCTransform)
DSC.human.mouse <- merge(DSC.human.mouse.list[[1]], y = DSC.human.mouse.list[2:3])
#DSC.human.mouse <- NormalizeData(DSC.human.mouse)
#DSC.human.mouse <- FindVariableFeatures(DSC.human.mouse, selection.method = "vst", nfeatures = 2000)
#DSC.human.mouse <- ScaleData(DSC.human.mouse, features = rownames(DSC.human.mouse))
DSC.human.mouse <- RunPCA(DSC.human.mouse, features = VariableFeatures(UE.FB), assay = "SCT")
DSC.human.mouse.harmony <- RunHarmony(DSC.human.mouse, group.by.vars = "literature", assay.use = "SCT")
ElbowPlot(DSC.human.mouse.harmony, ndims = 50, reduction = "harmony")
DSC.human.mouse.harmony <- RunUMAP(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindNeighbors(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindClusters(DSC.human.mouse.harmony, resolution = 0.5)
DimPlot(DSC.human.mouse.harmony, label = T)
DimPlot(DSC.human.mouse.harmony, label = T, group.by = "literature")
DimPlot(DSC.human.mouse.harmony, label = T, group.by = "cellname")

#20220727
#按照time进行harmony 整合
DSC.human.mouse <- merge(UE.FB, y = list(DSC.Nature, DSC.SD))
DSC.human.mouse <- NormalizeData(DSC.human.mouse)
DSC.human.mouse <- FindVariableFeatures(DSC.human.mouse, selection.method = "vst", nfeatures = 2000)
DSC.human.mouse <- ScaleData(DSC.human.mouse, features = rownames(DSC.human.mouse))
DSC.human.mouse <- RunPCA(DSC.human.mouse, features = VariableFeatures(DSC.human.mouse))
DSC.human.mouse.harmony <- RunHarmony(DSC.human.mouse, group.by.vars = "time")
ElbowPlot(DSC.human.mouse.harmony, ndims = 50, reduction = "harmony")
DSC.human.mouse.harmony <- RunUMAP(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindNeighbors(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindClusters(DSC.human.mouse.harmony, resolution = 0.5)
DimPlot(DSC.human.mouse.harmony, label = T)
DimPlot(DSC.human.mouse.harmony, label = T, group.by = "literature")
DimPlot(DSC.human.mouse.harmony, label = T, group.by = "cellname", split.by = "literature")

#SCTransform + harmony integration by time
DSC.human.mouse.list <- lapply(X = list(UE.FB, DSC.Nature, DSC.SD), FUN = SCTransform)
DSC.human.mouse <- merge(DSC.human.mouse.list[[1]], y = DSC.human.mouse.list[2:3])
DSC.human.mouse <- FindVariableFeatures(DSC.human.mouse, selection.method = "vst", nfeatures = 2000)
DSC.human.mouse <- RunPCA(DSC.human.mouse, features = VariableFeatures(DSC.human.mouse), assay = "SCT")
DSC.human.mouse.harmony <- RunHarmony(DSC.human.mouse, group.by.vars = "time", assay.use = "SCT")
ElbowPlot(DSC.human.mouse.harmony, ndims = 50, reduction = "harmony")
DSC.human.mouse.harmony <- RunUMAP(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindNeighbors(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindClusters(DSC.human.mouse.harmony, resolution = 0.5)
DimPlot(DSC.human.mouse.harmony, label = T)
DimPlot(DSC.human.mouse.harmony, label = T, group.by = "literature")
DimPlot(DSC.human.mouse.harmony, label = T, group.by = "cellname", split.by = "literature")

#RPCA
DSC.human.mouse.list <- lapply(X = list(UE.FB, DSC.Nature, DSC.SD), FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(DSC.human.mouse.list)
DSC.human.mouse.list <- lapply(X = list(UE.FB, DSC.Nature, DSC.SD), FUN = function(x){
  x <- ScaleData(x, features = features, verbose = F)
  x <- RunPCA(x, features = features)
})
DSC.human.mouse.anchors <- FindIntegrationAnchors(DSC.human.mouse.list, anchor.features = features, reduction = "rpca")
DSC.human.mouse.combined <- IntegrateData(anchorset = DSC.human.mouse.anchors)
DefaultAssay(DSC.human.mouse.combined) <- "integrated"

DSC.human.mouse.combined <- ScaleData(DSC.human.mouse.combined, verbose = FALSE)
DSC.human.mouse.combined <- RunPCA(DSC.human.mouse.combined, npcs = 50, verbose = FALSE)
ElbowPlot(DSC.human.mouse.combined, ndims = 50)
DSC.human.mouse.combined <- RunUMAP(DSC.human.mouse.combined, reduction = "pca", dims = 1:30)
DSC.human.mouse.combined <- FindNeighbors(DSC.human.mouse.combined, reduction = "pca", dims = 1:30)
DSC.human.mouse.combined <- FindClusters(DSC.human.mouse.combined, resolution = 0.5)
DimPlot(DSC.human.mouse.combined, label = T)
DimPlot(DSC.human.mouse.combined, label = T, group.by = "literature")
DimPlot(DSC.human.mouse.combined, label = T, group.by = "cellname", split.by = "literature")

#integrate 2 datasets of human - CCA
load("20220622-Comments/02.human-data/nature/Nature-human-DSC-recluster-CCA.RData")
load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-CCA.RData")
FB.DSC.decidua.combined$literature <- "SD"
FB.DSC.decidua.combined.sub <- subset(FB.DSC.decidua.combined, downsample = 500)
human.S.combined$literature <- "Nature"

features <- SelectIntegrationFeatures(object.list = c(FB.DSC.decidua.combined, human.S.combined))
human.decidua.anchors <- FindIntegrationAnchors(object.list = c(FB.DSC.decidua.combined, human.S.combined), anchor.features = features)
human.decidua.combined <- IntegrateData(anchorset = human.decidua.anchors)
DefaultAssay(human.decidua.combined) <- "integrated"
human.decidua.combined <- ScaleData(human.decidua.combined, features = rownames(human.decidua.combined))
human.decidua.combined <- RunPCA(human.decidua.combined, features = VariableFeatures(human.decidua.combined))
ElbowPlot(human.decidua.combined, ndims = 50)
human.decidua.combined <- RunUMAP(human.decidua.combined, reduction = "pca", dims = 1:30)
human.decidua.combined <- FindNeighbors(human.decidua.combined, reduction = "pca", dims = 1:30)
human.decidua.combined <- FindClusters(human.decidua.combined, resolution = 0.3)
DimPlot(human.decidua.combined, group.by = "literature")
DimPlot(human.decidua.combined, label = T)
DimPlot(human.decidua.combined, group.by = "cluster1", cols = cluster.colors[c(2:6,9:17)], split.by = "literature")
DimPlot(human.decidua.combined, group.by = "orig.cellname", cols = cluster.colors[c(2:6,9:17)], split.by = "literature")
human.decidua.combined.markers <- FindAllMarkers(human.decidua.combined, min.pct = 0.5, logfc.threshold = 0.25, only.pos = T)
human.decidua.combined.markers.top10 <- human.decidua.combined.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(human.decidua.combined, features = human.decidua.combined.markers.top10$gene)

#subset 500 cells/cluster
set.seed(0)
human.S.combined.sub <- subset(human.S.combined, downsample = 500)
features <- SelectIntegrationFeatures(object.list = c(FB.DSC.decidua.combined.sub, human.S.combined.sub))
human.decidua.anchors <- FindIntegrationAnchors(object.list = c(FB.DSC.decidua.combined.sub, human.S.combined.sub), anchor.features = features)
human.decidua.combined.sub <- IntegrateData(anchorset = human.decidua.anchors)
DefaultAssay(human.decidua.combined.sub) <- "integrated"
human.decidua.combined.sub <- ScaleData(human.decidua.combined.sub, features = rownames(human.decidua.combined.sub))
human.decidua.combined.sub <- RunPCA(human.decidua.combined.sub, features = VariableFeatures(human.decidua.combined.sub))
ElbowPlot(human.decidua.combined.sub, ndims = 50)
human.decidua.combined.sub <- RunUMAP(human.decidua.combined.sub, reduction = "pca", dims = 1:30)
human.decidua.combined.sub <- FindNeighbors(human.decidua.combined.sub, reduction = "pca", dims = 1:30)
human.decidua.combined.sub <- FindClusters(human.decidua.combined.sub, resolution = 0.3)
DimPlot(human.decidua.combined.sub, group.by = "literature")
DimPlot(human.decidua.combined.sub, label = T)
DimPlot(human.decidua.combined.sub, group.by = "cluster1", cols = cluster.colors[c(2:6,9:17)], split.by = "literature")
DimPlot(human.decidua.combined.sub, group.by = "orig.cellname", cols = cluster.colors[c(2:6,9:17)], split.by = "literature")
human.decidua.combined.sub.markers <- FindAllMarkers(human.decidua.combined.sub, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
human.decidua.combined.sub.markers.top10 <- human.decidua.combined.sub.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(human.decidua.combined.sub, features = human.decidua.combined.sub.markers.top10$gene)

#integrate 2 datasets of human - harmony
load("20220622-Comments/02.human-data/nature/Nature-human-DSC-recluster-CCA.RData")
load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-CCA.RData")
human.decidua <- merge(FB.DSC.decidua.combined, human.S.combined)
DefaultAssay(human.decidua) <- "RNA"
human.decidua <- NormalizeData(human.decidua)
human.decidua <- FindVariableFeatures(human.decidua, selection.method = "vst", nfeatures = 2000)
human.decidua <- ScaleData(human.decidua, features = rownames(human.decidua))
human.decidua <- RunPCA(human.decidua, features = VariableFeatures(human.decidua))
human.decidua.harmony <- RunHarmony(human.decidua, group.by.vars = "literature")
ElbowPlot(human.decidua.harmony, ndims = 50, reduction = "harmony")
human.decidua.harmony <- RunUMAP(human.decidua.harmony, reduction = "harmony", dims = 1:30)
human.decidua.harmony <- FindNeighbors(human.decidua.harmony, reduction = "harmony", dims = 1:30)
human.decidua.harmony <- FindClusters(human.decidua.harmony, resolution = 0.5)
DimPlot(human.decidua.harmony, label = T)
DimPlot(human.decidua.harmony, group.by = "literature")
DimPlot(human.decidua.harmony, group.by = "cluster1", split.by = "literature")
DimPlot(human.decidua.harmony, group.by = "orig.cellname", split.by = "literature")

####################################20221009##############################################
#Nature
load("20220622-Comments/02.human-data/nature/Nature-human-DSC-recluster-CCA.RData")
pdf("20220622-Comments/02.human-data/00.figures/01.Human-Nature-DSC-subclusters-UMAP.pdf", width = 10, height = 10)
DimPlot(human.S.combined, label = F, pt.size = 1)
dev.off()
human.S.combined.markers <- FindAllMarkers(human.S.combined, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
human.S.combined.markers.top20 <- human.S.combined.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
DefaultAssay(human.S.combined) <- "RNA"
human.S.combined <- ScaleData(human.S.combined, features = rownames(human.S.combined))
pdf("20220622-Comments/02.human-data/00.figures/02.Human-Nature-DSC-subclusters-DEGs-top20-heatmap.pdf", width = 10, height = 10)
DoHeatmap(human.S.combined, features = human.S.combined.markers.top20$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()

human.nature.cluster5.markers <- human.S.combined.markers[human.S.combined.markers$cluster == 5, ]
human.S.combined <- AddModuleScore(human.S.combined, features = list(human.nature.cluster5.markers$gene), name = "Immune_features")
pdf("20220622-Comments/02.human-data/00.figures/03.Human-Nature-iDSC-immune-features-avg-featureplot.pdf", width = 10, height = 10)
FeaturePlot(human.S.combined, features = "Immune_features1", min.cutoff = 0.1, max.cutoff = 0.6, pt.size = 1) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()
#examples
human.features <- c("ACTA2", "IGFBP1", "PDGFRA", "CCL4", "CD74", "NKG7", "GZMA", "GZMB", "CCL3", 
                    "HLA-DRA", "XCL2", "C1QA", "C1QB", "SPP1", "FCER1G")
pdf("20220622-Comments/02.human-data/00.figures/04.Human-Nature-iDSC-example-genes-featureplot.pdf", width = 20, height = 20)
FeaturePlot(human.S.combined, features = human.features, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()
#Science Advance
load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-CCA.RData")
pdf("20220622-Comments/02.human-data/00.figures/01.Human-SA-DSC-subclusters-UMAP.pdf", width = 10, height = 10)
DimPlot(FB.DSC.decidua.combined, label = F, pt.size = 2)
dev.off()
FB.DSC.decidua.combined.markers <- FindAllMarkers(FB.DSC.decidua.combined, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
FB.DSC.decidua.combined.markers.top20 <- FB.DSC.decidua.combined.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
DefaultAssay(FB.DSC.decidua.combined) <- "RNA"
FB.DSC.decidua.combined <- ScaleData(FB.DSC.decidua.combined, features = rownames(FB.DSC.decidua.combined))
FB.DSC.decidua.combined.cluster5.markers <- FB.DSC.decidua.combined.markers[FB.DSC.decidua.combined.markers$cluster == 5, ]
FB.DSC.decidua.combined <- AddModuleScore(FB.DSC.decidua.combined, features = list(FB.DSC.decidua.combined.cluster5.markers$gene), name = "Immune_features")
pdf("20220622-Comments/02.human-data/00.figures/03.Human-SA-iDSC-immune-features-avg-featureplot.pdf", width = 10, height = 10)
FeaturePlot(d, features = "Immune_features1", pt.size = 2, min.cutoff = 0.05, max.cutoff = 0.6) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

pdf("20220622-Comments/02.human-data/00.figures/02.Human-SA-DSC-subclusters-DEGs-top20-heatmap.pdf", width = 10, height = 10)
DoHeatmap(FB.DSC.decidua.combined, features = FB.DSC.decidua.combined.markers.top20$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()
pdf("20220622-Comments/02.human-data/00.figures/04.Human-SA-iDSC-example-genes-featureplot.pdf", width = 20, height = 20)
FeaturePlot(FB.DSC.decidua.combined, features = human.features, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

#integrate 2 human datasets, downsample 500 cells/cluster
load("20220622-Comments/02.human-data/nature/Nature-human-DSC-recluster-CCA.RData")
load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-CCA.RData")
FB.DSC.decidua.combined$literature <- "SD"
FB.DSC.decidua.combined$primary_cluster <- paste("SD_", FB.DSC.decidua.combined$seurat_clusters)
set.seed(0)
FB.DSC.decidua.combined.sub <- subset(FB.DSC.decidua.combined, downsample = 500)
human.S.combined$literature <- "Nature"
human.S.combined$primary_cluster <- paste("Nature_", human.S.combined$seurat_clusters)
set.seed(0)
human.S.combined.sub <- subset(human.S.combined, downsample = 500)

features <- SelectIntegrationFeatures(object.list = c(FB.DSC.decidua.combined.sub, human.S.combined.sub))
human.decidua.anchors <- FindIntegrationAnchors(object.list = c(FB.DSC.decidua.combined.sub, human.S.combined.sub), anchor.features = features)
human.decidua.combined.sub <- IntegrateData(anchorset = human.decidua.anchors, features.to.integrate = intersect(rownames(FB.DSC.decidua.combined.sub), rownames(human.S.combined.sub)))
DefaultAssay(human.decidua.combined.sub) <- "integrated"
human.decidua.combined.sub <- ScaleData(human.decidua.combined.sub, features = rownames(human.decidua.combined.sub))
human.decidua.combined.sub <- RunPCA(human.decidua.combined.sub, features = VariableFeatures(human.decidua.combined.sub))
ElbowPlot(human.decidua.combined.sub, ndims = 50)
human.decidua.combined.sub <- RunUMAP(human.decidua.combined.sub, reduction = "pca", dims = 1:50)
human.decidua.combined.sub <- FindNeighbors(human.decidua.combined.sub, reduction = "pca", dims = 1:50)
human.decidua.combined.sub <- FindClusters(human.decidua.combined.sub, resolution = 0.3)
pdf("20220622-Comments/02.human-data/00.figures/05.Human-Nature-SA-DSC-subclusters-integrate-UMAP-sub5000.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub, group.by = "literature", pt.size = 1, cols = c("#9370DB", "#3CB371"))
dev.off()
DimPlot(human.decidua.combined.sub, label = T)
DimPlot(human.decidua.combined.sub, group.by = "literature")
DefaultAssay(human.decidua.combined.sub) <- "RNA"
human.decidua.combined.sub <- ScaleData(human.decidua.combined.sub, features = rownames(human.decidua.combined.sub))
pdf("20220622-Comments/02.human-data/00.figures/06.Human-Nature-SA-DSC-subclusters-integrate-iDSC-example-genes-featureplot.pdf", width = 20, height = 20)
FeaturePlot(human.decidua.combined.sub, features = human.features, pt.size = 0.5, ncol = 4) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

DefaultAssay(human.decidua.combined.sub) <- "integrated"
human.decidua.combined.sub.markers <- FindAllMarkers(human.decidua.combined.sub, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
human.decidua.combined.sub.markers.top <- human.decidua.combined.sub.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
DoHeatmap(human.decidua.combined.sub, features = human.decidua.combined.sub.markers.top$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
DefaultAssay(human.decidua.combined.sub) <- "RNA"
human.decidua.combined.sub <- ScaleData(human.decidua.combined.sub, features = rownames(human.decidua.combined.sub))
human.decidua.combined.sub.cluster4.markers <- human.decidua.combined.sub.markers[human.decidua.combined.sub.markers$cluster == 4, ]
human.decidua.combined.sub <- AddModuleScore(human.decidua.combined.sub, features = list(human.decidua.combined.sub.cluster4.markers$gene), name = "Immune_features")
pdf("20220622-Comments/02.human-data/00.figures/03.Human-Nature-SD-iDSC-immune-features-avg-featureplot.pdf", width = 10, height = 10)
FeaturePlot(human.decidua.combined.sub, features = "Immune_features1", pt.size = 1, min.cutoff = 0.1, max.cutoff = 0.5) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()
save(human.decidua.combined.sub, file = "20220622-Comments/02.human-data/Human-Nature-SD-DSCs-downsample-500-integrate.RData")

############20221101 human and mouse integration DSCs and iDSCs###########################
#extract human DSC data
load("../../20220508-download-data/Nature-2018/Human-nature-2018-DSC-harmony-20220510.RData")#Nature data
human.S.harmony$cellname <- human.S.harmony$annotation
human.S.harmony$time <- human.S.harmony$Fetus
human.S.harmony$literature <- "Nature.human"
human.S.harmony.meta <- human.S.harmony@meta.data
human.S.harmony.data <- data.frame(as.matrix(GetAssayData(human.S.harmony, slot = "counts")), check.names = F)
human.S.harmony.data$human.gene <- rownames(human.S.harmony.data)
human.S.harmony.data <- merge(human.S.harmony.data, genesV2, by = "human.gene")
human.S.harmony.data <- human.S.harmony.data[!duplicated(human.S.harmony.data$mouse.gene), ]
rownames(human.S.harmony.data) <- human.S.harmony.data$mouse.gene
human.S.harmony.data$human.gene <- NULL
human.S.harmony.data$mouse.gene <- NULL
DSC.Nature <- CreateSeuratObject(human.S.harmony.data, meta.data = human.S.harmony.meta, project = "human.nature")

load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-harmony.RData")#Science Advance data
DimPlot(FB.DSC.decidua.harmony)
FB.DSC.decidua.harmony$literature <- "SD.human"
FB.DSC.decidua.harmony$cellname <- FB.DSC.decidua.harmony$orig.cellname
FB.DSC.decidua.harmony$orig.cellname <- NULL
FB.DSC.decidua.harmony$time <- FB.DSC.decidua.harmony$orig.ident
FB.DSC.decidua.harmony.meta <- FB.DSC.decidua.harmony@meta.data
FB.DSC.decidua.harmony.data <- data.frame(as.matrix(GetAssayData(FB.DSC.decidua.harmony, slot = "counts")), check.names = F)
FB.DSC.decidua.harmony.data$human.gene <- rownames(FB.DSC.decidua.harmony.data)
FB.DSC.decidua.harmony.data <- merge(FB.DSC.decidua.harmony.data, genesV2, by = "human.gene")
FB.DSC.decidua.harmony.data <- FB.DSC.decidua.harmony.data[!duplicated(FB.DSC.decidua.harmony.data$mouse.gene), ]
rownames(FB.DSC.decidua.harmony.data) <- FB.DSC.decidua.harmony.data$mouse.gene
FB.DSC.decidua.harmony.data$human.gene <- NULL
FB.DSC.decidua.harmony.data$mouse.gene <- NULL
DSC.SD <- CreateSeuratObject(FB.DSC.decidua.harmony.data, meta.data = FB.DSC.decidua.harmony.meta, project = "human.SD")

#harmony integration nature/SD/mouse normal samples
DSC.human.mouse <- merge(UE.FB, y = list(DSC.Nature, DSC.SD))
DSC.human.mouse <- NormalizeData(DSC.human.mouse)
DSC.human.mouse <- FindVariableFeatures(DSC.human.mouse, selection.method = "vst", nfeatures = 2000)
DSC.human.mouse <- ScaleData(DSC.human.mouse, features = rownames(DSC.human.mouse))
DSC.human.mouse <- RunPCA(DSC.human.mouse, features = VariableFeatures(DSC.human.mouse))
DSC.human.mouse.harmony <- RunHarmony(DSC.human.mouse, group.by.vars = "time")
ElbowPlot(DSC.human.mouse.harmony, ndims = 50, reduction = "harmony")
DSC.human.mouse.harmony <- RunUMAP(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindNeighbors(DSC.human.mouse.harmony, reduction = "harmony", dims = 1:30)
DSC.human.mouse.harmony <- FindClusters(DSC.human.mouse.harmony, resolution = 0.5)
literature.color <- c("#FF8C00", "#9370DB", "#3CB371")
pdf("20220622-Comments/02.human-data/00.figures/07.Human-mouse-harmony-integrate-UMAP-color-project.pdf", width = 11, height = 10)
DimPlot(DSC.human.mouse.harmony, label = F, group.by = "literature", pt.size = 0.4, cols = literature.color)
dev.off()
DSC.human.mouse.harmony$cellname <- factor(DSC.human.mouse.harmony$cellname, 
                                           levels = c(decidual.levels[1:7], "FB-immune", "dS1", "dS2", "dS3", "DSC", "FB1", 'FB2'))
pdf("20220622-Comments/02.human-data/00.figures/08.Human-mouse-harmony-integrate-UMAP-color-mouse-clusters.pdf", width = 11, height = 10)
DimPlot(DSC.human.mouse.harmony, label = F, pt.size = 0.4, group.by = "cellname", 
        cols = c(decidual.color[1:7], immune.color[10], rep("#9370DB", 3), rep("#3CB371", 3)))
dev.off()
save(DSC.human.mouse.harmony, file = "20220622-Comments/02.human-data/Human-mouse-DSC-integrate.RData")
#subset FB-immune to reclustering
load("20220622-Comments/02.human-data/nature/Nature-human-DSC-recluster-CCA.RData")
load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-CCA.RData")
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
DimPlot(FB.immune.SCT)
DimPlot(human.S.combined, label = F, pt.size = 1)
DimPlot(FB.DSC.decidua.combined, label = T, pt.size = 1)

human.S.combined$literature <- "Nature"
Nature.iDSC <- subset(human.S.combined, seurat_clusters %in% c(5))
Nature.iDSC.meta <- Nature.iDSC@meta.data
Nature.iDSC.data <- data.frame(as.matrix(Nature.iDSC@assays$RNA@counts), check.names = F)
Nature.iDSC.data$human.gene <- rownames(Nature.iDSC.data)
Nature.iDSC.data <- merge(Nature.iDSC.data, genesV2, by = "human.gene")
Nature.iDSC.data <- Nature.iDSC.data[!duplicated(Nature.iDSC.data$mouse.gene), ]
rownames(Nature.iDSC.data) <- Nature.iDSC.data$mouse.gene
Nature.iDSC.data$human.gene <- NULL
Nature.iDSC.data$mouse.gene <- NULL
iDSC.Nature <- CreateSeuratObject(Nature.iDSC.data, meta.data = Nature.iDSC.meta, project = "human.nature")

FB.DSC.decidua.combined$literature <- "SD"
SD.iDSC <- subset(FB.DSC.decidua.combined, seurat_clusters %in% c(5))
SD.iDSC.meta <- SD.iDSC@meta.data
SD.iDSC.data <- data.frame(as.matrix(SD.iDSC@assays$RNA@counts), check.names = F)
SD.iDSC.data$human.gene <- rownames(SD.iDSC.data)
SD.iDSC.data <- merge(SD.iDSC.data, genesV2, by = "human.gene")
SD.iDSC.data <- SD.iDSC.data[!duplicated(SD.iDSC.data$mouse.gene), ]
rownames(SD.iDSC.data) <- SD.iDSC.data$mouse.gene
SD.iDSC.data$human.gene <- NULL
SD.iDSC.data$mouse.gene <- NULL
iDSC.SD <- CreateSeuratObject(SD.iDSC.data, meta.data = SD.iDSC.meta, project = "human.SD")

#harmony integrate
FB.immune.SCT$literature <- "my.mouse"
FB.immune.SCT$cellname <- paste("iDSC_", FB.immune.SCT$seurat_clusters, sep = "")
iDSC.Nature$cellname <- "Nature"
iDSC.Nature$time <- "Nature"
iDSC.SD$cellname <- "SD"
iDSC.SD$time <- "SD"
Human.mouse.iDSC.list <- lapply(X = list(FB.immune.SCT, iDSC.Nature, iDSC.SD), FUN = SCTransform)
Human.mouse.iDSC <- merge(Human.mouse.iDSC.list[[1]], y = Human.mouse.iDSC.list[2:3])
Human.mouse.iDSC <- FindVariableFeatures(Human.mouse.iDSC, selection.method = "vst", nfeatures = 2000)
#Human.mouse.iDSC <- ScaleData(Human.mouse.iDSC, features = rownames(Human.mouse.iDSC))
Human.mouse.iDSC <- RunPCA(Human.mouse.iDSC, features = VariableFeatures(Human.mouse.iDSC), assay = "SCT")
Human.mouse.iDSC.harmony <- RunHarmony(Human.mouse.iDSC, group.by.vars = "time", assay.use = "SCT")
ElbowPlot(Human.mouse.iDSC.harmony, ndims = 50, reduction = "harmony")
Human.mouse.iDSC.harmony <- RunUMAP(Human.mouse.iDSC.harmony, reduction = "harmony", dims = 1:15)
Human.mouse.iDSC.harmony <- FindNeighbors(Human.mouse.iDSC.harmony, reduction = "harmony", dims = 1:15)
Human.mouse.iDSC.harmony <- FindClusters(Human.mouse.iDSC.harmony, resolution = 0.4)
Human.mouse.iDSC.harmony <- RenameIdents(Human.mouse.iDSC.harmony, "0" = "iDSC_0", "2" = "iDSC_1", "1" = "iDSC_2", "3" = "iDSC_1", "4" = "iDSC_0", "5" = "iDSC_0")
pdf("20220622-Comments/02.human-data/00.figures/09.Human-mouse-harmony-iDSC-integrate-UMAP-color-clusters.pdf", width = 6.5, height = 6)
DimPlot(Human.mouse.iDSC.harmony, pt.size = 1)
dev.off()
literature.color <- c("#FFA07A", "#9370DB", "#3CB371")
pdf("20220622-Comments/02.human-data/00.figures/10.Human-mouse-harmony-iDSC-integrate-UMAP-color-literature.pdf", width = 6.5, height = 6)
DimPlot(Human.mouse.iDSC.harmony, group.by = "literature", pt.size = 1, cols = literature.color)
dev.off()
save(Human.mouse.iDSC.harmony, file = "20220622-Comments/02.human-data/Human-mouse-iDSC-harmony-integration.RData")
DimPlot(Human.mouse.iDSC.harmony, group.by = "cellname")
FeaturePlot(Human.mouse.iDSC.harmony, features = c("Nkg7", "Gzma"))

#######################20221112-integrate two human datasets
load("20220622-Comments/02.human-data/nature/Nature-human-DSC-recluster-CCA.RData")
load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-CCA.RData")
FB.DSC.decidua.combined$literature <- "SD"
FB.DSC.decidua.combined$primary_cluster <- paste("SD_", FB.DSC.decidua.combined$seurat_clusters)
human.S.combined$literature <- "Nature"
human.S.combined$primary_cluster <- paste("Nature_", human.S.combined$seurat_clusters)
set.seed(0)
human.S.combined.sub <- human.S.combined
human.S.combined.sub@active.ident <- factor(human.S.combined.sub$orig.cellname)
human.S.combined.sub <- subset(human.S.combined.sub, subset = nCount_RNA>10000)
human.S.combined.sub <- subset(human.S.combined.sub, downsample = 1000)

features <- SelectIntegrationFeatures(object.list = c(FB.DSC.decidua.combined, human.S.combined.sub))
human.decidua.anchors <- FindIntegrationAnchors(object.list = c(FB.DSC.decidua.combined, human.S.combined.sub), anchor.features = features)
human.decidua.combined.sub <- IntegrateData(anchorset = human.decidua.anchors, 
                                            features.to.integrate = intersect(rownames(FB.DSC.decidua.combined), rownames(human.S.combined.sub)))
DefaultAssay(human.decidua.combined.sub) <- "integrated"
human.decidua.combined.sub <- ScaleData(human.decidua.combined.sub, features = rownames(human.decidua.combined.sub))
human.decidua.combined.sub <- RunPCA(human.decidua.combined.sub, features = VariableFeatures(human.decidua.combined.sub))
ElbowPlot(human.decidua.combined.sub, ndims = 50)
human.decidua.combined.sub <- RunUMAP(human.decidua.combined.sub, reduction = "pca", dims = 1:30)
human.decidua.combined.sub <- FindNeighbors(human.decidua.combined.sub, reduction = "pca", dims = 1:30)
human.decidua.combined.sub <- FindClusters(human.decidua.combined.sub, resolution = 0.6)
human.decidua.combined.sub <- RenameIdents(human.decidua.combined.sub, "0" = "cluster_0", "1" = "cluster_1", "2" = "cluster_2", "3" = "cluster_3", "5" = "cluster_5", 
                                           "4" = "cluster_0", "6" = "cluster_6", "7" = "cluster_7")
DimPlot(human.decidua.combined.sub, group.by = "literature")
pdf("20220622-Comments/02.human-data/00.figures/11.Human-DSC-integrate-UMAP-sub1000-color-cluster.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub, label = F, pt.size = 1)
dev.off()
pdf("20220622-Comments/02.human-data/00.figures/12.Human-DSC-integrate-UMAP-sub1000-color-literature.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub, label = F, pt.size = 1, group.by = "literature", cols = c("#9370DB", "#3CB371"))
dev.off()

FeaturePlot(human.decidua.combined.sub, features = c("nFeature_RNA"))
human.decidua.combined.sub.markers <- FindAllMarkers(human.decidua.combined.sub, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
human.decidua.combined.sub.markers.top10 <- human.decidua.combined.sub.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
pdf("20220622-Comments/02.human-data/00.figures/13.Human-DSC-integrate-sub1000-DEGs-top20-heatmap.pdf", width = 10, height = 10)
DoHeatmap(human.decidua.combined.sub, features = human.decidua.combined.sub.markers.top10$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()
save(human.decidua.combined.sub, file = "20220622-Comments/02.human-data/Human-Nature-SD-DSCs-downsample-500-integrate-20221112.RData")

#average expression of immune-related genes, using human Nature
human.S.meta <- read.table("../../20220508-download-data/Nature-2018/meta_10x.txt", sep = "\t", header = T, check.names = F)
human.S.data <- read.table("../../20220508-download-data/Nature-2018/raw_data_10x.txt", sep = "\t", header = T, check.names = F)
human.S.data.c <- human.S.data
human.S.data.c$Gene <- gsub("_ENSG.*", "", human.S.data.c$Gene)
human.S.data.c <- data.frame(human.S.data.c %>% group_by(Gene) %>% filter(!duplicated(Gene)), check.names = F)
rownames(human.S.data.c) <- human.S.data.c$Gene
human.S.data.c <- human.S.data.c[, -1]
human.S.seurat <- CreateSeuratObject(counts = human.S.data.c, min.cells = 3, min.features = 200, meta.data = human.S.meta)
save(human.S.seurat, file = "../../20220508-download-data/Nature-2018/Human-nature-2018-all-cells-20221113.RData")
human.S.seurat <- NormalizeData(human.S.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
human.S.seurat <- FindVariableFeatures(human.S.seurat, selection.method = "vst", nfeatures = 2000)
human.S.seurat <- ScaleData(human.S.seurat, features = VariableFeatures(human.S.seurat))
human.S.seurat <- RunPCA(human.S.seurat, features = VariableFeatures(human.S.seurat), npcs = 50)
human.S.seurat@active.ident <- factor(human.S.seurat$annotation)
human.S.seurat <- RenameIdents(human.S.seurat, "dS1" = "DSC", "dS2" = "DSC", "dS3" = 'DSC', "dNK2" = "Immune", "Tcells" = "Immune", "dM2" = "Immune", 
                               "dM1" = "Immune", "dNK3" = "Immune", "NK CD16+" = "Immune", "dM3" = "Immune", "MO" = "Immune", "dNK p" = "Immune", 
                               "DC1" = "Immune", "DC2" = "Immune", "Granulocytes" = "Immune", "NK CD16-" = "Immune", "ILC3" = "Immune", "dNK1" = "Immune")
human.DSC.DEGs <- FindMarkers(human.S.seurat, ident.1 = "DSC", only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
human.Immune.DEGs <- FindMarkers(human.S.seurat, ident.1 = "Immune", only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

human.decidua.combined.sub <- AddModuleScore(human.decidua.combined.sub, features = list(rownames(human.DSC.DEGs)), name = "DSC")
human.decidua.combined.sub <- AddModuleScore(human.decidua.combined.sub, features = list(rownames(human.Immune.DEGs)), name = "Immune")
pdf("20220622-Comments/02.human-data/00.figures/14.Human-DSC-module-DSC-immune-featureplot.pdf", width = 20, height = 10)
FeaturePlot(human.decidua.combined.sub, features = c("DSC1", "Immune1"), pt.size = 1) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

#20230410 human DSC gene expression and rename clusters
load("20220622-Comments/02.human-data/Human-Nature-SD-DSCs-downsample-500-integrate.RData")
DimPlot(human.decidua.combined.sub)
human.decidua.combined.sub <- RenameIdents(human.decidua.combined.sub, "0" = "hFB1", "3" = "hFB1", "1" = "hFB2", "2" = "hDSC", "4" = "hiDSC", "5" = "hFB2", "6" = "hFB1")
DimPlot(human.decidua.combined.sub)
VlnPlot(human.decidua.combined.sub, features = "IL1B", pt.size = 0)
pdf("20220622-Comments/02.human-data/00.figures/15.Human-DSC-IL1B-expression-featureplot.pdf", width = 10, height = 10)
FeaturePlot(human.decidua.combined.sub, features = "IL1B", pt.size = 1) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()
pdf("20220622-Comments/02.human-data/00.figures/15.Human-DSC-IL1B-expression-vlnplot.pdf", width = 10, height = 10)
VlnPlot(human.decidua.combined.sub, features = "IL1B", pt.size = 0, cols = c("#57BBC2", "#C0A335", "#006400", "#6A3D9A"))
dev.off()

pdf("20220622-Comments/02.human-data/00.figures/16.Human-DSC-integrate-cluster-dimplot.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub, pt.size = 1, cols = c("#57BBC2", "#C0A335", "#006400", "#6A3D9A"))
dev.off()

pdf("20220622-Comments/02.human-data/00.figures/16.Human-DSC-orig-cluster-dimplot.pdf", width = 10, height = 10)
DimPlot(human.decidua.combined.sub, pt.size = 1, group.by = "orig.cellname")
dev.off()

load("20220622-Comments/02.human-data/nature/Nature-human-DSC-recluster-CCA.RData")
DimPlot(human.S.combined, group.by = "annotation")
human.S.combined <- RenameIdents(human.S.combined, "1" = "hFB1", "0" = "hFB2", "2" = "hFB1", "3" = "hDSC", "4" = "hiDSC", "5" = "hFB2")
pdf("20220622-Comments/02.human-data/00.figures/17.Human-DSC-Nature-cluster-dimplot.pdf", width = 10, height = 10)
DimPlot(human.S.combined, pt.size = 1, cols = c("#57BBC2", "#C0A335", "#006400", "#6A3D9A"))
dev.off()
load("20220622-Comments/02.human-data/science_advance/decidua-subset-FB-CCA.RData")
DimPlot(FB.DSC.decidua.combined, group.by = "orig.cellname")
dev.off()
FB.DSC.decidua.combined <- RenameIdents(FB.DSC.decidua.combined, "0" = "hFB1", "2" = "hFB1", "3" = "hFB2", "4" = "hFB1", "1" = "hDSC", "5" = "hiDSC")
pdf("20220622-Comments/02.human-data/00.figures/17.Human-DSC-SA-cluster-dimplot.pdf", width = 10, height = 10)
DimPlot(FB.DSC.decidua.combined, pt.size = 2, cols = c("#57BBC2", "#C0A335", "#006400", "#6A3D9A"))
dev.off()

load("FB-immune-DSC-harmony-20221218.RData")
DimPlot(UE.FB)
pdf("20220622-Comments/02.human-data/00.figures/15.Mouse-DSC-IL1B-expression-featureplot.pdf", width = 10, height = 10)
FeaturePlot(UE.FB, features = "Il1b", pt.size = 1) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()
pdf("20220622-Comments/02.human-data/00.figures/15.Mouse-DSC-IL1B-expression-vlnplot.pdf", width = 10, height = 10)
VlnPlot(UE.FB, features = "Il1b", pt.size = 0, group.by = "cellname", cols = c(decidual.color[1:7], immune.color[10]))
dev.off()

#chemokine/Mdk/Lgals9 gene expression in human DSC and iDSC
human.DSC.iDSC.avg <- AverageExpression(human.decidua.combined.sub, slot = "data", assays = "RNA", use.counts = T) #use count matrix
human.DSC.iDSC.avg <- human.DSC.iDSC.avg$RNA
human.DSC.iDSC.avg.chemokine <- human.DSC.iDSC.avg[toupper(c(chemokine.ligands[1:7], "CCL27", chemokine.ligands[9:15], "CCL21", chemokine.ligands[17:18], "Mdk", "Lgals9")),]
color <- rev(brewer.pal(11, "Spectral"))
pdf("20220622-Comments/02.human-data/00.figures/18.human-DSC-iDSC-chemokines-MK-GALECTIN-L-heatmap.pdf", width = 10, height = 10)
pheatmap(log2(human.DSC.iDSC.avg.chemokine+1), scale = "none", cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(colors = color)(100))
dev.off()

#vascularization
vascular.gene <- c("CALCA", "ADM", "ANGPT2", "ANGPTL2", "ANGPTL4", "RARRES2", "SEMA3A", "SEMA3B", "SEMA3C", "SEMA3E", "VEGFA", "VEGFB", "VEGFD")
pdf("20220622-Comments/02.human-data/00.figures/19.Human-DSC-iDSC-Vascularization-expression-vlnplot.pdf", width = 10, height = 10)
VlnPlot(human.decidua.combined.sub, features = vascular.gene, pt.size = 0, cols = c("#57BBC2", "#C0A335", "#006400", "#6A3D9A"))
dev.off()

#human DSC/iDSC DEGs
DefaultAssay(human.decidua.combined.sub) <- "RNA"
human.decidua.combined.sub <- ScaleData(human.decidua.combined.sub, features = rownames(human.decidua.combined.sub))
human.decidua.combined.sub.markers <- FindAllMarkers(human.decidua.combined.sub, min.pct = 0.5, logfc.threshold = 0.25)
human.decidua.combined.sub.markers <- human.decidua.combined.sub.markers[human.decidua.combined.sub.markers$gene != "MALAT1", ]
human.decidua.combined.sub.markers.top20 <- human.decidua.combined.sub.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
pdf("20220622-Comments/02.human-data/00.figures/20.Human-DSC-subcluster-DEGs-top10.pdf", width = 10, height = 10)
DoHeatmap(human.decidua.combined.sub, features = human.decidua.combined.sub.markers.top20$gene)+ 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()
write.table(human.decidua.combined.sub.markers, file = "20220622-Comments/02.human-data/Human-DSC-subclusters-FB1-FB2-DSC-DEGs-fc0.25.csv", sep = ",", col.names = T, row.names = F, quote = F)

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
  enrich.GO <- enrich.GO@result[enrich.GO@result$pvalue<0.01, ]
  return(enrich.GO)
}
human.iDSC.GO.df <- data.frame()
for (c in unique(human.decidua.combined.sub.markers$cluster)) {
  print(c)
  geneset <- human.decidua.combined.sub.markers[human.decidua.combined.sub.markers$cluster == c, ]$gene
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
human.iDSC.GO.df.top5 <- human.iDSC.GO.df %>% group_by(cluster) %>% top_n(n = 5, -log(pvalue))
pdf("20220622-Comments/05.DBA/00.figures/06.DBA-deciuda-human.iDSC-DEGs-GO-terms-barplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(human.iDSC.GO.df.top5, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()
write.table(human.iDSC.GO.df, file = "20220622-Comments/02.human-data/Human-DSC-subclusters-FB1-FB2-DSC-GOs.csv", sep = ",", col.names = T, row.names = F, quote = F)

