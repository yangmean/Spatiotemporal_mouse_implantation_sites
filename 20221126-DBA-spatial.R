library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(future)
library(tidyr)

#prepare for spatial integration
load("20210827-figures/files/DBA.1-DSC-20221212.Rdata")
DimPlot(DBA.1.DSC, group.by = "cellcluster")
load("20210827-figures/files/UE-DSC-web.RData")
DimPlot(UE.decidual.susbet)
UE.decidual.susbet$cellcluster <- UE.decidual.susbet@active.ident
DSC.D1.D2 <- subset(UE.decidual.susbet, subset = cellcluster %in% c("D1_eSF", "D2_Pre-DSC"))
load("20220622-Comments/05.DBA/DBA-EC-clusters-20221215.RData")
DimPlot(DBA.1.EC, group.by = "cellcluster")
load("20210827-figures/files/DBA.1-immune-20221215.RData")
DimPlot(DBA.1.immune.m, group.by = "cellcluster")
load("../20210201-integrate/UE-trophoblast-hamrony.RData")
UE.tropho$cellcluster <- UE.tropho@active.ident
DimPlot(UE.tropho, group.by = "cellcluster")
load("../20210201-integrate/UE-harmony-pc30-UMAP-rename.RData")
DimPlot(UE.harmony.pc30)
UE.harmony.pc30$cellcluster <- UE.harmony.pc30@active.ident
UE.fetal <- subset(UE.harmony.pc30, subset = cellcluster %in% c("Mesenchymal stem cell", "Mesenchymal-epithelial transition", "Blood P", 
                                                                "Viseral endoderm", "Allantois mesodermal", "Fetal EC", "Erythroid"))
load("20210827-figures/files/DBA.1-20220912.Rdata")
DBA.1$cellcluster <- DBA.1@active.ident
DBA.1.maternal <- subset(DBA.1, subset = cellcluster %in% c("SMC", "Epithelial"))
DBA.1.ST.inte <- merge(DBA.1.maternal, y = list(DSC.D1.D2, UE.fetal, UE.tropho, DBA.1.DSC, DBA.1.EC, DBA.1.immune.m))
DBA.1.ST.inte <- SCTransform(DBA.1.ST.inte, ncells = 3000, verbose = FALSE)
DBA.1.ST.inte <- RunPCA(DBA.1.ST.inte)
DBA.1.ST.inte <- RunUMAP(DBA.1.ST.inte, dims = 1:30)
save(DBA.1.ST.inte, file = "20220622-Commentsx/05.DBA/DBA.1-merge-embryo-SCT.RData")

load("20220622-Comments/05.DBA/DBA.1-merge-embryo-SCT.RData")
DimPlot(DBA.1.ST.inte, group.by = "cellcluster")
cellsubset <- c("D1_eSF", "D2_Pre-DSC", "D3_Top2a", "D4_Ptn", "D5_Gatm", "D6_S100a8", "D7_Prl8a2", 
                "SMC", "Epithelial", "Monocyte", "Mac", "NK", "Neutrophil", "Proliferating Mac", 
                "DC-1", "T cell", "Proliterating EC", "Angiogenic EC", "EPC", "Venous EC1", "Venous EC2", 
                "Arterial EC", "tropho-0", "tropho-1", "tropho-2", "tropho-3", "tropho-4", "tropho-5", 
                "Viseral endoderm")

DBA.1.ST.inte.embryo <- subset(DBA.1.ST.inte, cellcluster %in% c("tropho-0", "tropho-1", "tropho-2", "tropho-3", "tropho-4", "tropho-5", 
                                                                 "Viseral endoderm", "Mesenchymal-epithelial transition", "Mesenchymal stem cell", 
                                                                 "Fetal EC", "Blood P", "Allantois mesodermal"))
DBA.1.ST.inte.maternal <- subset(DBA.1.ST.inte, cellcluster %in% c("D1_eSF", "D2_Pre-DSC", "D3_Top2a", "D4_Ptn", "D5_Gatm", "D6_S100a8", "D7_Prl8a2", 
                                                                   "SMC", "Epithelial", "Monocyte", "Mac", "NK", "Neutrophil", "Proliferating Mac", 
                                                                   "DC-1", "T cell", "Proliterating EC", "Angiogenic EC", "EPC", "Venous EC1", "Venous EC2", 
                                                                   "Arterial EC", "iDSC_DBA1"))
DBA.E3.data <- read.table("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/05.DBA/SS200000910TR_E3.TissueCut.gem.gz", sep = "\t", header = T, check.names = F)
DBA.E3.data.bin50 <- DBA.E3.data
DBA.E3.data.bin50$x <- trunc(DBA.E3.data.bin50$x/50) * 50
DBA.E3.data.bin50$y <- trunc(DBA.E3.data.bin50$y/50) * 50
DBA.E3.data.bin50$cellID <- paste(DBA.E3.data.bin50$x, "_", DBA.E3.data.bin50$y, sep = "")
DBA.E3.data.bin50 <- aggregate(DBA.E3.data.bin50$MIDCount, by = list(DBA.E3.data.bin50$cellID, DBA.E3.data.bin50$geneID), sum)
colnames(DBA.E3.data.bin50) <- c("cellID", "geneID", "MIDCount")
DBA.E3.data.bin50$cellInx <- match(DBA.E3.data.bin50$cellID, unique(DBA.E3.data.bin50$cellID))
DBA.E3.data.bin50$cellInx <- match(DBA.E3.data.bin50$cellID, unique(DBA.E3.data.bin50$cellID))
DBA.E3.data.bin50$geneInx <- match(DBA.E3.data.bin50$geneID, unique(DBA.E3.data.bin50$geneID))
mat <- sparseMatrix(i = DBA.E3.data.bin50$geneInx, j = DBA.E3.data.bin50$cellInx, x = DBA.E3.data.bin50$MIDCount, 
                    dimnames = list(unique(DBA.E3.data.bin50$geneID), unique(DBA.E3.data.bin50$cellID)))
DBA.E3.bin50.coord.df <- data.frame(cellname = colnames(mat))
rownames(DBA.E3.bin50.coord.df) <- DBA.E3.bin50.coord.df$cellname
DBA.E3.bin50.coord.df <- separate(DBA.E3.bin50.coord.df, col = cellname, sep = "_", into = c("x", "y"))
DBA.E3.bin50 <- CreateSeuratObject(mat, project = "DBA.E3.bin50", assay = "Spatial")
DBA.E3.bin50$slice <- 1
DBA.E3.bin50$region <- "DBA.E3.bin50"
colnames(DBA.E3.bin50.coord.df) <- c("imagerow", "imagecol")
DBA.E3.bin50.coord.df$imagerow <- as.numeric(DBA.E3.bin50.coord.df$imagerow)
DBA.E3.bin50.coord.df$imagecol <- as.numeric(DBA.E3.bin50.coord.df$imagecol)
DBA.E3.bin50@images$DBA.E3_bin50 <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = DBA.E3.bin50.coord.df)
#SCT transform
DBA.E3.bin50 <- SCTransform(DBA.E3.bin50, assay = "Spatial", verbose = FALSE)
DBA.E3.bin50 <- RunPCA(DBA.E3.bin50, assay = "SCT", verbose = F)
ElbowPlot(DBA.E3.bin50, ndims = 50)
DBA.E3.bin50 <- FindNeighbors(DBA.E3.bin50, reduction = "pca", dims = 1:30)
DBA.E3.bin50 <- FindClusters(DBA.E3.bin50, verbose = T, resolution = 0.3)
DBA.E3.bin50 <- RunUMAP(DBA.E3.bin50, reduction = "pca", dims = 1:30)
DimPlot(DBA.E3.bin50, label = T)
SpatialDimPlot(DBA.E3.bin50, label = T)
save(DBA.E3.bin50, file = "20220622-Comments/05.DBA/DBA-E3-data.Rdata")
rm(DBA.E3.bin50)
gc()
DBA.E3.bin50.markers <- FindAllMarkers(DBA.E3.bin50, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.E3.bin50.markers.top20 <- DBA.E3.bin50.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
write.table(DBA.E3.bin50.markers.top20, file = "20220622-Comments/05.DBA/DBA-E3-spatial-cluster-DEGs-top20.csv", sep = ",", quote = F)
#transfer E3.bin50 and single-cell RNA-seq
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte, query = DBA.E3.bin50, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E3.bin50[["pca"]], dims = 1:30)
DBA.E3.bin50[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E3.bin50) <- "predictions"
save(DBA.E3.bin50, file = "20220622-Comments/05.DBA/DBA-E3-data-map-20221212.RData")

#选择每个细胞分数最高的为cell type
DefaultAssay(DBA.E3.bin50) <-  "predictions"
DBA.E3.bin50.celltype <- GetAssayData(DBA.E3.bin50, slot = "data")
DBA.E3.bin50.celltype <- apply(DBA.E3.bin50.celltype, 2, function(x) rownames(DBA.E3.bin50.celltype)[which.max(x)])
DBA.E3.bin50.celltype <- data.frame(cell_type = DBA.E3.bin50.celltype, 
                                cell_name = names(DBA.E3.bin50.celltype), 
                                GetTissueCoordinates(DBA.E3.bin50), 
                                spatial_clusters = DBA.E3.bin50$seurat_clusters,
                                nFeature_RNA = DBA.E3.bin50$nFeature_Spatial, 
                                nCount_RNA = DBA.E3.bin50$nCount_Spatial)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#E6
DBA.E6.data <- read.table("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/05.DBA/SS200000910TR_E6.TissueCut.gem.gz", sep = "\t", header = T, check.names = F)
DBA.E6.data.bin50 <- DBA.E6.data
DBA.E6.data.bin50$x <- trunc(DBA.E6.data.bin50$x/50) * 50
DBA.E6.data.bin50$y <- trunc(DBA.E6.data.bin50$y/50) * 50
DBA.E6.data.bin50$cellID <- paste(DBA.E6.data.bin50$x, "_", DBA.E6.data.bin50$y, sep = "")
DBA.E6.data.bin50 <- aggregate(DBA.E6.data.bin50$MIDCount, by = list(DBA.E6.data.bin50$cellID, DBA.E6.data.bin50$geneID), sum)
colnames(DBA.E6.data.bin50) <- c("cellID", "geneID", "MIDCount")
DBA.E6.data.bin50$cellInx <- match(DBA.E6.data.bin50$cellID, unique(DBA.E6.data.bin50$cellID))
DBA.E6.data.bin50$cellInx <- match(DBA.E6.data.bin50$cellID, unique(DBA.E6.data.bin50$cellID))
DBA.E6.data.bin50$geneInx <- match(DBA.E6.data.bin50$geneID, unique(DBA.E6.data.bin50$geneID))
mat <- sparseMatrix(i = DBA.E6.data.bin50$geneInx, j = DBA.E6.data.bin50$cellInx, x = DBA.E6.data.bin50$MIDCount, 
                    dimnames = list(unique(DBA.E6.data.bin50$geneID), unique(DBA.E6.data.bin50$cellID)))
DBA.E6.bin50.coord.df <- data.frame(cellname = colnames(mat))
rownames(DBA.E6.bin50.coord.df) <- DBA.E6.bin50.coord.df$cellname
DBA.E6.bin50.coord.df <- separate(DBA.E6.bin50.coord.df, col = cellname, sep = "_", into = c("x", "y"))
DBA.E6.bin50 <- CreateSeuratObject(mat, project = "DBA.E6.bin50", assay = "Spatial")
DBA.E6.bin50$slice <- 1
DBA.E6.bin50$region <- "DBA.E6.bin50"
colnames(DBA.E6.bin50.coord.df) <- c("imagerow", "imagecol")
DBA.E6.bin50.coord.df$imagerow <- as.numeric(DBA.E6.bin50.coord.df$imagerow)
DBA.E6.bin50.coord.df$imagecol <- as.numeric(DBA.E6.bin50.coord.df$imagecol)
DBA.E6.bin50@images$DBA.E6_bin50 <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = DBA.E6.bin50.coord.df)
#SCT transform
DBA.E6.bin50 <- SCTransform(DBA.E6.bin50, assay = "Spatial", verbose = FALSE)
DBA.E6.bin50 <- RunPCA(DBA.E6.bin50, assay = "SCT", verbose = F)
ElbowPlot(DBA.E6.bin50, ndims = 50)
DBA.E6.bin50 <- FindNeighbors(DBA.E6.bin50, reduction = "pca", dims = 1:30)
DBA.E6.bin50 <- FindClusters(DBA.E6.bin50, verbose = T, resolution = 0.3)
DBA.E6.bin50 <- RunUMAP(DBA.E6.bin50, reduction = "pca", dims = 1:30)
DimPlot(DBA.E6.bin50, label = T)
SpatialDimPlot(DBA.E6.bin50, label = T)
save(DBA.E6.bin50, file = "20220622-Comments/05.DBA/DBA-E6-data.Rdata")
rm(DBA.E6.bin50)
gc()

#transfer E6.bin50 and single-cell RNA-seq
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte, query = DBA.E6.bin50, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E6.bin50[["pca"]], dims = 1:30)
DBA.E6.bin50[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E6.bin50) <- "predictions"
save(DBA.E6.bin50, file = "20220622-Comments/05.DBA/DBA-E6-data-map-20221212.RData")

#选择每个细胞分数最高的为cell type
DefaultAssay(DBA.E6.bin50) <-  "predictions"
DBA.E6.bin50.celltype <- GetAssayData(DBA.E6.bin50, slot = "data")
DBA.E6.bin50.celltype <- apply(DBA.E6.bin50.celltype, 2, function(x) rownames(DBA.E6.bin50.celltype)[which.max(x)])
DBA.E6.bin50.celltype <- data.frame(cell_type = DBA.E6.bin50.celltype, 
                                    cell_name = names(DBA.E6.bin50.celltype), 
                                    GetTissueCoordinates(DBA.E6.bin50), 
                                    spatial_clusters = DBA.E6.bin50$seurat_clusters,
                                    nFeature_RNA = DBA.E6.bin50$nFeature_Spatial, 
                                    nCount_RNA = DBA.E6.bin50$nCount_Spatial)
ggplot() + 
  geom_point(DBA.E6.bin50.celltype, mapping = aes(x = y, y = -x, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#F1
DBA.F1.data <- read.table("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/05.DBA/SS200000910TR_F1.TissueCut.gem.gz", sep = "\t", header = T, check.names = F)
DBA.F1.data.bin50 <- DBA.F1.data
DBA.F1.data.bin50$x <- trunc(DBA.F1.data.bin50$x/50) * 50
DBA.F1.data.bin50$y <- trunc(DBA.F1.data.bin50$y/50) * 50
DBA.F1.data.bin50$cellID <- paste(DBA.F1.data.bin50$x, "_", DBA.F1.data.bin50$y, sep = "")
DBA.F1.data.bin50 <- aggregate(DBA.F1.data.bin50$MIDCount, by = list(DBA.F1.data.bin50$cellID, DBA.F1.data.bin50$geneID), sum)
colnames(DBA.F1.data.bin50) <- c("cellID", "geneID", "MIDCount")
DBA.F1.data.bin50$cellInx <- match(DBA.F1.data.bin50$cellID, unique(DBA.F1.data.bin50$cellID))
DBA.F1.data.bin50$cellInx <- match(DBA.F1.data.bin50$cellID, unique(DBA.F1.data.bin50$cellID))
DBA.F1.data.bin50$geneInx <- match(DBA.F1.data.bin50$geneID, unique(DBA.F1.data.bin50$geneID))
mat <- sparseMatrix(i = DBA.F1.data.bin50$geneInx, j = DBA.F1.data.bin50$cellInx, x = DBA.F1.data.bin50$MIDCount, 
                    dimnames = list(unique(DBA.F1.data.bin50$geneID), unique(DBA.F1.data.bin50$cellID)))
DBA.F1.bin50.coord.df <- data.frame(cellname = colnames(mat))
rownames(DBA.F1.bin50.coord.df) <- DBA.F1.bin50.coord.df$cellname
DBA.F1.bin50.coord.df <- separate(DBA.F1.bin50.coord.df, col = cellname, sep = "_", into = c("x", "y"))
DBA.F1.bin50 <- CreateSeuratObject(mat, project = "DBA.F1.bin50", assay = "Spatial")
DBA.F1.bin50$slice <- 1
DBA.F1.bin50$region <- "DBA.F1.bin50"
colnames(DBA.F1.bin50.coord.df) <- c("imagerow", "imagecol")
DBA.F1.bin50.coord.df$imagerow <- as.numeric(DBA.F1.bin50.coord.df$imagerow)
DBA.F1.bin50.coord.df$imagecol <- as.numeric(DBA.F1.bin50.coord.df$imagecol)
DBA.F1.bin50@images$DBA.F1_bin50 <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = DBA.F1.bin50.coord.df)
#SCT transform
DBA.F1.bin50 <- SCTransform(DBA.F1.bin50, assay = "Spatial", verbose = FALSE)
DBA.F1.bin50 <- RunPCA(DBA.F1.bin50, assay = "SCT", verbose = F)
ElbowPlot(DBA.F1.bin50, ndims = 50)
DBA.F1.bin50 <- FindNeighbors(DBA.F1.bin50, reduction = "pca", dims = 1:30)
DBA.F1.bin50 <- FindClusters(DBA.F1.bin50, verbose = T, resolution = 0.3)
DBA.F1.bin50 <- RunUMAP(DBA.F1.bin50, reduction = "pca", dims = 1:30)
DimPlot(DBA.F1.bin50, label = T)
SpatialDimPlot(DBA.F1.bin50, label = T)
save(DBA.F1.bin50, file = "20220622-Comments/05.DBA/DBA-F1-data.Rdata")
rm(DBA.F1.bin50)
gc()

#transfer F1.bin50 and single-cell RNA-seq
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte, query = DBA.F1.bin50, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.F1.bin50[["pca"]], dims = 1:30)
DBA.F1.bin50[["predictions"]] <- predictions.assay
DefaultAssay(DBA.F1.bin50) <- "predictions"
save(DBA.F1.bin50, file = "20220622-Comments/05.DBA/DBA-F1-data-map-20221212.RData")

#选择每个细胞分数最高的为cell type
DefaultAssay(DBA.F1.bin50) <-  "predictions"
DBA.F1.bin50.celltype <- GetAssayData(DBA.F1.bin50, slot = "data")
DBA.F1.bin50.celltype <- apply(DBA.F1.bin50.celltype, 2, function(x) rownames(DBA.F1.bin50.celltype)[which.max(x)])
DBA.F1.bin50.celltype <- data.frame(cell_type = DBA.F1.bin50.celltype, 
                                    cell_name = names(DBA.F1.bin50.celltype), 
                                    GetTissueCoordinates(DBA.F1.bin50), 
                                    spatial_clusters = DBA.F1.bin50$seurat_clusters,
                                    nFeature_RNA = DBA.F1.bin50$nFeature_Spatial, 
                                    nCount_RNA = DBA.F1.bin50$nCount_Spatial)
ggplot() + 
  geom_point(DBA.F1.bin50.celltype, mapping = aes(x = -x, y = -y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#20221204 load normal mouse clusters
celltype.level <- c("Allantois mesodermal", "Blood P", "Lum-D", "Sfrp4-D", "Top2a-D", "Ptn-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "DC-1", "DC-2", "EC-0", "EC-1", 
                    "EC-2", "EC-3", "EC-4", "EC-5", "Epithelial", "Erythroid", "Fetal EC", "FB-immune", "Mac", "Mesenchymal stem cell", "Mesenchymal-epithelial transition", 
                    "Mono", "Neutrophil", "NK", "Pro-Mac", "SMC", "T", "tropho-0", "tropho-1", "tropho-2", "tropho-3", "tropho-4", "tropho-5", "Viseral endoderm", 
                    "cluster9-1", "B")
load("20210818-spatial/UE-harmony-pc30-UMAP-subname-SCT.RData")
load("20220622-Comments/05.DBA/DBA-E3-data.Rdata")
#E3
#transfer E3.bin50 and single-cell RNA-seq
anchors <- FindTransferAnchors(reference = UE.harmony.SCT, query = DBA.E3.bin50, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = UE.harmony.SCT$subname, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E3.bin50[["pca"]], dims = 1:30)
DBA.E3.bin50[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E3.bin50) <- "predictions"
save(DBA.E3.bin50, file = "20220622-Comments/05.DBA/DBA-E3-data-normal-sc-map.RData")

#选择每个细胞分数最高的为cell type
DefaultAssay(DBA.E3.bin50) <-  "predictions"
DBA.E3.bin50.celltype <- GetAssayData(DBA.E3.bin50, slot = "data")
DBA.E3.bin50.celltype <- apply(DBA.E3.bin50.celltype, 2, function(x) rownames(DBA.E3.bin50.celltype)[which.max(x)])
DBA.E3.bin50.celltype <- data.frame(cell_type = DBA.E3.bin50.celltype, 
                                    cell_name = names(DBA.E3.bin50.celltype), 
                                    GetTissueCoordinates(DBA.E3.bin50), 
                                    spatial_clusters = DBA.E3.bin50$seurat_clusters,
                                    nFeature_RNA = DBA.E3.bin50$nFeature_Spatial, 
                                    nCount_RNA = DBA.E3.bin50$nCount_Spatial)
DBA.E3.bin50.celltype$cell_type <- factor(DBA.E3.bin50.celltype$cell_type, levels = celltype.level)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#E6
load("20220622-Comments/05.DBA/DBA-E6-data.Rdata")
#transfer E6.bin50 and single-cell RNA-seq
anchors <- FindTransferAnchors(reference = UE.harmony.SCT, query = DBA.E6.bin50, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = UE.harmony.SCT$subname, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E6.bin50[["pca"]], dims = 1:30)
DBA.E6.bin50[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E6.bin50) <- "predictions"
save(DBA.E6.bin50, file = "20220622-Comments/05.DBA/DBA-E6-data-normal-sc-map.RData")

#选择每个细胞分数最高的为cell type
DefaultAssay(DBA.E6.bin50) <-  "predictions"
DBA.E6.bin50.celltype <- GetAssayData(DBA.E6.bin50, slot = "data")
DBA.E6.bin50.celltype <- apply(DBA.E6.bin50.celltype, 2, function(x) rownames(DBA.E6.bin50.celltype)[which.max(x)])
DBA.E6.bin50.celltype <- data.frame(cell_type = DBA.E6.bin50.celltype, 
                                    cell_name = names(DBA.E6.bin50.celltype), 
                                    GetTissueCoordinates(DBA.E6.bin50), 
                                    spatial_clusters = DBA.E6.bin50$seurat_clusters,
                                    nFeature_RNA = DBA.E6.bin50$nFeature_Spatial, 
                                    nCount_RNA = DBA.E6.bin50$nCount_Spatial)
DBA.E6.bin50.celltype$cell_type <- factor(DBA.E6.bin50.celltype$cell_type, levels = celltype.level)
ggplot() + 
  geom_point(DBA.E6.bin50.celltype, mapping = aes(x = y, y = -x, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#F1
load("20220622-Comments/05.DBA/DBA-F1-data.Rdata")
#transfer F1.bin50 and single-cell RNA-seq
anchors <- FindTransferAnchors(reference = UE.harmony.SCT, query = DBA.F1.bin50, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = UE.harmony.SCT$subname, 
                                  prediction.assay = TRUE, weight.reduction = DBA.F1.bin50[["pca"]], dims = 1:30)
DBA.F1.bin50[["predictions"]] <- predictions.assay
DefaultAssay(DBA.F1.bin50) <- "predictions"
save(DBA.F1.bin50, file = "20220622-Comments/05.DBA/DBA-F1-data-normal-sc-map.RData")

#选择每个细胞分数最高的为cell type
DefaultAssay(DBA.F1.bin50) <-  "predictions"
DBA.F1.bin50.celltype <- GetAssayData(DBA.F1.bin50, slot = "data")
DBA.F1.bin50.celltype <- apply(DBA.F1.bin50.celltype, 2, function(x) rownames(DBA.F1.bin50.celltype)[which.max(x)])
DBA.F1.bin50.celltype <- data.frame(cell_type = DBA.F1.bin50.celltype, 
                                    cell_name = names(DBA.F1.bin50.celltype), 
                                    GetTissueCoordinates(DBA.F1.bin50), 
                                    spatial_clusters = DBA.F1.bin50$seurat_clusters,
                                    nFeature_RNA = DBA.F1.bin50$nFeature_Spatial, 
                                    nCount_RNA = DBA.F1.bin50$nCount_Spatial)
DBA.F1.bin50.celltype$cell_type <- factor(DBA.F1.bin50.celltype$cell_type, levels = celltype.level)
ggplot() + 
  geom_point(DBA.F1.bin50.celltype, mapping = aes(x = x, y = -y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#20221207 merge normal and DBA single cell data clusters and integrate to DBA spatial 
load("20210818-spatial/UE-harmony-pc30-UMAP-subname-SCT.RData")
UE.harmony.SCT$cellcluster <- paste("N_", UE.harmony.SCT$subname, sep = "")
load("20220622-Comments/05.DBA/DBA.1-merge-embryo-SCT.RData")
DBA.1.ST.inte$cellcluster <- paste("DBA_", DBA.1.ST.inte$cellcluster, sep = "")
UE.Normal.DBA <- merge(UE.harmony.SCT, DBA.1.ST.inte)
UE.Normal.DBA$sample <- gsub("_.*", "", UE.Normal.DBA$cellcluster)
UE.Normal.DBA <- SCTransform(UE.Normal.DBA, ncells = 3000, verbose = FALSE)
UE.Normal.DBA@active.ident <- factor(UE.Normal.DBA$cellcluster)
UE.Normal.DBA <- RunPCA(UE.Normal.DBA)
UE.Normal.DBA <- RunUMAP(UE.Normal.DBA, dims = 1:30)

#E3
load("20220622-Comments/05.DBA/DBA-E3-data.Rdata")
#transfer E3.bin50 and single-cell RNA-seq
anchors <- FindTransferAnchors(reference = UE.Normal.DBA, query = DBA.E3.bin50, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = UE.Normal.DBA$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E3.bin50[["pca"]], dims = 1:30)
DBA.E3.bin50[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E3.bin50) <- "predictions"
save(DBA.E3.bin50, file = "20220622-Comments/05.DBA/DBA-E3-data-map-NormalandDBA.RData")
#选择每个细胞分数最高的为cell type
DefaultAssay(DBA.E3.bin50) <-  "predictions"
DBA.E3.bin50.celltype <- GetAssayData(DBA.E3.bin50, slot = "data")
DBA.E3.bin50.celltype <- apply(DBA.E3.bin50.celltype, 2, function(x) rownames(DBA.E3.bin50.celltype)[which.max(x)])
DBA.E3.bin50.celltype <- data.frame(cell_type = DBA.E3.bin50.celltype, 
                                    cell_name = names(DBA.E3.bin50.celltype), 
                                    GetTissueCoordinates(DBA.E3.bin50), 
                                    spatial_clusters = DBA.E3.bin50$seurat_clusters,
                                    nFeature_RNA = DBA.E3.bin50$nFeature_Spatial, 
                                    nCount_RNA = DBA.E3.bin50$nCount_Spatial)
DBA.E3.bin50.celltype$cell_type <- factor(DBA.E3.bin50.celltype$cell_type, levels = celltype.level)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

load("20220622-Comments/05.DBA/DBA-E6-data.Rdata")
#transfer E6.bin50 and single-cell RNA-seq
anchors <- FindTransferAnchors(reference = UE.Normal.DBA, query = DBA.E6.bin50, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = UE.Normal.DBA$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E6.bin50[["pca"]], dims = 1:30)
DBA.E6.bin50[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E6.bin50) <- "predictions"
save(DBA.E6.bin50, file = "20220622-Comments/05.DBA/DBA-E6-data-map-NormalandDBA.RData")
#选择每个细胞分数最高的为cell type
DefaultAssay(DBA.E6.bin50) <-  "predictions"
DBA.E6.bin50.celltype <- GetAssayData(DBA.E6.bin50, slot = "data")
DBA.E6.bin50.celltype <- apply(DBA.E6.bin50.celltype, 2, function(x) rownames(DBA.E6.bin50.celltype)[which.max(x)])
DBA.E6.bin50.celltype <- data.frame(cell_type = DBA.E6.bin50.celltype, 
                                    cell_name = names(DBA.E6.bin50.celltype), 
                                    GetTissueCoordinates(DBA.E6.bin50), 
                                    spatial_clusters = DBA.E6.bin50$seurat_clusters,
                                    nFeature_RNA = DBA.E6.bin50$nFeature_Spatial, 
                                    nCount_RNA = DBA.E6.bin50$nCount_Spatial)
DBA.E6.bin50.celltype$cell_type <- factor(DBA.E6.bin50.celltype$cell_type, levels = celltype.level)
ggplot() + 
  geom_point(DBA.E6.bin50.celltype, mapping = aes(x = y, y = -x, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

load("20220622-Comments/05.DBA/DBA-F1-data.Rdata")
#transfer F1.bin50 and single-cell RNA-seq
anchors <- FindTransferAnchors(reference = UE.Normal.DBA, query = DBA.F1.bin50, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = UE.Normal.DBA$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.F1.bin50[["pca"]], dims = 1:30)
DBA.F1.bin50[["predictions"]] <- predictions.assay
DefaultAssay(DBA.F1.bin50) <- "predictions"
save(DBA.F1.bin50, file = "20220622-Comments/05.DBA/DBA-F1-data-map-NormalandDBA.RData")

DefaultAssay(DBA.F1.bin50) <-  "predictions"
DBA.F1.bin50.celltype <- GetAssayData(DBA.F1.bin50, slot = "data")
DBA.F1.bin50.celltype <- apply(DBA.F1.bin50.celltype, 2, function(x) rownames(DBA.F1.bin50.celltype)[which.max(x)])
DBA.F1.bin50.celltype <- data.frame(cell_type = DBA.F1.bin50.celltype, 
                                    cell_name = names(DBA.F1.bin50.celltype), 
                                    GetTissueCoordinates(DBA.F1.bin50), 
                                    spatial_clusters = DBA.F1.bin50$seurat_clusters,
                                    nFeature_RNA = DBA.F1.bin50$nFeature_Spatial, 
                                    nCount_RNA = DBA.F1.bin50$nCount_Spatial)
DBA.F1.bin50.celltype$cell_type <- factor(DBA.F1.bin50.celltype$cell_type, levels = celltype.level)
ggplot() + 
  geom_point(DBA.F1.bin50.celltype, mapping = aes(x = x, y = -y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#iDSC subcluster integration
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221212.RData")
DimPlot(DBA.1.iDSC, label = T)
#E3
load("20220622-Comments/05.DBA/DBA-E3-data-map-20221212.RData")
DBA.E3.bin50.celltype <- GetAssayData(DBA.E3.bin50, slot = "data")
DBA.E3.bin50.celltype <- apply(DBA.E3.bin50.celltype, 2, function(x) rownames(DBA.E3.bin50.celltype)[which.max(x)])
DBA.E3.bin50$predictions_cell <- DBA.E3.bin50.celltype
DBA.E3.bin50.iDSC <- subset(DBA.E3.bin50, subset = predictions_cell %in% c("iDSC-DBA1"))
DefaultAssay(DBA.E3.bin50.iDSC) <- "SCT"
DBA.E3.bin50.iDSC <- SCTransform(DBA.E3.bin50.iDSC, assay = "Spatial", verbose = FALSE)
DBA.E3.bin50.iDSC <- RunPCA(DBA.E3.bin50.iDSC, assay = "SCT", verbose = F)
ElbowPlot(DBA.E3.bin50.iDSC, ndims = 50)
DBA.E3.bin50.iDSC <- FindNeighbors(DBA.E3.bin50.iDSC, reduction = "pca", dims = 1:30)
DBA.E3.bin50.iDSC <- FindClusters(DBA.E3.bin50.iDSC, verbose = T, resolution = 0.3)
DBA.E3.bin50.iDSC <- RunUMAP(DBA.E3.bin50.iDSC, reduction = "pca", dims = 1:30)
anchors <- FindTransferAnchors(reference = DBA.1.iDSC, query = DBA.E3.bin50.iDSC, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.iDSC$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E3.bin50.iDSC[["pca"]], dims = 1:30)
DBA.E3.bin50.iDSC[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E3.bin50.iDSC) <- "predictions"
save(DBA.E3.bin50.iDSC, file = "20220622-Comments/05.DBA/DBA-E3-iDSC-subclusters-spatial.RData")
DBA.E3.bin50.iDSC.celltype <- GetAssayData(DBA.E3.bin50.iDSC, slot = "data")
DBA.E3.bin50.iDSC.celltype <- apply(DBA.E3.bin50.iDSC.celltype, 2, function(x) rownames(DBA.E3.bin50.iDSC.celltype)[which.max(x)])
DBA.E3.bin50.iDSC.celltype <- data.frame(cell_type = DBA.E3.bin50.iDSC.celltype, 
                                    cell_name = names(DBA.E3.bin50.iDSC.celltype), 
                                    GetTissueCoordinates(DBA.E3.bin50.iDSC), 
                                    spatial_clusters = DBA.E3.bin50.iDSC$seurat_clusters,
                                    nFeature_RNA = DBA.E3.bin50.iDSC$nFeature_Spatial, 
                                    nCount_RNA = DBA.E3.bin50.iDSC$nCount_Spatial)
DBA.E3.bin50.iDSC.celltype$cell_type <- factor(DBA.E3.bin50.iDSC.celltype$cell_type, levels = celltype.level)
ggplot() + 
  geom_point(DBA.E3.bin50.iDSC.celltype, mapping = aes(x = x, y = y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#E6
load("20220622-Comments/05.DBA/DBA-E6-data-map-20221212.RData")
DBA.E6.bin50.celltype <- GetAssayData(DBA.E6.bin50, slot = "data")
DBA.E6.bin50.celltype <- apply(DBA.E6.bin50.celltype, 2, function(x) rownames(DBA.E6.bin50.celltype)[which.max(x)])
DBA.E6.bin50$predictions_cell <- DBA.E6.bin50.celltype
DBA.E6.bin50.iDSC <- subset(DBA.E6.bin50, subset = predictions_cell %in% c("iDSC-DBA1"))
DefaultAssay(DBA.E6.bin50.iDSC) <- "SCT"
DBA.E6.bin50.iDSC <- SCTransform(DBA.E6.bin50.iDSC, assay = "Spatial", verbose = FALSE)
DBA.E6.bin50.iDSC <- RunPCA(DBA.E6.bin50.iDSC, assay = "SCT", verbose = F)
ElbowPlot(DBA.E6.bin50.iDSC, ndims = 50)
DBA.E6.bin50.iDSC <- FindNeighbors(DBA.E6.bin50.iDSC, reduction = "pca", dims = 1:30)
DBA.E6.bin50.iDSC <- FindClusters(DBA.E6.bin50.iDSC, verbose = T, resolution = 0.3)
DBA.E6.bin50.iDSC <- RunUMAP(DBA.E6.bin50.iDSC, reduction = "pca", dims = 1:30)
anchors <- FindTransferAnchors(reference = DBA.1.iDSC, query = DBA.E6.bin50.iDSC, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.iDSC$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E6.bin50.iDSC[["pca"]], dims = 1:30)
DBA.E6.bin50.iDSC[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E6.bin50.iDSC) <- "predictions"
save(DBA.E6.bin50.iDSC, file = "20220622-Comments/05.DBA/DBA-E6-iDSC-subclusters-spatial.RData")
DBA.E6.bin50.iDSC.celltype <- GetAssayData(DBA.E6.bin50.iDSC, slot = "data")
DBA.E6.bin50.iDSC.celltype <- apply(DBA.E6.bin50.iDSC.celltype, 2, function(x) rownames(DBA.E6.bin50.iDSC.celltype)[which.max(x)])
DBA.E6.bin50.iDSC.celltype <- data.frame(cell_type = DBA.E6.bin50.iDSC.celltype, 
                                         cell_name = names(DBA.E6.bin50.iDSC.celltype), 
                                         GetTissueCoordinates(DBA.E6.bin50.iDSC), 
                                         spatial_clusters = DBA.E6.bin50.iDSC$seurat_clusters,
                                         nFeature_RNA = DBA.E6.bin50.iDSC$nFeature_Spatial, 
                                         nCount_RNA = DBA.E6.bin50.iDSC$nCount_Spatial)
ggplot() + 
  geom_point(DBA.E6.bin50.iDSC.celltype, mapping = aes(x = y, y = -x, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#F1
load("20220622-Comments/05.DBA/DBA-F1-data-map-20221212.RData")
DBA.F1.bin50.celltype <- GetAssayData(DBA.F1.bin50, slot = "data")
DBA.F1.bin50.celltype <- apply(DBA.F1.bin50.celltype, 2, function(x) rownames(DBA.F1.bin50.celltype)[which.max(x)])
DBA.F1.bin50$predictions_cell <- DBA.F1.bin50.celltype
DBA.F1.bin50.iDSC <- subset(DBA.F1.bin50, subset = predictions_cell %in% c("iDSC-DBA1"))
DefaultAssay(DBA.F1.bin50.iDSC) <- "SCT"
DBA.F1.bin50.iDSC <- SCTransform(DBA.F1.bin50.iDSC, assay = "Spatial", verbose = FALSE)
DBA.F1.bin50.iDSC <- RunPCA(DBA.F1.bin50.iDSC, assay = "SCT", verbose = F)
ElbowPlot(DBA.F1.bin50.iDSC, ndims = 50)
DBA.F1.bin50.iDSC <- FindNeighbors(DBA.F1.bin50.iDSC, reduction = "pca", dims = 1:30)
DBA.F1.bin50.iDSC <- FindClusters(DBA.F1.bin50.iDSC, verbose = T, resolution = 0.3)
DBA.F1.bin50.iDSC <- RunUMAP(DBA.F1.bin50.iDSC, reduction = "pca", dims = 1:30)
anchors <- FindTransferAnchors(reference = DBA.1.iDSC, query = DBA.F1.bin50.iDSC, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.iDSC$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.F1.bin50.iDSC[["pca"]], dims = 1:30)
DBA.F1.bin50.iDSC[["predictions"]] <- predictions.assay
DefaultAssay(DBA.F1.bin50.iDSC) <- "predictions"
save(DBA.F1.bin50.iDSC, file = "20220622-Comments/05.DBA/DBA-F1-iDSC-subclusters-spatial.RData")
DBA.F1.bin50.iDSC.celltype <- GetAssayData(DBA.F1.bin50.iDSC, slot = "data")
DBA.F1.bin50.iDSC.celltype <- apply(DBA.F1.bin50.iDSC.celltype, 2, function(x) rownames(DBA.F1.bin50.iDSC.celltype)[which.max(x)])
DBA.F1.bin50.iDSC.celltype <- data.frame(cell_type = DBA.F1.bin50.iDSC.celltype, 
                                         cell_name = names(DBA.F1.bin50.iDSC.celltype), 
                                         GetTissueCoordinates(DBA.F1.bin50.iDSC), 
                                         spatial_clusters = DBA.F1.bin50.iDSC$seurat_clusters,
                                         nFeature_RNA = DBA.F1.bin50.iDSC$nFeature_Spatial, 
                                         nCount_RNA = DBA.F1.bin50.iDSC$nCount_Spatial)
ggplot() + 
  geom_point(DBA.F1.bin50.iDSC.celltype, mapping = aes(x = x, y = -y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#embryo and maternal integration sperately
#E3
load("20220622-Comments/05.DBA/DBA-E3-data.Rdata")
E3.embryo.cellnames <- DBA.E3.bin50.celltype[DBA.E3.bin50.celltype$spatial_clusters %in% c(3), ]$cell_name
E3.maternal.cellnames <- DBA.E3.bin50.celltype[DBA.E3.bin50.celltype$spatial_clusters != 3, ]$cell_name
DBA.E3.bin50.embryo <- subset(DBA.E3.bin50, subset = seurat_clusters %in% c(3))
DBA.E3.bin50.maternal <- subset(DBA.E3.bin50, subset = seurat_clusters != 3)
#transfer - embryo
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte.embryo, query = DBA.E3.bin50.embryo, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte.embryo$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E3.bin50.embryo[["pca"]], dims = 1:30)
DBA.E3.bin50.embryo[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E3.bin50.embryo) <- "predictions"
save(DBA.E3.bin50.embryo, file = "20220622-Comments/05.DBA/DBA-E3-data-map-embryo-20221216.RData")
#transfer - maternal
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte.maternal, query = DBA.E3.bin50.maternal, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte.maternal$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E3.bin50.maternal[["pca"]], dims = 1:30)
DBA.E3.bin50.maternal[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E3.bin50.maternal) <- "predictions"
save(DBA.E3.bin50.maternal, file = "20220622-Comments/05.DBA/DBA-E3-data-map-maternal-20221216.RData")

DBA.E3.bin50.embryo.celltype <- GetAssayData(DBA.E3.bin50.embryo, slot = "data")
DBA.E3.bin50.embryo.celltype <- apply(DBA.E3.bin50.embryo.celltype, 2, function(x) rownames(DBA.E3.bin50.embryo.celltype)[which.max(x)])
DBA.E3.bin50.embryo.celltype <- data.frame(cell_type = DBA.E3.bin50.embryo.celltype, 
                                         cell_name = names(DBA.E3.bin50.embryo.celltype), 
                                         GetTissueCoordinates(DBA.E3.bin50.embryo), 
                                         spatial_clusters = DBA.E3.bin50.embryo$seurat_clusters,
                                         nFeature_RNA = DBA.E3.bin50.embryo$nFeature_Spatial, 
                                         nCount_RNA = DBA.E3.bin50.embryo$nCount_Spatial)
DBA.E3.bin50.maternal.celltype <- GetAssayData(DBA.E3.bin50.maternal, slot = "data")
DBA.E3.bin50.maternal.celltype <- apply(DBA.E3.bin50.maternal.celltype, 2, function(x) rownames(DBA.E3.bin50.maternal.celltype)[which.max(x)])
DBA.E3.bin50.maternal.celltype <- data.frame(cell_type = DBA.E3.bin50.maternal.celltype, 
                                           cell_name = names(DBA.E3.bin50.maternal.celltype), 
                                           GetTissueCoordinates(DBA.E3.bin50.maternal), 
                                           spatial_clusters = DBA.E3.bin50.maternal$seurat_clusters,
                                           nFeature_RNA = DBA.E3.bin50.maternal$nFeature_Spatial, 
                                           nCount_RNA = DBA.E3.bin50.maternal$nCount_Spatial)
DBA.E3.bin50.celltype <- rbind(DBA.E3.bin50.maternal.celltype, DBA.E3.bin50.embryo.celltype)
write.table(DBA.E3.bin50.celltype, file = "20220622-Comments/05.DBA/DBA-E3-bin50-celltype.csv", sep = ",", quote = F, row.names = F, col.names = T)
ggplot() + 
  geom_point(DBA.E3.bin50.maternal.celltype, mapping = aes(x = x, y = y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#E6
load("20220622-Comments/05.DBA/DBA-E6-data.Rdata")
DBA.E6.bin50.embryo <- subset(DBA.E6.bin50, subset = seurat_clusters %in% c(5,6,8,9,12))
DBA.E6.bin50.maternal <- subset(DBA.E6.bin50, subset = seurat_clusters %in% c(0,1,2,3,4,7,10,11))
#transfer - embryo
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte.embryo, query = DBA.E6.bin50.embryo, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte.embryo$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E6.bin50.embryo[["pca"]], dims = 1:30)
DBA.E6.bin50.embryo[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E6.bin50.embryo) <- "predictions"
save(DBA.E6.bin50.embryo, file = "20220622-Comments/05.DBA/DBA-E6-data-map-embryo-20221216.RData")
#transfer - maternal
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte.maternal, query = DBA.E6.bin50.maternal, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte.maternal$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E6.bin50.maternal[["pca"]], dims = 1:30)
DBA.E6.bin50.maternal[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E6.bin50.maternal) <- "predictions"
save(DBA.E6.bin50.maternal, file = "20220622-Comments/05.DBA/DBA-E6-data-map-maternal-20221216.RData")

DBA.E6.bin50.embryo.celltype <- GetAssayData(DBA.E6.bin50.embryo, slot = "data")
DBA.E6.bin50.embryo.celltype <- apply(DBA.E6.bin50.embryo.celltype, 2, function(x) rownames(DBA.E6.bin50.embryo.celltype)[which.max(x)])
DBA.E6.bin50.embryo.celltype <- data.frame(cell_type = DBA.E6.bin50.embryo.celltype, 
                                           cell_name = names(DBA.E6.bin50.embryo.celltype), 
                                           GetTissueCoordinates(DBA.E6.bin50.embryo), 
                                           spatial_clusters = DBA.E6.bin50.embryo$seurat_clusters,
                                           nFeature_RNA = DBA.E6.bin50.embryo$nFeature_Spatial, 
                                           nCount_RNA = DBA.E6.bin50.embryo$nCount_Spatial)
DBA.E6.bin50.maternal.celltype <- GetAssayData(DBA.E6.bin50.maternal, slot = "data")
DBA.E6.bin50.maternal.celltype <- apply(DBA.E6.bin50.maternal.celltype, 2, function(x) rownames(DBA.E6.bin50.maternal.celltype)[which.max(x)])
DBA.E6.bin50.maternal.celltype <- data.frame(cell_type = DBA.E6.bin50.maternal.celltype, 
                                             cell_name = names(DBA.E6.bin50.maternal.celltype), 
                                             GetTissueCoordinates(DBA.E6.bin50.maternal), 
                                             spatial_clusters = DBA.E6.bin50.maternal$seurat_clusters,
                                             nFeature_RNA = DBA.E6.bin50.maternal$nFeature_Spatial, 
                                             nCount_RNA = DBA.E6.bin50.maternal$nCount_Spatial)
DBA.E6.bin50.celltype <- rbind(DBA.E6.bin50.maternal.celltype, DBA.E6.bin50.embryo.celltype)
write.table(DBA.E6.bin50.celltype, file = "20220622-Comments/05.DBA/DBA-E6-bin50-celltype.csv", sep = ",", quote = F, row.names = F, col.names = T)
ggplot() + 
  geom_point(DBA.E6.bin50.celltype, mapping = aes(x = y, y = -x, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#F1
load("20220622-Comments/05.DBA/DBA-F1-data.Rdata")
DBA.F1.bin50.embryo <- subset(DBA.F1.bin50, subset = seurat_clusters %in% c(9))
DBA.F1.bin50.maternal <- subset(DBA.F1.bin50, subset = seurat_clusters %in% c(0:8,10:13))
#transfer - embryo
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte.embryo, query = DBA.F1.bin50.embryo, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte.embryo$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.F1.bin50.embryo[["pca"]], dims = 1:30)
DBA.F1.bin50.embryo[["predictions"]] <- predictions.assay
DefaultAssay(DBA.F1.bin50.embryo) <- "predictions"
save(DBA.F1.bin50.embryo, file = "20220622-Comments/05.DBA/DBA-F1-data-map-embryo-20221216.RData")
#transfer - maternal
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte.maternal, query = DBA.F1.bin50.maternal, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte.maternal$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.F1.bin50.maternal[["pca"]], dims = 1:30)
DBA.F1.bin50.maternal[["predictions"]] <- predictions.assay
DefaultAssay(DBA.F1.bin50.maternal) <- "predictions"
save(DBA.F1.bin50.maternal, file = "20220622-Comments/05.DBA/DBA-F1-data-map-maternal-20221216.RData")

DBA.F1.bin50.embryo.celltype <- GetAssayData(DBA.F1.bin50.embryo, slot = "data")
DBA.F1.bin50.embryo.celltype <- apply(DBA.F1.bin50.embryo.celltype, 2, function(x) rownames(DBA.F1.bin50.embryo.celltype)[which.max(x)])
DBA.F1.bin50.embryo.celltype <- data.frame(cell_type = DBA.F1.bin50.embryo.celltype, 
                                           cell_name = names(DBA.F1.bin50.embryo.celltype), 
                                           GetTissueCoordinates(DBA.F1.bin50.embryo), 
                                           spatial_clusters = DBA.F1.bin50.embryo$seurat_clusters,
                                           nFeature_RNA = DBA.F1.bin50.embryo$nFeature_Spatial, 
                                           nCount_RNA = DBA.F1.bin50.embryo$nCount_Spatial)
DBA.F1.bin50.maternal.celltype <- GetAssayData(DBA.F1.bin50.maternal, slot = "data")
DBA.F1.bin50.maternal.celltype <- apply(DBA.F1.bin50.maternal.celltype, 2, function(x) rownames(DBA.F1.bin50.maternal.celltype)[which.max(x)])
DBA.F1.bin50.maternal.celltype <- data.frame(cell_type = DBA.F1.bin50.maternal.celltype, 
                                             cell_name = names(DBA.F1.bin50.maternal.celltype), 
                                             GetTissueCoordinates(DBA.F1.bin50.maternal), 
                                             spatial_clusters = DBA.F1.bin50.maternal$seurat_clusters,
                                             nFeature_RNA = DBA.F1.bin50.maternal$nFeature_Spatial, 
                                             nCount_RNA = DBA.F1.bin50.maternal$nCount_Spatial)
DBA.F1.bin50.celltype <- rbind(DBA.F1.bin50.maternal.celltype, DBA.F1.bin50.embryo.celltype)
write.table(DBA.F1.bin50.celltype, file = "20220622-Comments/05.DBA/DBA-F1-bin50-celltype.csv", sep = ",", quote = F, row.names = F, col.names = T)
ggplot() + 
  geom_point(DBA.F1.bin50.celltype, mapping = aes(x = -x, y = -y, color = cell_type), size = 0.001) +
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type)

#iDSC 4 subclusters integration
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221212.RData")
DimPlot(DBA.1.iDSC, label = T)
#E3
load("20220622-Comments/05.DBA/DBA-E3-data-map-maternal-20221216.RData")
DBA.E3.bin50.maternal.celltype <- GetAssayData(DBA.E3.bin50.maternal, slot = "data")
DBA.E3.bin50.maternal.celltype <- apply(DBA.E3.bin50.maternal.celltype, 2, function(x) rownames(DBA.E3.bin50.maternal.celltype)[which.max(x)])
DBA.E3.bin50.maternal$predictions_cell <- DBA.E3.bin50.maternal.celltype
DBA.E3.bin50.iDSC <- subset(DBA.E3.bin50.maternal, subset = predictions_cell %in% c("iDSC-DBA1"))
DefaultAssay(DBA.E3.bin50.iDSC) <- "SCT"
DBA.E3.bin50.iDSC <- SCTransform(DBA.E3.bin50.iDSC, assay = "Spatial", verbose = FALSE)
DBA.E3.bin50.iDSC <- RunPCA(DBA.E3.bin50.iDSC, assay = "SCT", verbose = F)
ElbowPlot(DBA.E3.bin50.iDSC, ndims = 50)
DBA.E3.bin50.iDSC <- FindNeighbors(DBA.E3.bin50.iDSC, reduction = "pca", dims = 1:30)
DBA.E3.bin50.iDSC <- FindClusters(DBA.E3.bin50.iDSC, verbose = T, resolution = 0.3)
DBA.E3.bin50.iDSC <- RunUMAP(DBA.E3.bin50.iDSC, reduction = "pca", dims = 1:30)
anchors <- FindTransferAnchors(reference = DBA.1.iDSC, query = DBA.E3.bin50.iDSC, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.iDSC$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E3.bin50.iDSC[["pca"]], dims = 1:30)
DBA.E3.bin50.iDSC[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E3.bin50.iDSC) <- "predictions"
save(DBA.E3.bin50.iDSC, file = "20220622-Comments/05.DBA/DBA-E3-iDSC-subclusters-spatial.RData")

#E6
load("20220622-Comments/05.DBA/DBA-E6-data-map-maternal-20221216.RData")
DBA.E6.bin50.maternal.celltype <- GetAssayData(DBA.E6.bin50.maternal, slot = "data")
DBA.E6.bin50.maternal.celltype <- apply(DBA.E6.bin50.maternal.celltype, 2, function(x) rownames(DBA.E6.bin50.maternal.celltype)[which.max(x)])
DBA.E6.bin50.maternal$predictions_cell <- DBA.E6.bin50.maternal.celltype
DBA.E6.bin50.iDSC <- subset(DBA.E6.bin50.maternal, subset = predictions_cell %in% c("iDSC-DBA1"))
DefaultAssay(DBA.E6.bin50.iDSC) <- "SCT"
DBA.E6.bin50.iDSC <- SCTransform(DBA.E6.bin50.iDSC, assay = "Spatial", verbose = FALSE)
DBA.E6.bin50.iDSC <- RunPCA(DBA.E6.bin50.iDSC, assay = "SCT", verbose = F)
ElbowPlot(DBA.E6.bin50.iDSC, ndims = 50)
DBA.E6.bin50.iDSC <- FindNeighbors(DBA.E6.bin50.iDSC, reduction = "pca", dims = 1:30)
DBA.E6.bin50.iDSC <- FindClusters(DBA.E6.bin50.iDSC, verbose = T, resolution = 0.3)
DBA.E6.bin50.iDSC <- RunUMAP(DBA.E6.bin50.iDSC, reduction = "pca", dims = 1:30)
anchors <- FindTransferAnchors(reference = DBA.1.iDSC, query = DBA.E6.bin50.iDSC, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.iDSC$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E6.bin50.iDSC[["pca"]], dims = 1:30)
DBA.E6.bin50.iDSC[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E6.bin50.iDSC) <- "predictions"
save(DBA.E6.bin50.iDSC, file = "20220622-Comments/05.DBA/DBA-E6-iDSC-subclusters-spatial.RData")
DBA.E6.bin50.iDSC.celltype <- GetAssayData(DBA.E6.bin50.iDSC, slot = "data")

#F1
load("20220622-Comments/05.DBA/DBA-F1-data-map-maternal-20221216.RData")
DBA.F1.bin50.maternal.celltype <- GetAssayData(DBA.F1.bin50.maternal, slot = "data")
DBA.F1.bin50.maternal.celltype <- apply(DBA.F1.bin50.maternal.celltype, 2, function(x) rownames(DBA.F1.bin50.maternal.celltype)[which.max(x)])
DBA.F1.bin50.maternal$predictions_cell <- DBA.F1.bin50.maternal.celltype
DBA.F1.bin50.iDSC <- subset(DBA.F1.bin50.maternal, subset = predictions_cell %in% c("iDSC-DBA1"))
DefaultAssay(DBA.F1.bin50.iDSC) <- "SCT"
DBA.F1.bin50.iDSC <- SCTransform(DBA.F1.bin50.iDSC, assay = "Spatial", verbose = FALSE)
DBA.F1.bin50.iDSC <- RunPCA(DBA.F1.bin50.iDSC, assay = "SCT", verbose = F)
ElbowPlot(DBA.F1.bin50.iDSC, ndims = 50)
DBA.F1.bin50.iDSC <- FindNeighbors(DBA.F1.bin50.iDSC, reduction = "pca", dims = 1:30)
DBA.F1.bin50.iDSC <- FindClusters(DBA.F1.bin50.iDSC, verbose = T, resolution = 0.3)
DBA.F1.bin50.iDSC <- RunUMAP(DBA.F1.bin50.iDSC, reduction = "pca", dims = 1:30)
anchors <- FindTransferAnchors(reference = DBA.1.iDSC, query = DBA.F1.bin50.iDSC, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.iDSC$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.F1.bin50.iDSC[["pca"]], dims = 1:30)
DBA.F1.bin50.iDSC[["predictions"]] <- predictions.assay
DefaultAssay(DBA.F1.bin50.iDSC) <- "predictions"
save(DBA.F1.bin50.iDSC, file = "20220622-Comments/05.DBA/DBA-F1-iDSC-subclusters-spatial.RData")

#Data arrange + figures
major.color <- c("#5F9EA0", "#D2691E", "#20B2AA", "#6B8E23", "#006400", "#6495ED", "#8B0000", "#4B0082", "#CD5C5C", "#DAA520", "#9370DB", "#8B4513", "#FFD700")
decidual.color <- c("#6A5ACD","#9ACD32","#FB8072","#8DD3C7","#08519C","#FDB462","#228B22")
immune.color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#2E8B57", "#CAB2D6", "#FF7F00", "#6A3D9A")
EC.color <- c("#7B68EE", "#1B9E77", "#D95F02", "#00BFFF", "#00FA9A")
tropho.color <- c("#3CB371", "#DA70D6", "#4682B4", "#483D8B", "#008B8B", "#FF8C00")

major.levels <- c("Decidual", "Immune", "EC", "SMC", "Epithelial", "TP", "Erythroid", "Mesenchymal stem cell", "Mesenchymal-epithelial transition", "Blood P", "Viseral endoderm", "Allantois mesodermal", "Fetal EC")
decidual.levels <- c("D1-eSF", "D2-Pre-DSC", "D3-Top2a", "D4-Ptn", "D5-Gatm", "D6-S100a8", "D7-Prl8a2")
immune.levels <- c("Monocyte", "Mac", "Proliferating Mac", "DC-1", "Neutrophil", "NK", "T cell", "iDSC-DBA1")
EC.levels <- c("Angiogenic EC", "Venous EC1", "Proliterating EC", "Venous EC2", "EPC")
tropho.levels <- c("tropho-0", "tropho-1", "tropho-2", "tropho-3", "tropho-4", "tropho-5")

#E3
celltype.level <- c(decidual.levels, immune.levels[c(1:8)], EC.levels, major.levels[c(4,5,8:13)], tropho.levels[1:5])
celltype.color <- c(decidual.color, immune.color[c(1:8)], EC.color, major.color[c(4,5,8:13)], tropho.color[1:5])
DBA.E3.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin50.celltype$cell_type <- factor(DBA.E3.bin50.celltype$cell_type, levels = celltype.level)
pdf("20220622-Comments/05.DBA/00.figures/15-DBA-E3-bin50-all-celltypes.pdf", width = 10, height = 9)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y, color = cell_type), size = 1) +
  scale_color_manual(values = celltype.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/15-DBA-E3-bin50-gray-background.pdf", width = 6, height = 7)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y), size = 1, color = "#E7E7E7") +
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
plots <- list()
for (c in seq(1, 33, 1)) {
  print(celltype.level[c])
  c.df <- DBA.E3.bin50.celltype[DBA.E3.bin50.celltype$cell_type == celltype.level[c], ]
  plot <- ggplot() +
    geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size= 0.1) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = celltype.color[c], size= 0.6) + 
    theme_classic() + 
    theme(axis.text = element_blank()) + 
    ggtitle(celltype.level[c])
  plots[[c]] <- plot
}
pdf("20220622-Comments/05.DBA/00.figures/16-DBA-E3-bin50-sub-celltypes.pdf", width = 15, height = 20)
wrap_plots(plots, ncol = 6)
dev.off()

#E6
celltype.level <- c(decidual.levels, immune.levels[c(1:8)], EC.levels, major.levels[c(4,5,8:13)], tropho.levels[1:5])
celltype.color <- c(decidual.color, immune.color[c(1:8)], EC.color, major.color[c(4,5,8:13)], tropho.color[1:5])
DBA.E6.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E6.bin50.celltype$cell_type <- factor(DBA.E6.bin50.celltype$cell_type, levels = celltype.level)
pdf("20220622-Comments/05.DBA/00.figures/15-DBA-E6-bin50-all-celltypes.pdf", width = 11, height = 10)
ggplot() + 
  geom_point(DBA.E6.bin50.celltype, mapping = aes(x = y, y = -x, color = cell_type), size = 1) +
  scale_color_manual(values = celltype.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/15-DBA-E6-bin50-gray-background.pdf", width = 9, height = 9)
ggplot() + 
  geom_point(DBA.E6.bin50.celltype, mapping = aes(x = y, y = -x), size = 1, color = "#E7E7E7") +
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
plots <- list()
for (c in seq(1, 33, 1)) {
  print(celltype.level[c])
  c.df <- DBA.E6.bin50.celltype[DBA.E6.bin50.celltype$cell_type == celltype.level[c], ]
  plot <- ggplot() +
    geom_point(DBA.E6.bin50.celltype, mapping = aes(x = y, y = -x), color = "#E7E7E7", size= 0.1) + 
    geom_point(c.df, mapping = aes(x = y, y = -x), color = celltype.color[c], size= 0.6) + 
    theme_classic() + 
    theme(axis.text = element_blank()) + 
    ggtitle(celltype.level[c])
  plots[[c]] <- plot
}
pdf("20220622-Comments/05.DBA/00.figures/16-DBA-E6-bin50-sub-celltypes.pdf", width = 15, height = 20)
wrap_plots(plots, ncol = 6)
dev.off()
#F1
celltype.level <- c(decidual.levels, immune.levels[c(1:6,8)], EC.levels, major.levels[c(4,5,8:12)], tropho.levels[2:5])
celltype.color <- c(decidual.color, immune.color[c(1:6,8)], EC.color, major.color[c(4,5,8:12)], tropho.color[2:5])
DBA.F1.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.F1.bin50.celltype$cell_type <- factor(DBA.F1.bin50.celltype$cell_type, levels = celltype.level)
pdf("20220622-Comments/05.DBA/00.figures/15-DBA-F1-bin50-all-celltypes.pdf", width = 9, height = 6)
ggplot() + 
  geom_point(DBA.F1.bin50.celltype, mapping = aes(x = -x, y = -y, color = cell_type), size = 0.5) +
  scale_color_manual(values = celltype.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/15-DBA-F1-bin50-gray-background.pdf", width = 5, height = 4.5)
ggplot() + 
  geom_point(DBA.F1.bin50.celltype, mapping = aes(x = -x, y = -y), size = 0.5, color = "#E7E7E7") +
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
plots <- list()
for (c in seq(1, 30, 1)) {
  print(celltype.level[c])
  c.df <- DBA.F1.bin50.celltype[DBA.F1.bin50.celltype$cell_type == celltype.level[c], ]
  plot <- ggplot() +
    geom_point(DBA.F1.bin50.celltype, mapping = aes(x = -x, y = -y), color = "#E7E7E7", size= 0.1) + 
    geom_point(c.df, mapping = aes(x = -x, y = -y), color = celltype.color[c], size= 0.6) + 
    theme_classic() + 
    theme(axis.text = element_blank()) + 
    ggtitle(celltype.level[c])
  plots[[c]] <- plot
}
pdf("20220622-Comments/05.DBA/00.figures/16-DBA-F1-bin50-sub-celltypes.pdf", width = 15, height = 15)
wrap_plots(plots, ncol = 6)
dev.off()

#heatmap score
annotation.color <- list(cell_type = c("D1-eSF"="#6A5ACD", "D2-Pre-DSC"="#9ACD32", "D3-Top2a"="#FB8072","D4-Ptn"="#8DD3C7","D5-Gatm"="#08519C",
                                       "D6-S100a8"="#FDB462","D7-Prl8a2"="#228B22", "Monocyte"="#A6CEE3",
                                       "Mac"="#1F78B4","Proliferating Mac"="#B2DF8A","DC-1"="#33A02C", "Neutrophil"="#2E8B57", "NK"="#CAB2D6",
                                       "T cell"="#FDBF6F","iDSC-DBA1"="#6A3D9A","Angiogenic EC"="#7B68EE", "Venous EC1"="#1B9E77",
                                       "Proliterating EC"="#D95F02","Venous EC2"="#00BFFF", "EPC"="#00FA9A", "SMC"="#6B8E23","Epithelial"="#006400",
                                       "Mesenchymal stem cell"="#4B0082","Mesenchymal-epithelial transition"="#CD5C5C",
                                       "Blood P"="#DAA520", "Viseral endoderm"="#9370DB", "Allantois mesodermal"="#8B4513",
                                       "Fetal EC" = "#FFD700", "tropho-0"="#3CB371",
                                       "tropho-1"="#DA70D6", "tropho-2"="#4682B4","tropho-3"="#483D8B","tropho-4"="#008B8B"))
#E3
celltype.level <- c(decidual.levels, immune.levels[c(1:8)], EC.levels, major.levels[c(4,5,8:13)], tropho.levels[1:5])
celltype.color <- c(decidual.color, immune.color[c(1:8)], EC.color, major.color[c(4,5,8:13)], tropho.color[1:5])
DBA.E3.bin50.maternal.score <- data.frame(as.matrix(GetAssayData(DBA.E3.bin50.maternal, slot = "data")), check.names = F)
DBA.E3.bin50.maternal.score$celltype <- rownames(DBA.E3.bin50.maternal.score)
DBA.E3.bin50.maternal.score <- melt(DBA.E3.bin50.maternal.score, id.vars = "celltype")
DBA.E3.bin50.embryo.score <- data.frame(as.matrix(GetAssayData(DBA.E3.bin50.embryo, slot = "data")), check.names = F)
DBA.E3.bin50.embryo.score$celltype <- rownames(DBA.E3.bin50.embryo.score)
DBA.E3.bin50.embryo.score <- melt(DBA.E3.bin50.embryo.score, id.vars = "celltype")
DBA.E3.bin50.score <- rbind(DBA.E3.bin50.maternal.score, DBA.E3.bin50.embryo.score)
DBA.E3.bin50.score <- dcast(DBA.E3.bin50.score, celltype ~ variable)
DBA.E3.bin50.score[is.na(DBA.E3.bin50.score)] <- 0
rownames(DBA.E3.bin50.score) <- DBA.E3.bin50.score$celltype
DBA.E3.bin50.score$celltype <- NULL
DBA.E3.bin50.score <- DBA.E3.bin50.score[-18, ]
DBA.E3.bin50.celltype.max <- apply(DBA.E3.bin50.score, 2, function(x) rownames(DBA.E3.bin50.score)[which.max(x)])
DBA.E3.bin50.celltype.sort <- DBA.E3.bin50.score[, names(sort(DBA.E3.bin50.celltype.max))]
col_anno <- data.frame(cell_type = DBA.E3.bin50.celltype.max)
rownames(col_anno) <- names(DBA.E3.bin50.celltype.max)
col_anno$cell_type <- factor(col_anno$cell_type, levels = celltype.level)
col_anno$cell_name <- rownames(col_anno)
col_anno.sort <- data.frame()
for (n in celltype.level) {
  col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
}
col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type)
rownames(col_anno.df) <- rownames(col_anno.sort)
color <- brewer.pal(n =11, name = "RdBu")
color <- color[c(1:4, 7:11)]
color = colorRampPalette(rev(color))(100)
DBA.E3.bin50.celltype.sort.s <- DBA.E3.bin50.celltype.sort[celltype.level, rownames(col_anno.df)]
pdf("20220622-Comments/05.DBA/00.figures/17-DBA-E3-bin50-celltype-score.pdf", width = 20, height = 10)
pheatmap(DBA.E3.bin50.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
         scale = "column", show_colnames = F, annotation_col = col_anno.df, 
         color = color)
dev.off()

#E6
celltype.level <- c(decidual.levels, immune.levels[c(1:8)], EC.levels, major.levels[c(4,5,8:13)], tropho.levels[1:5])
celltype.color <- c(decidual.color, immune.color[c(1:8)], EC.color, major.color[c(4,5,8:13)], tropho.color[1:5])
DBA.E6.bin50.maternal.score <- data.frame(as.matrix(GetAssayData(DBA.E6.bin50.maternal, slot = "data")), check.names = F)
DBA.E6.bin50.maternal.score$celltype <- rownames(DBA.E6.bin50.maternal.score)
DBA.E6.bin50.maternal.score <- melt(DBA.E6.bin50.maternal.score, id.vars = "celltype")
DBA.E6.bin50.embryo.score <- data.frame(as.matrix(GetAssayData(DBA.E6.bin50.embryo, slot = "data")), check.names = F)
DBA.E6.bin50.embryo.score$celltype <- rownames(DBA.E6.bin50.embryo.score)
DBA.E6.bin50.embryo.score <- melt(DBA.E6.bin50.embryo.score, id.vars = "celltype")
DBA.E6.bin50.score <- rbind(DBA.E6.bin50.maternal.score, DBA.E6.bin50.embryo.score)
DBA.E6.bin50.score <- dcast(DBA.E6.bin50.score, celltype ~ variable)
DBA.E6.bin50.score[is.na(DBA.E6.bin50.score)] <- 0
rownames(DBA.E6.bin50.score) <- DBA.E6.bin50.score$celltype
DBA.E6.bin50.score$celltype <- NULL
DBA.E6.bin50.score <- DBA.E6.bin50.score[-18, ]
DBA.E6.bin50.celltype.max <- apply(DBA.E6.bin50.score, 2, function(x) rownames(DBA.E6.bin50.score)[which.max(x)])
DBA.E6.bin50.celltype.sort <- DBA.E6.bin50.score[, names(sort(DBA.E6.bin50.celltype.max))]
col_anno <- data.frame(cell_type = DBA.E6.bin50.celltype.max)
rownames(col_anno) <- names(DBA.E6.bin50.celltype.max)
col_anno$cell_type <- factor(col_anno$cell_type, levels = celltype.level)
col_anno$cell_name <- rownames(col_anno)
col_anno.sort <- data.frame()
for (n in celltype.level) {
  col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
}
col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type)
rownames(col_anno.df) <- rownames(col_anno.sort)
color <- brewer.pal(n =11, name = "RdBu")
color <- color[c(1:4, 7:11)]
color = colorRampPalette(rev(color))(100)
DBA.E6.bin50.celltype.sort.s <- DBA.E6.bin50.celltype.sort[celltype.level, rownames(col_anno.df)]
pdf("20220622-Comments/05.DBA/00.figures/17-DBA-E6-bin50-celltype-score.pdf", width = 20, height = 10)
pheatmap(DBA.E6.bin50.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
         scale = "column", show_colnames = F, annotation_col = col_anno.df, 
         color = color)
dev.off()

#F1
celltype.level <- c(decidual.levels, immune.levels[c(1:6,8)], EC.levels, major.levels[c(4,5,8:12)], tropho.levels[2:5])
celltype.color <- c(decidual.color, immune.color[c(1:6,8)], EC.color, major.color[c(4,5,8:12)], tropho.color[2:5])
DBA.F1.bin50.maternal.score <- data.frame(as.matrix(GetAssayData(DBA.F1.bin50.maternal, slot = "data")), check.names = F)
DBA.F1.bin50.maternal.score$celltype <- rownames(DBA.F1.bin50.maternal.score)
DBA.F1.bin50.maternal.score <- melt(DBA.F1.bin50.maternal.score, id.vars = "celltype")
DBA.F1.bin50.embryo.score <- data.frame(as.matrix(GetAssayData(DBA.F1.bin50.embryo, slot = "data")), check.names = F)
DBA.F1.bin50.embryo.score$celltype <- rownames(DBA.F1.bin50.embryo.score)
DBA.F1.bin50.embryo.score <- melt(DBA.F1.bin50.embryo.score, id.vars = "celltype")
DBA.F1.bin50.score <- rbind(DBA.F1.bin50.maternal.score, DBA.F1.bin50.embryo.score)
DBA.F1.bin50.score <- dcast(DBA.F1.bin50.score, celltype ~ variable)
DBA.F1.bin50.score[is.na(DBA.F1.bin50.score)] <- 0
rownames(DBA.F1.bin50.score) <- DBA.F1.bin50.score$celltype
DBA.F1.bin50.score$celltype <- NULL
DBA.F1.bin50.score <- DBA.F1.bin50.score[-18, ]
DBA.F1.bin50.celltype.max <- apply(DBA.F1.bin50.score, 2, function(x) rownames(DBA.F1.bin50.score)[which.max(x)])
DBA.F1.bin50.celltype.sort <- DBA.F1.bin50.score[, names(sort(DBA.F1.bin50.celltype.max))]
col_anno <- data.frame(cell_type = DBA.F1.bin50.celltype.max)
rownames(col_anno) <- names(DBA.F1.bin50.celltype.max)
col_anno$cell_type <- factor(col_anno$cell_type, levels = celltype.level)
col_anno$cell_name <- rownames(col_anno)
col_anno.sort <- data.frame()
for (n in celltype.level) {
  col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
}
col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type)
rownames(col_anno.df) <- rownames(col_anno.sort)
color <- brewer.pal(n =11, name = "RdBu")
color <- color[c(1:4, 7:11)]
color = colorRampPalette(rev(color))(100)
DBA.F1.bin50.celltype.sort.s <- DBA.F1.bin50.celltype.sort[celltype.level, rownames(col_anno.df)]
pdf("20220622-Comments/05.DBA/00.figures/17-DBA-F1-bin50-celltype-score.pdf", width = 20, height = 10)
pheatmap(DBA.F1.bin50.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
         scale = "column", show_colnames = F, annotation_col = col_anno.df, 
         color = color)
dev.off()

#iDSC subclusters
iDSC.color <- c("#DB7093", "#BD7FF8", "#87AD34", "#56BCC2")
#E3
load("20220622-Comments/05.DBA/DBA-E3-iDSC-subclusters-spatial.RData")
DBA.E3.bin50.iDSC.celltype <- GetAssayData(DBA.E3.bin50.iDSC, slot = "data")
DBA.E3.bin50.iDSC.celltype <- apply(DBA.E3.bin50.iDSC.celltype, 2, function(x) rownames(DBA.E3.bin50.iDSC.celltype)[which.max(x)])
DBA.E3.bin50.iDSC.celltype <- data.frame(cell_type = DBA.E3.bin50.iDSC.celltype, 
                                         cell_name = names(DBA.E3.bin50.iDSC.celltype), 
                                         GetTissueCoordinates(DBA.E3.bin50.iDSC), 
                                         spatial_clusters = DBA.E3.bin50.iDSC$seurat_clusters,
                                         nFeature_RNA = DBA.E3.bin50.iDSC$nFeature_Spatial, 
                                         nCount_RNA = DBA.E3.bin50.iDSC$nCount_Spatial)
write.table(DBA.E3.bin50.iDSC.celltype, file = "20220622-Comments/05.DBA/DBA-E3-bin50-iDSC-celltype.csv", sep = ",", quote = F, row.names = F, col.names = T)
pdf("20220622-Comments/05.DBA/00.figures/18-DBA-E3-bin50-iDSC-subclusters.pdf", width = 5, height = 5)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.E3.bin50.iDSC.celltype, mapping = aes(x = x, y = y, color = cell_type), size = 0.3) +
  scale_color_manual(values = iDSC.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#E6
load("20220622-Comments/05.DBA/DBA-E6-iDSC-subclusters-spatial.RData")
DBA.E6.bin50.iDSC.celltype <- GetAssayData(DBA.E6.bin50.iDSC, slot = "data")
DBA.E6.bin50.iDSC.celltype <- apply(DBA.E6.bin50.iDSC.celltype, 2, function(x) rownames(DBA.E6.bin50.iDSC.celltype)[which.max(x)])
DBA.E6.bin50.iDSC.celltype <- data.frame(cell_type = DBA.E6.bin50.iDSC.celltype, 
                                         cell_name = names(DBA.E6.bin50.iDSC.celltype), 
                                         GetTissueCoordinates(DBA.E6.bin50.iDSC), 
                                         spatial_clusters = DBA.E6.bin50.iDSC$seurat_clusters,
                                         nFeature_RNA = DBA.E6.bin50.iDSC$nFeature_Spatial, 
                                         nCount_RNA = DBA.E6.bin50.iDSC$nCount_Spatial)
write.table(DBA.E6.bin50.iDSC.celltype, file = "20220622-Comments/05.DBA/DBA-E6-bin50-iDSC-celltype.csv", sep = ",", quote = F, row.names = F, col.names = T)
pdf("20220622-Comments/05.DBA/00.figures/18-DBA-E6-bin50-iDSC-subclusters.pdf", width = 6.2, height = 6)
ggplot() + 
  geom_point(DBA.E6.bin50.celltype, mapping = aes(x = y, y = -x), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.E6.bin50.iDSC.celltype, mapping = aes(x = y, y = -x, color = cell_type), size = 0.3) +
  scale_color_manual(values = iDSC.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
#F1
load("20220622-Comments/05.DBA/DBA-F1-iDSC-subclusters-spatial.RData")
DBA.F1.bin50.iDSC.celltype <- GetAssayData(DBA.F1.bin50.iDSC, slot = "data")
DBA.F1.bin50.iDSC.celltype <- apply(DBA.F1.bin50.iDSC.celltype, 2, function(x) rownames(DBA.F1.bin50.iDSC.celltype)[which.max(x)])
DBA.F1.bin50.iDSC.celltype <- data.frame(cell_type = DBA.F1.bin50.iDSC.celltype, 
                                         cell_name = names(DBA.F1.bin50.iDSC.celltype), 
                                         GetTissueCoordinates(DBA.F1.bin50.iDSC), 
                                         spatial_clusters = DBA.F1.bin50.iDSC$seurat_clusters,
                                         nFeature_RNA = DBA.F1.bin50.iDSC$nFeature_Spatial, 
                                         nCount_RNA = DBA.F1.bin50.iDSC$nCount_Spatial)
write.table(DBA.F1.bin50.iDSC.celltype, file = "20220622-Comments/05.DBA/DBA-F1-bin50-iDSC-celltype.csv", sep = ",", quote = F, row.names = F, col.names = T)
pdf("20220622-Comments/05.DBA/00.figures/18-DBA-F1-bin50-iDSC-subclusters.pdf", width = 5.4, height = 5)
ggplot() + 
  geom_point(DBA.F1.bin50.celltype, mapping = aes(x = -x, y = -y), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.F1.bin50.iDSC.celltype, mapping = aes(x = -x, y = -y, color = cell_type), size = 0.3) +
  scale_color_manual(values = iDSC.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()


################################ bin20  ##############################
load("20220622-Comments/05.DBA/DBA.1-merge-embryo-SCT.RData")
#E3
DBA.E3.data <- read.table("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/05.DBA/SS200000910TR_E3.TissueCut.gem.gz", sep = "\t", header = T, check.names = F)
DBA.E3.data.bin20 <- DBA.E3.data
DBA.E3.data.bin20$x <- trunc(DBA.E3.data.bin20$x/20) * 20
DBA.E3.data.bin20$y <- trunc(DBA.E3.data.bin20$y/20) * 20
DBA.E3.data.bin20$cellID <- paste(DBA.E3.data.bin20$x, "_", DBA.E3.data.bin20$y, sep = "")
DBA.E3.data.bin20 <- aggregate(DBA.E3.data.bin20$MIDCount, by = list(DBA.E3.data.bin20$cellID, DBA.E3.data.bin20$geneID), sum)
colnames(DBA.E3.data.bin20) <- c("cellID", "geneID", "MIDCount")
DBA.E3.data.bin20$cellInx <- match(DBA.E3.data.bin20$cellID, unique(DBA.E3.data.bin20$cellID))
DBA.E3.data.bin20$cellInx <- match(DBA.E3.data.bin20$cellID, unique(DBA.E3.data.bin20$cellID))
DBA.E3.data.bin20$geneInx <- match(DBA.E3.data.bin20$geneID, unique(DBA.E3.data.bin20$geneID))
mat <- sparseMatrix(i = DBA.E3.data.bin20$geneInx, j = DBA.E3.data.bin20$cellInx, x = DBA.E3.data.bin20$MIDCount, 
                    dimnames = list(unique(DBA.E3.data.bin20$geneID), unique(DBA.E3.data.bin20$cellID)))
DBA.E3.bin20.coord.df <- data.frame(cellname = colnames(mat))
rownames(DBA.E3.bin20.coord.df) <- DBA.E3.bin20.coord.df$cellname
DBA.E3.bin20.coord.df <- separate(DBA.E3.bin20.coord.df, col = cellname, sep = "_", into = c("x", "y"))
DBA.E3.bin20 <- CreateSeuratObject(mat, project = "DBA.E3.bin20", assay = "Spatial")
DBA.E3.bin20$slice <- 1
DBA.E3.bin20$region <- "DBA.E3.bin20"
colnames(DBA.E3.bin20.coord.df) <- c("imagerow", "imagecol")
DBA.E3.bin20.coord.df$imagerow <- as.numeric(DBA.E3.bin20.coord.df$imagerow)
DBA.E3.bin20.coord.df$imagecol <- as.numeric(DBA.E3.bin20.coord.df$imagecol)
DBA.E3.bin20@images$DBA.E3_bin20 <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = DBA.E3.bin20.coord.df)
#SCT transform
DBA.E3.bin20 <- SCTransform(DBA.E3.bin20, assay = "Spatial", verbose = FALSE, ncells = 5000)
DBA.E3.bin20 <- RunPCA(DBA.E3.bin20, assay = "SCT", verbose = F)
ElbowPlot(DBA.E3.bin20, ndims = 50)
DBA.E3.bin20 <- FindNeighbors(DBA.E3.bin20, reduction = "pca", dims = 1:30)
DBA.E3.bin20 <- FindClusters(DBA.E3.bin20, verbose = T, resolution = 0.3)
DBA.E3.bin20 <- RunUMAP(DBA.E3.bin20, reduction = "pca", dims = 1:30)
DimPlot(DBA.E3.bin20, label = T)
SpatialDimPlot(DBA.E3.bin20, label = T)
save(DBA.E3.bin20, file = "20220622-Comments/05.DBA/DBA-E3-data-bin20.Rdata")
rm(DBA.E3.bin20)
gc()
#transfer E3.bin20 and single-cell RNA-seq
load("20220622-Comments/05.DBA/DBA-E3-data-bin20.Rdata")
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte, query = DBA.E3.bin20, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E3.bin20[["pca"]], dims = 1:30)
DBA.E3.bin20[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E3.bin20) <- "predictions"
save(DBA.E3.bin20, file = "20220622-Comments/05.DBA/DBA-E3-data-map-bin20-20221222.RData")

#E6
DBA.E6.data <- read.table("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/05.DBA/SS200000910TR_E6.TissueCut.gem.gz", sep = "\t", header = T, check.names = F)
DBA.E6.data.bin20 <- DBA.E6.data
DBA.E6.data.bin20$x <- trunc(DBA.E6.data.bin20$x/20) * 20
DBA.E6.data.bin20$y <- trunc(DBA.E6.data.bin20$y/20) * 20
DBA.E6.data.bin20$cellID <- paste(DBA.E6.data.bin20$x, "_", DBA.E6.data.bin20$y, sep = "")
DBA.E6.data.bin20 <- aggregate(DBA.E6.data.bin20$MIDCount, by = list(DBA.E6.data.bin20$cellID, DBA.E6.data.bin20$geneID), sum)
colnames(DBA.E6.data.bin20) <- c("cellID", "geneID", "MIDCount")
DBA.E6.data.bin20$cellInx <- match(DBA.E6.data.bin20$cellID, unique(DBA.E6.data.bin20$cellID))
DBA.E6.data.bin20$cellInx <- match(DBA.E6.data.bin20$cellID, unique(DBA.E6.data.bin20$cellID))
DBA.E6.data.bin20$geneInx <- match(DBA.E6.data.bin20$geneID, unique(DBA.E6.data.bin20$geneID))
mat <- sparseMatrix(i = DBA.E6.data.bin20$geneInx, j = DBA.E6.data.bin20$cellInx, x = DBA.E6.data.bin20$MIDCount, 
                    dimnames = list(unique(DBA.E6.data.bin20$geneID), unique(DBA.E6.data.bin20$cellID)))
DBA.E6.bin20.coord.df <- data.frame(cellname = colnames(mat))
rownames(DBA.E6.bin20.coord.df) <- DBA.E6.bin20.coord.df$cellname
DBA.E6.bin20.coord.df <- separate(DBA.E6.bin20.coord.df, col = cellname, sep = "_", into = c("x", "y"))
DBA.E6.bin20 <- CreateSeuratObject(mat, project = "DBA.E6.bin20", assay = "Spatial")
DBA.E6.bin20$slice <- 1
DBA.E6.bin20$region <- "DBA.E6.bin20"
colnames(DBA.E6.bin20.coord.df) <- c("imagerow", "imagecol")
DBA.E6.bin20.coord.df$imagerow <- as.numeric(DBA.E6.bin20.coord.df$imagerow)
DBA.E6.bin20.coord.df$imagecol <- as.numeric(DBA.E6.bin20.coord.df$imagecol)
DBA.E6.bin20@images$DBA.E6_bin20 <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = DBA.E6.bin20.coord.df)
#SCT transform
DBA.E6.bin20 <- SCTransform(DBA.E6.bin20, assay = "Spatial", verbose = FALSE, ncells = 5000)
DBA.E6.bin20 <- RunPCA(DBA.E6.bin20, assay = "SCT", verbose = F)
ElbowPlot(DBA.E6.bin20, ndims = 50)
DBA.E6.bin20 <- FindNeighbors(DBA.E6.bin20, reduction = "pca", dims = 1:30)
DBA.E6.bin20 <- FindClusters(DBA.E6.bin20, verbose = T, resolution = 0.3)
DBA.E6.bin20 <- RunUMAP(DBA.E6.bin20, reduction = "pca", dims = 1:30)
DimPlot(DBA.E6.bin20, label = T)
SpatialDimPlot(DBA.E6.bin20, label = T)
save(DBA.E6.bin20, file = "20220622-Comments/05.DBA/DBA-E6-data-bin20.Rdata")
#transfer E6.bin20 and single-cell RNA-seq
load("20220622-Comments/05.DBA/DBA-E6-data-bin20.Rdata")
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte, query = DBA.E6.bin20, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.E6.bin20[["pca"]], dims = 1:30)
DBA.E6.bin20[["predictions"]] <- predictions.assay
DefaultAssay(DBA.E6.bin20) <- "predictions"
save(DBA.E6.bin20, file = "20220622-Comments/05.DBA/DBA-E6-data-map-bin20-20221222.RData")

#F1
DBA.F1.data <- read.table("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/05.DBA/SS200000910TR_F1.TissueCut.gem.gz", sep = "\t", header = T, check.names = F)
DBA.F1.data.bin20 <- DBA.F1.data
DBA.F1.data.bin20$x <- trunc(DBA.F1.data.bin20$x/20) * 20
DBA.F1.data.bin20$y <- trunc(DBA.F1.data.bin20$y/20) * 20
DBA.F1.data.bin20$cellID <- paste(DBA.F1.data.bin20$x, "_", DBA.F1.data.bin20$y, sep = "")
DBA.F1.data.bin20 <- aggregate(DBA.F1.data.bin20$MIDCount, by = list(DBA.F1.data.bin20$cellID, DBA.F1.data.bin20$geneID), sum)
colnames(DBA.F1.data.bin20) <- c("cellID", "geneID", "MIDCount")
DBA.F1.data.bin20$cellInx <- match(DBA.F1.data.bin20$cellID, unique(DBA.F1.data.bin20$cellID))
DBA.F1.data.bin20$cellInx <- match(DBA.F1.data.bin20$cellID, unique(DBA.F1.data.bin20$cellID))
DBA.F1.data.bin20$geneInx <- match(DBA.F1.data.bin20$geneID, unique(DBA.F1.data.bin20$geneID))
mat <- sparseMatrix(i = DBA.F1.data.bin20$geneInx, j = DBA.F1.data.bin20$cellInx, x = DBA.F1.data.bin20$MIDCount, 
                    dimnames = list(unique(DBA.F1.data.bin20$geneID), unique(DBA.F1.data.bin20$cellID)))
DBA.F1.bin20.coord.df <- data.frame(cellname = colnames(mat))
rownames(DBA.F1.bin20.coord.df) <- DBA.F1.bin20.coord.df$cellname
DBA.F1.bin20.coord.df <- separate(DBA.F1.bin20.coord.df, col = cellname, sep = "_", into = c("x", "y"))
DBA.F1.bin20 <- CreateSeuratObject(mat, project = "DBA.F1.bin20", assay = "Spatial")
DBA.F1.bin20$slice <- 1
DBA.F1.bin20$region <- "DBA.F1.bin20"
colnames(DBA.F1.bin20.coord.df) <- c("imagerow", "imagecol")
DBA.F1.bin20.coord.df$imagerow <- as.numeric(DBA.F1.bin20.coord.df$imagerow)
DBA.F1.bin20.coord.df$imagecol <- as.numeric(DBA.F1.bin20.coord.df$imagecol)
DBA.F1.bin20@images$DBA.F1_bin20 <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = DBA.F1.bin20.coord.df)
#SCT transform
DBA.F1.bin20 <- SCTransform(DBA.F1.bin20, assay = "Spatial", verbose = FALSE, ncells = 5000)
DBA.F1.bin20 <- RunPCA(DBA.F1.bin20, assay = "SCT", verbose = F)
ElbowPlot(DBA.F1.bin20, ndims = 50)
DBA.F1.bin20 <- FindNeighbors(DBA.F1.bin20, reduction = "pca", dims = 1:30)
DBA.F1.bin20 <- FindClusters(DBA.F1.bin20, verbose = T, resolution = 0.3)
DBA.F1.bin20 <- RunUMAP(DBA.F1.bin20, reduction = "pca", dims = 1:30)
DimPlot(DBA.F1.bin20, label = T)
SpatialDimPlot(DBA.F1.bin20, label = T)
save(DBA.F1.bin20, file = "20220622-Comments/05.DBA/DBA-F1-data-bin20.Rdata")
#transfer F1.bin20 and single-cell RNA-seq
load("20220622-Comments/05.DBA/DBA-F1-data-bin20.Rdata")
anchors <- FindTransferAnchors(reference = DBA.1.ST.inte, query = DBA.F1.bin20, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = DBA.1.ST.inte$cellcluster, 
                                  prediction.assay = TRUE, weight.reduction = DBA.F1.bin20[["pca"]], dims = 1:30)
DBA.F1.bin20[["predictions"]] <- predictions.assay
DefaultAssay(DBA.F1.bin20) <- "predictions"
save(DBA.F1.bin20, file = "20220622-Comments/05.DBA/DBA-F1-data-map-bin20-20221222.RData")

#E3
#bin20 immune cell plots
immune.color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#2E8B57", "#CAB2D6", "#FF7F00")
immune.levels <- c("Monocyte", "Mac", "Proliferating Mac", "DC-1", "Neutrophil", "NK", "T cell")
load("20220622-Comments/05.DBA/DBA-E3-data-map-bin20-20221222.RData")
DBA.E3.bin20.celltype <- GetAssayData(DBA.E3.bin20, slot = "data")
DBA.E3.bin20.celltype <- apply(DBA.E3.bin20.celltype, 2, function(x) rownames(DBA.E3.bin20.celltype)[which.max(x)])
DBA.E3.bin20.celltype <- data.frame(cell_type = DBA.E3.bin20.celltype, 
                                    cell_name = names(DBA.E3.bin20.celltype), 
                                    GetTissueCoordinates(DBA.E3.bin20), 
                                    spatial_clusters = DBA.E3.bin20$seurat_clusters,
                                    nFeature_RNA = DBA.E3.bin20$nFeature_Spatial, 
                                    nCount_RNA = DBA.E3.bin20$nCount_Spatial)
write.table(DBA.E3.bin20.celltype, file = "20220622-Comments/05.DBA/DBA-E3-bin20-celltype.csv", sep = ",", col.names = T, row.names = F, quote = F)
DBA.E3.bin20.celltype.immune <- DBA.E3.bin20.celltype[DBA.E3.bin20.celltype$cell_type %in% immune.levels, ]
DBA.E3.bin20.celltype.immune$cell_type <- factor(DBA.E3.bin20.celltype.immune$cell_type, levels = immune.levels)
pdf("20220622-Comments/05.DBA/00.figures/35.DBA-E3-bin20-immune-all.pdf", width = 10, height = 10)
ggplot() + 
  geom_point(DBA.E3.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.E3.bin20.celltype.immune, mapping = aes(x = x, y = y, color = cell_type), size = 0.1) +
  scale_color_manual(values = immune.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
plots <- list()
for (c in seq(1, 6, 1)) {
  print(immune.levels[c])
  c.df <- DBA.E3.bin20.celltype.immune[DBA.E3.bin20.celltype.immune$cell_type == immune.levels[c], ]
  plot <- ggplot() +
    geom_point(DBA.E3.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size= 0.1) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = immune.color[c], size= 0.1) + 
    theme_classic() + 
    theme(axis.text = element_blank()) + 
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20220622-Comments/05.DBA/00.figures/35.DBA-E3-bin20-immune-subcluster.pdf", width = 10, height = 10)
wrap_plots(plots, ncol = 3)
dev.off()

#E6
load("20220622-Comments/05.DBA/DBA-E6-data-map-bin20-20221222.RData")
DBA.E6.bin20.celltype <- GetAssayData(DBA.E6.bin20, slot = "data")
DBA.E6.bin20.celltype <- apply(DBA.E6.bin20.celltype, 2, function(x) rownames(DBA.E6.bin20.celltype)[which.max(x)])
DBA.E6.bin20.celltype <- data.frame(cell_type = DBA.E6.bin20.celltype, 
                                    cell_name = names(DBA.E6.bin20.celltype), 
                                    GetTissueCoordinates(DBA.E6.bin20), 
                                    spatial_clusters = DBA.E6.bin20$seurat_clusters,
                                    nFeature_RNA = DBA.E6.bin20$nFeature_Spatial, 
                                    nCount_RNA = DBA.E6.bin20$nCount_Spatial)
write.table(DBA.E6.bin20.celltype, file = "20220622-Comments/05.DBA/DBA-E6-bin20-celltype.csv", sep = ",", col.names = T, row.names = F, quote = F)
DBA.E6.bin20.celltype.immune <- DBA.E6.bin20.celltype[DBA.E6.bin20.celltype$cell_type %in% immune.levels, ]
DBA.E6.bin20.celltype.immune$cell_type <- factor(DBA.E6.bin20.celltype.immune$cell_type, levels = immune.levels)
pdf("20220622-Comments/05.DBA/00.figures/35.DBA-E6-bin20-immune-all.pdf", width = 10, height = 10)
ggplot() + 
  geom_point(DBA.E6.bin20.celltype, mapping = aes(x = -y, y = -x), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.E6.bin20.celltype.immune, mapping = aes(x = -y, y = -x, color = cell_type), size = 0.1) +
  scale_color_manual(values = immune.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
plots <- list()
for (c in seq(1, 6, 1)) {
  print(immune.levels[c])
  c.df <- DBA.E6.bin20.celltype.immune[DBA.E6.bin20.celltype.immune$cell_type == immune.levels[c], ]
  plot <- ggplot() +
    geom_point(DBA.E6.bin20.celltype, mapping = aes(x = -y, y = -x), color = "#E7E7E7", size= 0.1) + 
    geom_point(c.df, mapping = aes(x = -y, y = -x), color = immune.color[c], size= 0.1) + 
    theme_classic() + 
    theme(axis.text = element_blank()) + 
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20220622-Comments/05.DBA/00.figures/35.DBA-E6-bin20-immune-subcluster.pdf", width = 13, height = 13)
wrap_plots(plots, ncol = 3)
dev.off()

#Immune cell spatial distribution --- bin50
#E3
DBA.E3.bin50.celltype.immune <- DBA.E3.bin50.celltype.iDSCsub[DBA.E3.bin50.celltype.iDSCsub$cell_type %in% immune.levels, ]
DBA.E3.bin50.celltype.immune$cell_type <- factor(DBA.E3.bin50.celltype.immune$cell_type, levels = immune.levels)
pdf("20220622-Comments/05.DBA/00.figures/35.DBA-E3-bin50-immune-all.pdf", width = 6, height = 7)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype.iDSCsub, mapping = aes(x = x, y = y), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.E3.bin50.celltype.immune, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = immune.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#E6
DBA.E6.bin50.celltype.immune <- DBA.E6.bin50.celltype.iDSCsub[DBA.E6.bin50.celltype.iDSCsub$cell_type %in% immune.levels, ]
DBA.E6.bin50.celltype.immune$cell_type <- factor(DBA.E6.bin50.celltype.immune$cell_type, levels = immune.levels)
pdf("20220622-Comments/05.DBA/00.figures/35.DBA-E6-bin50-immune-all.pdf", width = 8, height = 7)
ggplot() + 
  geom_point(DBA.E6.bin50.celltype.iDSCsub, mapping = aes(x = -y, y = -x), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.E6.bin50.celltype.immune, mapping = aes(x = -y, y = -x, color = cell_type), size = 0.5) +
  scale_color_manual(values = immune.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#F1
DBA.F1.bin50.celltype.immune <- DBA.F1.bin50.celltype.iDSCsub[DBA.F1.bin50.celltype.iDSCsub$cell_type %in% immune.levels, ]
DBA.F1.bin50.celltype.immune$cell_type <- factor(DBA.F1.bin50.celltype.immune$cell_type, levels = immune.levels)
pdf("20220622-Comments/05.DBA/00.figures/35.DBA-F1-bin50-immune-all.pdf", width = 7, height = 5)
ggplot() + 
  geom_point(DBA.F1.bin50.celltype.iDSCsub, mapping = aes(x = -y, y = -x), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.F1.bin50.celltype.immune, mapping = aes(x = -y, y = -x, color = cell_type), size = 0.5) +
  scale_color_manual(values = immune.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()


load("20220622-Comments/05.DBA/DBA-F1-data-map-bin20-20221222.RData")
DBA.F1.bin20.celltype <- GetAssayData(DBA.F1.bin20, slot = "data")
DBA.F1.bin20.celltype <- apply(DBA.F1.bin20.celltype, 2, function(x) rownames(DBA.F1.bin20.celltype)[which.max(x)])
DBA.F1.bin20.celltype <- data.frame(cell_type = DBA.F1.bin20.celltype, 
                                    cell_name = names(DBA.F1.bin20.celltype), 
                                    GetTissueCoordinates(DBA.F1.bin20), 
                                    spatial_clusters = DBA.F1.bin20$seurat_clusters,
                                    nFeature_RNA = DBA.F1.bin20$nFeature_Spatial, 
                                    nCount_RNA = DBA.F1.bin20$nCount_Spatial)
write.table(DBA.F1.bin20.celltype, file = "20220622-Comments/05.DBA/DBA-F1-bin20-celltype.csv", sep = ",", col.names = T, row.names = F, quote = F)
DBA.F1.bin20.celltype.immune <- DBA.F1.bin20.celltype[DBA.F1.bin20.celltype$cell_type %in% immune.levels, ]
DBA.F1.bin20.celltype.immune$cell_type <- factor(DBA.F1.bin20.celltype.immune$cell_type, levels = immune.levels)
pdf("20220622-Comments/05.DBA/00.figures/35.DBA-F1-bin20-immune-all.pdf", width = 10, height = 10)
ggplot() + 
  geom_point(DBA.F1.bin20.celltype, mapping = aes(x = -x, y = -y), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.F1.bin20.celltype.immune, mapping = aes(x = -x, y = -y, color = cell_type), size = 0.1) +
  scale_color_manual(values = immune.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
plots <- list()
for (c in seq(1, 6, 1)) {
  print(immune.levels[c])
  c.df <- DBA.F1.bin20.celltype.immune[DBA.F1.bin20.celltype.immune$cell_type == immune.levels[c], ]
  plot <- ggplot() +
    geom_point(DBA.F1.bin20.celltype, mapping = aes(x = -x, y = -y), color = "#E7E7E7", size= 0.1) + 
    geom_point(c.df, mapping = aes(x = -x, y = -y), color = immune.color[c], size= 0.1) + 
    theme_classic() + 
    theme(axis.text = element_blank()) + 
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20220622-Comments/05.DBA/00.figures/35.DBA-F1-bin20-immune-subcluster.pdf", width = 13, height = 11)
wrap_plots(plots, ncol = 3)
dev.off()

#计算空间上相邻的细胞
CalcSlideClosedCells <- function(df, bin = 50){
  df$x <- abs(as.numeric(df$x))
  df$y <- abs(as.numeric(df$y))
  celltype <- as.character(unique(df$cell_type))
  result.df <- data.frame(cell_type = celltype)
  for (c in seq(1, length(celltype), 1)) {
    print(c)
    c.df <- df[df$cell_type == celltype[c], ]
    names1 <- paste(c.df$x + bin, c.df$y, sep = "_")
    names2 <- paste(c.df$x + bin, c.df$y + bin, sep = "_")
    names3 <- paste(c.df$x + bin, c.df$y - bin, sep = "_")
    names4 <- paste(c.df$x - bin, c.df$y, sep = "_")
    names5 <- paste(c.df$x - bin, c.df$y + bin, sep = "_")
    names6 <- paste(c.df$x - bin, c.df$y - bin, sep = "_")
    names7 <- paste(c.df$x, c.df$y + bin, sep = "_")
    names8 <- paste(c.df$x, c.df$y - bin, sep = "_")
    names <- unique(c(names1, names2, names3, names4, names5, names6, names7, names8))
    closed.df <- df[df$cell_name %in% names, ]
    closed.df.tb <- data.frame(table(closed.df$cell_type))
    colnames(closed.df.tb)[1] <- "cell_type"
    result.df <- merge(result.df, closed.df.tb, by = "cell_type", all = T)
  }
  rownames(result.df) <- result.df[, 1]
  result.df <- result.df[, -1]
  colnames(result.df) <- c(celltype)
  result.df[is.na(result.df)] <- 0
  return(result.df)
}
DBA.E3.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin50.celltype <- DBA.E3.bin50.celltype[, -c(5)]
DBA.E3.bin50.iDSC.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-iDSC-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin50.iDSC.celltype <- DBA.E3.bin50.iDSC.celltype[, -c(5,6)]
DBA.E3.bin50.celltype <- DBA.E3.bin50.celltype[DBA.E3.bin50.celltype$cell_type != "iDSC-DBA1", ]
DBA.E3.bin50.celltype.iDSCsub <- rbind(DBA.E3.bin50.celltype, DBA.E3.bin50.iDSC.celltype)
E3.closed <- CalcSlideClosedCells(DBA.E3.bin50.celltype.iDSCsub)
E3.closed.r <- apply(E3.closed, 2, function(x){x/rowSums(E3.closed)})
E3.closed.r <- melt(E3.closed.r)
colnames(E3.closed.r) <- c("source", "bin50.closed", "value")
E3.closed.r$sample <- "E3"

DBA.E6.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E6.bin50.celltype <- DBA.E6.bin50.celltype[, -c(5)]
DBA.E6.bin50.iDSC.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-iDSC-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E6.bin50.iDSC.celltype <- DBA.E6.bin50.iDSC.celltype[, -c(5,6)]
DBA.E6.bin50.celltype <- DBA.E6.bin50.celltype[DBA.E6.bin50.celltype$cell_type != "iDSC-DBA1", ]
DBA.E6.bin50.celltype.iDSCsub <- rbind(DBA.E6.bin50.celltype, DBA.E6.bin50.iDSC.celltype)
E6.closed <- CalcSlideClosedCells(DBA.E6.bin50.celltype.iDSCsub)
E6.closed.r <- apply(E6.closed, 2, function(x){x/rowSums(E6.closed)})
E6.closed.r <- melt(E6.closed.r)
colnames(E6.closed.r) <- c("source", "bin50.closed", "value")
E6.closed.r$sample <- "E6"

DBA.F1.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.F1.bin50.celltype <- DBA.F1.bin50.celltype[, -c(5)]
DBA.F1.bin50.iDSC.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin50-iDSC-celltype.csv", sep = ",", check.names = F, header = T)
DBA.F1.bin50.iDSC.celltype <- DBA.F1.bin50.iDSC.celltype[, -c(5,6)]
DBA.F1.bin50.celltype <- DBA.F1.bin50.celltype[DBA.F1.bin50.celltype$cell_type != "iDSC-DBA1", ]
DBA.F1.bin50.celltype.iDSCsub <- rbind(DBA.F1.bin50.celltype, DBA.F1.bin50.iDSC.celltype)
F1.closed <- CalcSlideClosedCells(DBA.F1.bin50.celltype.iDSCsub)
F1.closed.r <- apply(F1.closed, 2, function(x){x/rowSums(F1.closed)})
F1.closed.r <- melt(F1.closed.r)
colnames(F1.closed.r) <- c("source", "bin50.closed", "value")
F1.closed.r$sample <- "F1"
DBA.bin50.closed <- Reduce(rbind, list(E3.closed.r, E6.closed.r, F1.closed.r))
write.table(DBA.bin50.closed, file = "20220622-Comments/05.DBA/DBA-slides-bin50-closed-cells-calc.csv", sep = ",", row.names = F, quote = F)

DBA.bin50.closed.iDSC <- DBA.bin50.closed[DBA.bin50.closed$source %in% c("DBA.1.iDSC0-1","DBA.1.iDSC0-2", "DBA.1.iDSC1","DBA.1.iDSC2") &
                                            DBA.bin50.closed$bin50.closed %in% c("D6-S100a8", immune.levels[1:6]), ]
DBA.bin50.closed.iDSC$source <- factor(DBA.bin50.closed.iDSC$source, levels = c("DBA.1.iDSC0-1","DBA.1.iDSC0-2", "DBA.1.iDSC1","DBA.1.iDSC2"))
DBA.bin50.closed.iDSC$bin50.closed <- factor(DBA.bin50.closed.iDSC$bin50.closed, levels = c("D6-S100a8", immune.levels[1:6]))
pdf("20220622-Comments/05.DBA/00.figures/36.DBA-iDSC-4clusters-ajacent-cell-D6-immune-barplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(DBA.bin50.closed.iDSC, mapping = aes(x = bin50.closed, y = value, fill = bin50.closed), stat = "identity", width = 1) +
  scale_fill_manual(values = c(decidual.color[6], immune.color[1:6])) +
  facet_grid(sample ~ source, space = "free") +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", color = "black"))
dev.off()

#merge all slides celltype information
DBA.E3.bin50.celltype.iDSCsub$sample <- "E3"
DBA.E6.bin50.celltype.iDSCsub$sample <- "E6"
DBA.F1.bin50.celltype.iDSCsub$sample <- "F1"
DBA.slides.def <- Reduce(rbind, list(DBA.E3.bin50.celltype.iDSCsub, DBA.E6.bin50.celltype.iDSCsub, DBA.F1.bin50.celltype.iDSCsub))
write.table(DBA.slides.def, file = "20220622-Comments/05.DBA/DBA-all-slides-bin50-celltype-df.csv", sep = ",", col.names = T, row.names = F, quote = F)
#MET
MET.cell <- c("D1-eSF", "D2-Pre-DSC")
DBA.slides.def.MET <- DBA.slides.def[DBA.slides.def$cell_type %in% MET.cell, ]
DBA.slides.def.MET.df <- table(DBA.slides.def.MET$cell_type, DBA.slides.def.MET$sample)
DBA.slides.def.MET.df <- DBA.slides.def.MET.df[MET.cell, ]
DBA.slides.def.MET.df <- apply(DBA.slides.def.MET.df, 1, function(x){x/table(DBA.slides.def$sample) * 100})
DBA.slides.def.MET.df <- melt(DBA.slides.def.MET.df)
DBA.slides.def.MET.df$zone <- "MET"
#cycling
cycling.cell <- c("D3-Top2a")
DBA.slides.def.cycling <- DBA.slides.def[DBA.slides.def$cell_type %in% cycling.cell, ]
DBA.slides.def.cycling.df <- table(DBA.slides.def.cycling$cell_type, DBA.slides.def.cycling$sample)
DBA.slides.def.cycling.df <- DBA.slides.def.cycling.df[cycling.cell, ]
DBA.slides.def.cycling.df <- DBA.slides.def.cycling.df/table(DBA.slides.def$sample) * 100
DBA.slides.def.cycling.df <- data.frame(Var1 = names(DBA.slides.def.cycling.df), Var2 = "D3-Top2a", 
                                               value = as.numeric(DBA.slides.def.cycling.df), zone = "cycling DSC")
#nourish
nourish.cell <- c("D5-Gatm", "D7-Prl8a2")
DBA.slides.def.nourish <- DBA.slides.def[DBA.slides.def$cell_type %in% nourish.cell, ]
DBA.slides.def.nourish.df <- table(DBA.slides.def.nourish$cell_type, DBA.slides.def.nourish$sample)
DBA.slides.def.nourish.df <- DBA.slides.def.nourish.df[nourish.cell, ]
DBA.slides.def.nourish.df <- apply(DBA.slides.def.nourish.df, 1, function(x){x/table(DBA.slides.def$sample) * 100})
DBA.slides.def.nourish.df <- melt(DBA.slides.def.nourish.df)
DBA.slides.def.nourish.df$zone <- "nourish"
#Immune
immune.cell <- c("Monocyte", "Mac", "Proliferating Mac", "DC-1", "Neutrophil", "NK")
DBA.slides.def.immune <- DBA.slides.def[DBA.slides.def$cell_type %in% immune.cell, ]
DBA.slides.def.immune.df <- table(DBA.slides.def.immune$cell_type, DBA.slides.def.immune$sample)
DBA.slides.def.immune.df <- DBA.slides.def.immune.df[immune.cell, ]
DBA.slides.def.immune.df <- apply(DBA.slides.def.immune.df, 1, function(x){x/table(DBA.slides.def$sample) * 100})
DBA.slides.def.immune.df <- melt(DBA.slides.def.immune.df)
DBA.slides.def.immune.df$zone <- "immune"
#EC
EC.cell <- c("Proliterating EC", "Angiogenic EC", "Venous EC1", "Venous EC2", "EPC", "D6-S100a8")
DBA.slides.def.EC <- DBA.slides.def[DBA.slides.def$cell_type %in% EC.cell, ]
DBA.slides.def.EC.df <- table(DBA.slides.def.EC$cell_type, DBA.slides.def.EC$sample)
DBA.slides.def.EC.df <- DBA.slides.def.EC.df[EC.cell, ]
DBA.slides.def.EC.df <- apply(DBA.slides.def.EC.df, 1, function(x){x/table(DBA.slides.def$sample) * 100})
DBA.slides.def.EC.df <- melt(DBA.slides.def.EC.df)
DBA.slides.def.EC.df$zone <- "EC"
#iDSC
iDSC.cell <- c("DBA.1.iDSC0-1", "DBA.1.iDSC0-2", "DBA.1.iDSC1", "DBA.1.iDSC2")
DBA.slides.def.iDSC <- DBA.slides.def[DBA.slides.def$cell_type %in% iDSC.cell, ]
DBA.slides.def.iDSC.df <- table(DBA.slides.def.iDSC$cell_type, DBA.slides.def.iDSC$sample)
DBA.slides.def.iDSC.df <- DBA.slides.def.iDSC.df[iDSC.cell, ]
DBA.slides.def.iDSC.df <- apply(DBA.slides.def.iDSC.df, 1, function(x){x/table(DBA.slides.def$sample) * 100})
DBA.slides.def.iDSC.df <- melt(DBA.slides.def.iDSC.df)
DBA.slides.def.iDSC.df$zone <- "iDSC"
#merge all
DBA.slides.def.zone.df <- Reduce(rbind, list(DBA.slides.def.MET.df, DBA.slides.def.cycling.df,
                                             DBA.slides.def.nourish.df, DBA.slides.def.immune.df, 
                                             DBA.slides.def.EC.df, DBA.slides.def.iDSC.df))
DBA.slides.def.zone.df$Var1 <- factor(DBA.slides.def.zone.df$Var1, levels = c("F1", "E3", "E6"))
DBA.slides.def.zone.df$Var2 <- factor(DBA.slides.def.zone.df$Var2, levels = c(decidual.levels, immune.levels, EC.cell[1:5], iDSC.cell))
pdf("20220622-Comments/05.DBA/00.figures/37.DBA-all-slides-zone-fraction-bin50.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(DBA.slides.def.zone.df, mapping = aes(x = Var1, y = value, fill = Var2), stat = "identity", position = "stack") + 
  facet_wrap(~ zone, scales = "free") + 
  scale_fill_manual(values = c(decidual.color[c(1:3,5:7)], immune.color[c(1:6)], EC.color[c(1:2,4:6)], rep("#6A3D9A", 4))) + 
  theme_classic()
dev.off()

#Immune zone bin20
DBA.E3.bin20.celltype$sample <- "E3"
DBA.E6.bin20.celltype$sample <- "E6"
DBA.F1.bin20.celltype$sample <- "F1"
DBA.slides.def.bin20 <- Reduce(rbind, list(DBA.E3.bin20.celltype, DBA.E6.bin20.celltype, DBA.F1.bin20.celltype))
#Immune
immune.cell <- c("Monocyte", "Mac", "Proliferating Mac", "DC-1", "Neutrophil", "NK")
DBA.slides.def.bin20.immune <- DBA.slides.def.bin20[DBA.slides.def.bin20$cell_type %in% immune.cell, ]
DBA.slides.def.bin20.immune.df <- table(DBA.slides.def.bin20.immune$cell_type, DBA.slides.def.bin20.immune$sample)
DBA.slides.def.bin20.immune.df <- DBA.slides.def.bin20.immune.df[immune.cell, ]
DBA.slides.def.bin20.immune.df <- apply(DBA.slides.def.bin20.immune.df, 1, function(x){x/table(DBA.slides.def.bin20$sample) * 100})
DBA.slides.def.bin20.immune.df <- melt(DBA.slides.def.bin20.immune.df)
DBA.slides.def.bin20.immune.df$zone <- "immune"
DBA.slides.def.bin20.immune.df$Var1 <- factor(DBA.slides.def.bin20.immune.df$Var1, levels = c("F1", "E3", "E6"))
DBA.slides.def.bin20.immune.df$Var2 <- factor(DBA.slides.def.bin20.immune.df$Var2, levels = immune.levels)
pdf("20220622-Comments/05.DBA/00.figures/37.DBA-all-slides-immune-zone-fraction-bin20.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(DBA.slides.def.bin20.immune.df, mapping = aes(x = Var1, y = value, fill = Var2), stat = "identity", position = "stack") + 
  facet_wrap(~ zone, scales = "free") + 
  scale_fill_manual(values = c(immune.color[c(1:6)])) + 
  theme_classic()
dev.off()

#EC/iDSC0-1/D6 紧邻的比例，将DBA和Normal 进行compare
#处理normal的数据，E8.5-B6 和 E9.5
E85.B6.bin50.celltype.df.sub <- read.table("20210827-figures/figure5/E85-B6-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F, header = T)
E85.B6.closed <- CalcSlideClosedCells(E85.B6.bin50.celltype.df.sub, bin = 50)
E85.B6.closed.r <- apply(E85.B6.closed, 2, function(x){x/rowSums(E85.B6.closed)})
E85.B6.closed.r <- melt(E85.B6.closed.r)
colnames(E85.B6.closed.r) <- c("source", "bin50.closed", "value")
E85.B6.closed.r$sample <- "E85.B6"

E95.bin50.celltype.df.sub <- read.table("20210827-figures/figure5/E95-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F, header = T)
E95.closed <- CalcSlideClosedCells(E95.bin50.celltype.df.sub, bin = 50)
E95.closed.r <- apply(E95.closed, 2, function(x){x/rowSums(E95.closed)})
E95.closed.r <- melt(E95.closed.r)
colnames(E95.closed.r) <- c("source", "bin50.closed", "value")
E95.closed.r$sample <- "E95"
Normal.E85.E95.closed <- rbind(E85.B6.closed.r, E95.closed.r)
write.table(Normal.E85.E95.closed, file = "20220622-Comments/05.DBA/Normal-E85-E95-slides-bin50-closed-cells-calc.csv", sep = ",", quote = F, row.names = F, col.names = T)

#rename some cell types to the same (E85.B6, E95, DBA.E3, DBA.E6, DBA.F1)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("Lum-D", "D1-eSF", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("Sfrp4-D", "D2-Pre-DSC", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("Top2a-D", "D3-Top2a", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("Ptn-D", "D4-Ptn", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("Gatm-D", "D5-Gatm", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("S100a8-D", "D6-S100a8", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("Prl8a2-D", "D7-Prl8a2", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("Mono", "Monocyte", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("Pro-Mac", "Proliferating Mac", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("^T$", "T cell", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("EC-0", "Angiogenic EC", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("EC-1", "Venous EC1", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("EC-2", "Proliferating EC", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("EC-3", "Venous EC2", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("EC-4", "Arterial EC", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("FB_immune_0", "iDSC0", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("FB_immune_1", "iDSC1", E85.B6.bin50.celltype.df.sub$cell_type)
E85.B6.bin50.celltype.df.sub$cell_type <- gsub("FB_immune_2", "iDSC2", E85.B6.bin50.celltype.df.sub$cell_type)

E95.bin50.celltype.df.sub$cell_type <- gsub("Lum-D", "D1-eSF", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("Sfrp4-D", "D2-Pre-DSC", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("Top2a-D", "D3-Top2a", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("Ptn-D", "D4-Ptn", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("Gatm-D", "D5-Gatm", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("S100a8-D", "D6-S100a8", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("Prl8a2-D", "D7-Prl8a2", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("Mono", "Monocyte", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("Pro-Mac", "Proliferating Mac", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("^T$", "T cell", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("EC-0", "Angiogenic EC", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("EC-1", "Venous EC1", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("EC-2", "Proliferating EC", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("EC-3", "Venous EC2", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("EC-4", "Arterial EC", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("FB_immune_0", "iDSC0", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("FB_immune_1", "iDSC1", E95.bin50.celltype.df.sub$cell_type)
E95.bin50.celltype.df.sub$cell_type <- gsub("FB_immune_2", "iDSC2", E95.bin50.celltype.df.sub$cell_type)

DBA.E3.bin50.celltype.iDSCsub$cell_type <- gsub("Proliterating EC", "Proliferating EC", DBA.E3.bin50.celltype.iDSCsub$cell_type)
DBA.E3.bin50.celltype.iDSCsub$cell_type <- gsub("DBA.1.iDSC0-1", "iDSC0", DBA.E3.bin50.celltype.iDSCsub$cell_type)
DBA.E3.bin50.celltype.iDSCsub$cell_type <- gsub("DBA.1.iDSC1", "iDSC1", DBA.E3.bin50.celltype.iDSCsub$cell_type)
DBA.E3.bin50.celltype.iDSCsub$cell_type <- gsub("DBA.1.iDSC2", "iDSC2", DBA.E3.bin50.celltype.iDSCsub$cell_type)

DBA.E6.bin50.celltype.iDSCsub$cell_type <- gsub("Proliterating EC", "Proliferating EC", DBA.E6.bin50.celltype.iDSCsub$cell_type)
DBA.E6.bin50.celltype.iDSCsub$cell_type <- gsub("DBA.1.iDSC0-1", "iDSC0", DBA.E6.bin50.celltype.iDSCsub$cell_type)
DBA.E6.bin50.celltype.iDSCsub$cell_type <- gsub("DBA.1.iDSC1", "iDSC1", DBA.E6.bin50.celltype.iDSCsub$cell_type)
DBA.E6.bin50.celltype.iDSCsub$cell_type <- gsub("DBA.1.iDSC2", "iDSC2", DBA.E6.bin50.celltype.iDSCsub$cell_type)

DBA.F1.bin50.celltype.iDSCsub$cell_type <- gsub("Proliterating EC", "Proliferating EC", DBA.F1.bin50.celltype.iDSCsub$cell_type)
DBA.F1.bin50.celltype.iDSCsub$cell_type <- gsub("DBA.1.iDSC0-1", "iDSC0", DBA.F1.bin50.celltype.iDSCsub$cell_type)
DBA.F1.bin50.celltype.iDSCsub$cell_type <- gsub("DBA.1.iDSC1", "iDSC1", DBA.F1.bin50.celltype.iDSCsub$cell_type)
DBA.F1.bin50.celltype.iDSCsub$cell_type <- gsub("DBA.1.iDSC2", "iDSC2", DBA.F1.bin50.celltype.iDSCsub$cell_type)

#统计EC各个细胞周围细胞的比例
CalcSlideClosedCells_single <- function(df, celltype = NULL, bin = 50){
  df$x <- abs(as.numeric(df$x))
  df$y <- abs(as.numeric(df$y))
  results.df <- data.frame()
  celltype <- intersect(unique(df$cell_type), celltype)
  for (ce in celltype) {
    results.tmp.df <- data.frame()
    cellname <- df[df$cell_type == ce, ]$cell_name
    for (c in cellname) {
      print(c)
      x = df[df$cell_name == c, ]$x
      y = df[df$cell_name == c, ]$y
      names1 <- paste(x + bin, y, sep = "_")
      names2 <- paste(x + bin, y + bin, sep = "_")
      names3 <- paste(x + bin, y - bin, sep = "_")
      names4 <- paste(x - bin, y, sep = "_")
      names5 <- paste(x - bin, y + bin, sep = "_")
      names6 <- paste(x - bin, y - bin, sep = "_")
      names7 <- paste(x, y + bin, sep = "_")
      names8 <- paste(x, y - bin, sep = "_")
      closed.df <- df[df$cell_name %in% c(names1, names2, names3, names4, names5, names6, names7, names8), ]
      closed.celltype <- table(closed.df$cell_type)/8
      results.tmp.df <- rbind(results.tmp.df, data.frame(bin50.closed = names(closed.celltype), value = as.numeric(closed.celltype)))
      if("D6-S100a8" %in% names(closed.celltype) == FALSE){
        results.tmp.df <- rbind(results.tmp.df, c("D6-S100a8", 0))
      }
      if("iDSC0" %in% names(closed.celltype) == FALSE){
        results.tmp.df <- rbind(results.tmp.df, c("iDSC0", 0))
      }
    }
    results.tmp.df$source <- ce
    results.df <- rbind(results.tmp.df, results.df)
  }
  results.df$value <- as.numeric(results.df$value)
  return(results.df)
}
EC.cell <- c("Proliferating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2", "EPC")
Normal.E85.B6 <- CalcSlideClosedCells_single(E85.B6.bin50.celltype.df.sub, celltype = EC.cell, bin = 50)
Normal.E95 <- CalcSlideClosedCells_single(E95.bin50.celltype.df.sub, celltype = EC.cell, bin = 50)
DBA.E3 <- CalcSlideClosedCells_single(DBA.E3.bin50.celltype.iDSCsub, celltype = EC.cell, bin = 50)
DBA.E6 <- CalcSlideClosedCells_single(DBA.E6.bin50.celltype.iDSCsub, celltype = EC.cell, bin = 50)
DBA.F1 <- CalcSlideClosedCells_single(DBA.F1.bin50.celltype.iDSCsub, celltype = EC.cell, bin = 50)
Normal.E85.B6$sample <- "E85.B6"
Normal.E95$sample <- "E95"
DBA.E3$sample <- "E3"
DBA.E6$sample <- "E6"
DBA.F1$sample <- "F1"

EC.single.closed <- rbind(Normal.E85.B6, Normal.E95, DBA.E3, DBA.E6, DBA.F1)
EC.single.closed <- EC.single.closed[EC.single.closed$bin50.closed %in% c("D6-S100a8", "iDSC0"), ]
ggplot() +
  geom_violin(EC.single.closed, mapping = aes(x = bin50.closed, y = value, fill = bin50.closed), scale = "width") +
  facet_grid(sample ~ source)

E3.EC <- CalcSlideClosedCells_single(DBA.E3.bin50.celltype.iDSCsub, c("Proliterating EC", "Angiogenic EC", "Venous EC1", "Venous EC2", "EPC"), bin = 50)
#E3.EC$value <- as.numeric(E3.EC$value)
E3.EC$sample <- "E3"
E3.EC <- E3.EC[E3.EC$bin50.closed %in% c("D6-S100a8", "DBA.1.iDSC0-1"), ]
ggplot() +
  geom_boxplot(E3.EC, mapping = aes(x = bin50.closed, y = value, fill = bin50.closed)) +
  facet_wrap(~ source)

#比较Normal和DBA的iDSC周围的免疫细胞的比例变化
DBA.bin50.closed <- read.table("20220622-Comments/05.DBA/DBA-slides-bin50-closed-cells-calc.csv", sep = ",", header = T, check.names = F)
Normal.E85.E95.closed <- read.table("20220622-Comments/05.DBA/Normal-E85-E95-slides-bin50-closed-cells-calc.csv", sep = ",", header = T, check.names = F)
DBA.bin50.closed$source <- gsub("Proliterating EC", "Proliferating EC", DBA.bin50.closed$source)
DBA.bin50.closed$source <- gsub("DBA.1.iDSC0-1", "iDSC0", DBA.bin50.closed$source)
DBA.bin50.closed$source <- gsub("DBA.1.iDSC1", "iDSC1", DBA.bin50.closed$source)
DBA.bin50.closed$source <- gsub("DBA.1.iDSC2", "iDSC2", DBA.bin50.closed$source)
DBA.bin50.closed$bin50.closed <- gsub("Proliterating EC", "Proliferating EC", DBA.bin50.closed$bin50.closed)
DBA.bin50.closed$bin50.closed <- gsub("DBA.1.iDSC0-1", "iDSC0", DBA.bin50.closed$bin50.closed)
DBA.bin50.closed$bin50.closed <- gsub("DBA.1.iDSC1", "iDSC1", DBA.bin50.closed$bin50.closed)
DBA.bin50.closed$bin50.closed <- gsub("DBA.1.iDSC2", "iDSC2", DBA.bin50.closed$bin50.closed)

Normal.DBA.bin50.closed <- rbind(Normal.E85.E95.closed, DBA.bin50.closed)
Normal.DBA.bin50.closed.iDSC <- Normal.DBA.bin50.closed[Normal.DBA.bin50.closed$source %in% c("iDSC0","iDSC1", "iDSC2","DBA.1.iDSC0-2") &
                                                     Normal.DBA.bin50.closed$bin50.closed %in% c(immune.levels[1:6]), ]
Normal.DBA.bin50.closed.iDSC$source <- factor(Normal.DBA.bin50.closed.iDSC$source, levels = c("iDSC0","iDSC1", "iDSC2","DBA.1.iDSC0-2"))
Normal.DBA.bin50.closed.iDSC$bin50.closed <- factor(Normal.DBA.bin50.closed.iDSC$bin50.closed, levels = c(immune.levels[1:6]))
pdf("20220622-Comments/05.DBA/00.figures/41.Normal-DBA-iDSC-4clusters-ajacent-cell-immune-barplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(Normal.DBA.bin50.closed.iDSC, mapping = aes(x = bin50.closed, y = value, fill = bin50.closed), stat = "identity", width = 1) +
  scale_fill_manual(values = c(immune.color[1:6])) +
  facet_grid(sample ~ source, space = "free") +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", color = "black"))
dev.off()

#EC/D6/iDSC 亚群在Normal E85/E95 和 DBA几个样本中的分布
EC.hub <- c("Proliferating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2", "EPC", "iDSC0", "D6-S100a8")
EC.hub.color <- c("#D95F02", "#7B68EE", "#E7298A", "#1B9E77", "#00BFFF", "#00FA9A", "#6A3D9A", "#FDB462")
E85.B6.EC <- E85.B6.bin50.celltype.df.sub[E85.B6.bin50.celltype.df.sub$cell_type %in% EC.hub, ]
E85.B6.EC$cell_type <- factor(E85.B6.EC$cell_type, levels = EC.hub)
pdf("20220622-Comments/05.DBA/00.figures/42.EC-iDSC0-D6-E85.B6-celltype.pdf", width = 5, height = 4.5)
ggplot() + 
  geom_point(E85.B6.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) +
  geom_point(E85.B6.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.3) +
  scale_color_manual(values = c(EC.hub.color[c(1,2,4,5,7,8)])) +
  theme_classic()
dev.off()

E95.EC <- E95.bin50.celltype.df.sub[E95.bin50.celltype.df.sub$cell_type %in% EC.hub, ]
E95.EC$cell_type <- factor(E95.EC$cell_type, levels = EC.hub)
pdf("20220622-Comments/05.DBA/00.figures/42.EC-iDSC0-D6-E95-celltype.pdf", width = 7, height = 6.5)
ggplot() + 
  geom_point(E95.bin50.celltype.df.sub, mapping = aes(x = y, y = x), color = "#E7E7E7", size = 0.5) +
  geom_point(E95.EC, mapping = aes(x = y, y = x, color = cell_type), size = 0.3) +
  scale_color_manual(values = c(EC.hub.color[c(1,2,3,4,5,7,8)])) +
  theme_classic()
dev.off()

DBA.E3.EC <- DBA.E3.bin50.celltype.iDSCsub[DBA.E3.bin50.celltype.iDSCsub$cell_type %in% EC.hub, ]
DBA.E3.EC$cell_type <- factor(DBA.E3.EC$cell_type, levels = EC.hub)
pdf("20220622-Comments/05.DBA/00.figures/42.EC-iDSC0-D6-DBA.E3-celltype.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype.iDSCsub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) +
  geom_point(DBA.E3.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.3) +
  scale_color_manual(values = c(EC.hub.color[c(1,2,4,5,6,7,8)])) +
  theme_classic()
dev.off()

DBA.E6.EC <- DBA.E6.bin50.celltype.iDSCsub[DBA.E6.bin50.celltype.iDSCsub$cell_type %in% EC.hub, ]
DBA.E6.EC$cell_type <- factor(DBA.E6.EC$cell_type, levels = EC.hub)
pdf("20220622-Comments/05.DBA/00.figures/42.EC-iDSC0-D6-DBA.E6-celltype.pdf", width = 7, height = 6.5)
ggplot() + 
  geom_point(DBA.E6.bin50.celltype.iDSCsub, mapping = aes(x = y, y = x), color = "#E7E7E7", size = 0.5) +
  geom_point(DBA.E6.EC, mapping = aes(x = y, y = x, color = cell_type), size = 0.3) +
  scale_color_manual(values = c(EC.hub.color[c(1,2,4,5,6,7,8)])) +
  theme_classic()
dev.off()

DBA.F1.EC <- DBA.F1.bin50.celltype.iDSCsub[DBA.F1.bin50.celltype.iDSCsub$cell_type %in% EC.hub, ]
DBA.F1.EC$cell_type <- factor(DBA.F1.EC$cell_type, levels = EC.hub)
pdf("20220622-Comments/05.DBA/00.figures/42.EC-iDSC0-D6-DBA.F1-celltype.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(DBA.F1.bin50.celltype.iDSCsub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) +
  geom_point(DBA.F1.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.3) +
  scale_color_manual(values = c(EC.hub.color[c(1,2,4,5,6,7,8)])) +
  theme_classic()
dev.off()

#Normal and DBA 空间一起统计免疫和EC hub的细胞比例
DBA.slides.def <- read.table("20220622-Comments/05.DBA/DBA-all-slides-bin50-celltype-df.csv", sep = ",", header = T, check.names = F)
DBA.slides.def <- DBA.slides.def[, c(1,2,3,4,9)]
E85.B6.bin50.celltype.df.sub$sample <- "E85.B6"
E95.bin50.celltype.df.sub$sample <- "E95"
DBA.Normal.E85.E95.df.bin50 <- rbind(DBA.slides.def, E85.B6.bin50.celltype.df.sub, E95.bin50.celltype.df.sub)
immune.cell <- c("Monocyte", "Mac", "Proliferating Mac", "DC-1", "Neutrophil", "NK", "T cell", "B")
DBA.Normal.E85.E95.df.bin50.immune <- DBA.Normal.E85.E95.df.bin50[DBA.Normal.E85.E95.df.bin50$cell_type %in% immune.cell, ]
DBA.Normal.E85.E95.df.bin50.immune.df <- table(DBA.Normal.E85.E95.df.bin50.immune$cell_type, DBA.Normal.E85.E95.df.bin50.immune$sample)
DBA.Normal.E85.E95.df.bin50.immune.df <- DBA.Normal.E85.E95.df.bin50.immune.df[immune.cell, ]
DBA.Normal.E85.E95.df.bin50.immune.df <- apply(DBA.Normal.E85.E95.df.bin50.immune.df, 1, function(x){x/table(DBA.Normal.E85.E95.df.bin50$sample) * 100})
DBA.Normal.E85.E95.df.bin50.immune.df <- melt(DBA.Normal.E85.E95.df.bin50.immune.df)
DBA.Normal.E85.E95.df.bin50.immune.df$zone <- "immune"
#EC
EC.cell <- c("Proliterating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2", "EPC")
DBA.Normal.E85.E95.df.bin50.EC <- DBA.Normal.E85.E95.df.bin50[DBA.Normal.E85.E95.df.bin50$cell_type %in% EC.cell, ]
DBA.Normal.E85.E95.df.bin50.EC.df <- table(DBA.Normal.E85.E95.df.bin50.EC$cell_type, DBA.Normal.E85.E95.df.bin50.EC$sample)
DBA.Normal.E85.E95.df.bin50.EC.df <- DBA.Normal.E85.E95.df.bin50.EC.df[EC.cell, ]
DBA.Normal.E85.E95.df.bin50.EC.df <- apply(DBA.Normal.E85.E95.df.bin50.EC.df, 1, function(x){x/table(DBA.Normal.E85.E95.df.bin50$sample) * 100})
DBA.Normal.E85.E95.df.bin50.EC.df <- melt(DBA.Normal.E85.E95.df.bin50.EC.df)
DBA.Normal.E85.E95.df.bin50.EC.df$zone <- "EC"
#merge all
immune.color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#2E8B57", "#CAB2D6", "#FF7F00", "#FDBF6F")
EC.color <- c("#D95F02", "#7B68EE", "#E7298A", "#1B9E77", "#00BFFF", "#00FA9A")
DBA.Normal.E85.E95.df.bin50.zone.df <- Reduce(rbind, list(DBA.Normal.E85.E95.df.bin50.immune.df, DBA.Normal.E85.E95.df.bin50.EC.df))
DBA.Normal.E85.E95.df.bin50.zone.df$Var1 <- factor(DBA.Normal.E85.E95.df.bin50.zone.df$Var1, levels = c("E85.B6", "F1", "E95", "E3", "E6"))
DBA.Normal.E85.E95.df.bin50.zone.df$Var2 <- factor(DBA.Normal.E85.E95.df.bin50.zone.df$Var2, levels = c(immune.cell, EC.cell))
pdf("20220622-Comments/05.DBA/00.figures/44.DBA-Normal-E85-E95-immune-EC-spatial-proportion-bin50.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(DBA.Normal.E85.E95.df.bin50.zone.df, mapping = aes(x = Var1, y = value, fill = Var2), stat = "identity", position = "stack") + 
  facet_wrap(~ zone, scales = "free") + 
  scale_fill_manual(values = c(immune.color, EC.color)) + 
  theme_classic()
dev.off()

#Immune hub bin20
E85.B6.bin20.celltype.df <- read.table("20210827-figures/figure5/E85-B6-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F, header = T)
E95.bin20.celltype.df <- read.table("20210827-figures/figure5/E95-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F, header = T)
E85.B6.bin20.celltype.df$sample <- "E85.B6"
E95.bin20.celltype.df$sample <- "E95"
DBA.E3.bin20.celltype$sample <- "E3"
DBA.E6.bin20.celltype$sample <- "E6"
DBA.F1.bin20.celltype$sample <- "F1"
DBA.E3.bin20.celltype <- DBA.E3.bin20.celltype[, c(1:4,9)]
DBA.E6.bin20.celltype <- DBA.E6.bin20.celltype[, c(1:4,9)]
DBA.F1.bin20.celltype <- DBA.F1.bin20.celltype[, c(1:4,9)]
Normal.DBA.slides.def.bin20 <- Reduce(rbind, list(E85.B6.bin20.celltype.df, E95.bin20.celltype.df, DBA.E3.bin20.celltype, DBA.E6.bin20.celltype, DBA.F1.bin20.celltype))
Normal.DBA.slides.def.bin20$cell_type <- gsub("Monocyte", "AAAAA", Normal.DBA.slides.def.bin20$cell_type)
Normal.DBA.slides.def.bin20$cell_type <- gsub("Mono", "Monocyte", Normal.DBA.slides.def.bin20$cell_type)
Normal.DBA.slides.def.bin20$cell_type <- gsub("AAAAA", "Monocyte", Normal.DBA.slides.def.bin20$cell_type)
Normal.DBA.slides.def.bin20$cell_type <- gsub("T", "T cell", Normal.DBA.slides.def.bin20$cell_type)
Normal.DBA.slides.def.bin20$cell_type <- gsub("Pro-Mac", "Proliferating Mac", Normal.DBA.slides.def.bin20$cell_type)

Normal.DBA.slides.def.bin20.immune <- Normal.DBA.slides.def.bin20[Normal.DBA.slides.def.bin20$cell_type %in% immune.cell, ]
Normal.DBA.slides.def.bin20.immune.df <- table(Normal.DBA.slides.def.bin20.immune$cell_type, Normal.DBA.slides.def.bin20.immune$sample)
Normal.DBA.slides.def.bin20.immune.df <- Normal.DBA.slides.def.bin20.immune.df[immune.cell, ]
Normal.DBA.slides.def.bin20.immune.df <- apply(Normal.DBA.slides.def.bin20.immune.df, 1, function(x){x/table(Normal.DBA.slides.def.bin20$sample) * 100})
Normal.DBA.slides.def.bin20.immune.df <- melt(Normal.DBA.slides.def.bin20.immune.df)
Normal.DBA.slides.def.bin20.immune.df$zone <- "immune"
Normal.DBA.slides.def.bin20.immune.df$Var1 <- factor(Normal.DBA.slides.def.bin20.immune.df$Var1, levels = c("E85.B6", "F1", "E95", "E3", "E6"))
Normal.DBA.slides.def.bin20.immune.df$Var2 <- factor(Normal.DBA.slides.def.bin20.immune.df$Var2, levels = immune.cell)
pdf("20220622-Comments/05.DBA/00.figures/44.DBA-Normal-E85-E95-immune-spatial-proportion-bin20.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(Normal.DBA.slides.def.bin20.immune.df, mapping = aes(x = Var1, y = value, fill = Var2), stat = "identity", position = "stack") + 
  facet_wrap(~ zone, scales = "free") + 
  scale_fill_manual(values = c(immune.color)) + 
  theme_classic()
dev.off()
  
#Quality
load("20220622-Comments/05.DBA/DBA-E3-data.Rdata")
load("20220622-Comments/05.DBA/DBA-E6-data.Rdata")
load("20220622-Comments/05.DBA/DBA-F1-data.Rdata")
pdf("20220622-Comments/05.DBA/00.figures/56.QC-CBA-E3.pdf", width = 10, height = 10)
SpatialFeaturePlot(DBA.E3.bin50, features = c("nFeature_Spatial", "nCount_Spatial"), stroke = NA)
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/56.QC-CBA-E6.pdf", width = 10, height = 10)
SpatialFeaturePlot(DBA.E6.bin50, features = c("nFeature_Spatial", "nCount_Spatial"), stroke = NA)
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/56.QC-CBA-F1.pdf", width = 10, height = 10)
SpatialFeaturePlot(DBA.F1.bin50, features = c("nFeature_Spatial", "nCount_Spatial"), stroke = NA)
dev.off()

#all slide prediction score
load("20220622-Comments/05.DBA/DBA-E3-data-map-embryo-20221216.RData")
DBA.E3.bin50.embryo <- subset(DBA.E3.bin50.embryo, features = c(rownames(DBA.E3.bin50.embryo)[1:12]))
load("20220622-Comments/05.DBA/DBA-E6-data-map-embryo-20221216.RData")
DBA.E6.bin50.embryo <- subset(DBA.E6.bin50.embryo, features = c(rownames(DBA.E6.bin50.embryo)[1:12]))
load("20220622-Comments/05.DBA/DBA-F1-data-map-embryo-20221216.RData")
DBA.F1.bin50.embryo <- subset(DBA.F1.bin50.embryo, features = c(rownames(DBA.F1.bin50.embryo)[1:12]))
load("20220622-Comments/05.DBA/DBA-E3-data-map-maternal-20221216.RData")
DBA.E3.bin50.maternal <- subset(DBA.E3.bin50.maternal, features = c(rownames(DBA.E3.bin50.maternal)[1:23]))
load("20220622-Comments/05.DBA/DBA-E6-data-map-maternal-20221216.RData")
DBA.E6.bin50.maternal <- subset(DBA.E6.bin50.maternal, features = c(rownames(DBA.E6.bin50.maternal)[1:23]))
load("20220622-Comments/05.DBA/DBA-F1-data-map-maternal-20221216.RData")
DBA.F1.bin50.maternal <- subset(DBA.F1.bin50.maternal, features = c(rownames(DBA.F1.bin50.maternal)[1:23]))

DBA.slide.score <- merge(DBA.E3.bin50.embryo, y = list(DBA.E6.bin50.embryo, DBA.F1.bin50.embryo, DBA.E3.bin50.maternal, DBA.E6.bin50.maternal, DBA.F1.bin50.maternal),
                           add.cell.ids = c("DBA.E3.bin50.embryo", "DBA.E6.bin50.embryo", "DBA.F1.bin50.embryo", "DBA.E3.bin50.maternal", "DBA.E6.bin50.maternal", "DBA.F1.bin50.maternal"))
DBA.slide.score.levels <- c(decidual.levels, immune.levels[c(8,1:7)], EC.levels[c(1,3,2,4,5)], major.levels[c(4,5,11)], tropho.levels[c(3,2,5,1,4)])
#celltype.color <- c(decidual.color, immune.color[c(8,1:7)], EC.color, major.color[c(4,5,11)], tropho.color[c(3,2,5,1,4)])
annotation.color <- list(cell_type = c("D1-eSF"="#6A5ACD", "D2-Pre-DSC"="#9ACD32", "D3-Top2a"="#FB8072","D4-Ptn"="#8DD3C7","D5-Gatm"="#08519C",
                                       "D6-S100a8"="#FDB462","D7-Prl8a2"="#228B22", "Monocyte"="#A6CEE3",
                                       "Mac"="#1F78B4","Proliferating Mac"="#B2DF8A","DC-1"="#33A02C", "Neutrophil"="#2E8B57", "NK"="#CAB2D6",
                                       "T cell"="#FDBF6F","iDSC-DBA1"="#6A3D9A","Angiogenic EC"="#7B68EE", "Venous EC1"="#1B9E77",
                                       "Proliterating EC"="#D95F02","Venous EC2"="#00BFFF", "EPC"="#00FA9A", "SMC"="#6B8E23","Epithelial"="#006400",
                                       "Mesenchymal stem cell"="#4B0082","Mesenchymal-epithelial transition"="#CD5C5C",
                                       "Blood P"="#DAA520", "Viseral endoderm"="#9370DB", "Allantois mesodermal"="#8B4513",
                                       "Fetal EC" = "#FFD700", "tropho-0"="#3CB371",
                                       "tropho-1"="#DA70D6", "tropho-2"="#4682B4","tropho-3"="#483D8B","tropho-4"="#008B8B"))

DBA.slide.score.celltype.score <- GetAssayData(DBA.slide.score, slot = "data")
DBA.slide.score.celltype.max <- apply(DBA.slide.score.celltype.score, 2, 
                                        function(x) rownames(DBA.slide.score.celltype.score)[which.max(x)])
DBA.slide.score.celltype.max.sub <- DBA.slide.score.celltype.max[DBA.slide.score.celltype.max %in% DBA.slide.score.levels]

DBA.slide.score.celltype.sort <- DBA.slide.score.celltype.score[, names(sort(DBA.slide.score.celltype.max.sub))]
col_anno <- data.frame(cell_type = DBA.slide.score.celltype.max.sub)
rownames(col_anno) <- names(DBA.slide.score.celltype.max.sub)
col_anno$cell_type <- factor(col_anno$cell_type, levels = DBA.slide.score.levels)
col_anno$cell_name <- rownames(col_anno)
col_anno.sort <- data.frame()
for (n in DBA.slide.score.levels) {
  col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
}
col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type)
rownames(col_anno.df) <- rownames(col_anno.sort)
color <- brewer.pal(n =11, name = "RdBu")
color <- color[c(1:4, 7:11)]
color = colorRampPalette(rev(color))(100)
DBA.slide.score.celltype.sort.s <- DBA.slide.score.celltype.sort[DBA.slide.score.levels, rownames(col_anno.df)]
pdf("20220622-Comments/05.DBA/00.figures/57.CBA-all-slides-score-heatmap.pdf", width = 20, height = 10)
pheatmap(DBA.slide.score.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
         scale = "column", show_colnames = F, annotation_col = col_anno.df, 
         color = color)
dev.off()

#spatial clustering
pdf("20220622-Comments/05.DBA/00.figures/58.DBA-E3-bin50-clustering-spatial.pdf", width = 6, height = 7)
SpatialDimPlot(DBA.E3.bin50, stroke = NA)
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/58.DBA-E6-bin50-clustering-spatial.pdf", width = 8, height = 7)
SpatialDimPlot(DBA.E6.bin50, stroke = NA)
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/58.DBA-F1-bin50-clustering-spatial.pdf", width = 7, height = 5)
SpatialDimPlot(DBA.F1.bin50, stroke = NA)
dev.off()
  
DBA.E3.bin50.markers <- FindAllMarkers(DBA.E3.bin50, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.E3.bin50.markers.top20 <- DBA.E3.bin50.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
write.table(DBA.E3.bin50.markers.top20, file = "20220622-Comments/05.DBA/DBA-E3-spatial-cluster-DEGs-top20.csv", sep = ",", quote = F)

DBA.E6.bin50.markers <- FindAllMarkers(DBA.E6.bin50, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.E6.bin50.markers.top20 <- DBA.E6.bin50.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
write.table(DBA.E6.bin50.markers.top20, file = "20220622-Comments/05.DBA/DBA-E6-spatial-cluster-DEGs-top20.csv", sep = ",", quote = F)

DBA.F1.bin50.markers <- FindAllMarkers(DBA.F1.bin50, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.F1.bin50.markers.top20 <- DBA.F1.bin50.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
write.table(DBA.F1.bin50.markers.top20, file = "20220622-Comments/05.DBA/DBA-F1-spatial-cluster-DEGs-top20.csv", sep = ",", quote = F)
