#Author: Min Yang

library("Seurat") #version > 4
library("ggplot2")

#DSC integration, keep the same UMAP reduction as the normal samples, showed in figure 1
load("DSC.Normal.RData")
load("DSC.CBA.RData")
DSC.Normal.UMAP_1 <- DSC.Normal@reductions$umap@cell.embeddings[, 1]
DSC.Normal.UMAP_2 <- DSC.Normal@reductions$umap@cell.embeddings[, 2]
DSC.Normal <- ScaleData(DSC.Normal, features = rownames(DSC.Normal))
DSC.Normal <- RunPCA(DSC.Normal, npcs = 50)
DSC.Normal <- RunUMAP(DSC.Normal, dims = 1:30, return.model = TRUE)
DSC.Normal@reductions$umap@misc$model$embedding[, 1] <- DSC.Normal.UMAP_1
DSC.Normal@reductions$umap@misc$model$embedding[, 2] <- DSC.Normal.UMAP_2
DSC.Normal@reductions$umap@cell.embeddings[, 1] <- DSC.Normal.UMAP_1
DSC.Normal@reductions$umap@cell.embeddings[, 2] <- DSC.Normal.UMAP_2

CBA.DSC <- NormalizeData(CBA.DSC)
CBA.DSC <- FindVariableFeatures(CBA.DSC, selection.method = "vst", nfeatures = 2000)
DSC.anchors <- FindTransferAnchors(reference = UE.decidual.susbet, query = CBA.DSC, dims = 1:30, reference.reduction = "pca")
DSC.predictions <- TransferData(anchorset = DSC.anchors, refdata = UE.decidual.susbet$sub_clusters, dims = 1:30)
CBA.DSC <- AddMetaData(CBA.DSC, metadata = DSC.predictions)
CBA.DSC <- MapQuery(anchorset = DSC.anchors, reference = UE.decidual.susbet, query = CBA.DSC, 
                      refdata = list(celltype = "sub_clusters"), reference.reduction = "pca", reduction.model = "umap")
Normal.DSC.umap <- data.frame(UMAP_1 = DSC.Normal.UMAP_1, 
                              UMAP_2 = DSC.Normal.UMAP_2,
                              celltype = DSC.Normal$subname)
CBA.DSC.umap <- data.frame(UMAP_1 = CBA.DSC@reductions$ref.umap@cell.embeddings[, 1], 
                             UMAP_2 = CBA.DSC@reductions$ref.umap@cell.embeddings[, 2], 
                             celltype = CBA.DSC$subname)
ggplot() + 
  geom_point(Normal.DSC.umap, mapping = aes(x = UMAP_1, y = UMAP_2, color = celltype), alpha = 0.1, shape = 4, zie = 1.5) + 
  geom_point(CBA.DSC.umap, mapping = aes(x = UMAP_1, y = UMAP_2, color = celltype), shape = 16, size = 1) +
  scale_color_manual(values = color) +
  theme_classic()

#The same strategy used in immune cell/EC/iDSC integration






