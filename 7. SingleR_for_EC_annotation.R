#Author: Min Yang


# SingleR annoatation
library("SingleR")
cell.ref <- read.csv("Cell-EC-normalize/Data.csv", header = T, row.names = 1)
cell.ref.t <- data.frame(t(as.matrix(cell.ref)), check.names = F)
cell.ref.t$Observation <- rownames(cell.ref.t)
cell.ref.meta <- read.csv("Cell-EC-normalize/Metadata.csv")
cell.ref.m <- merge(cell.ref.meta, cell.ref.t, by = "Observation", all = T)
cell.ref.m <- aggregate(cell.ref.m[, -c(1:3)], by = list(cell.ref.m$Cluster), mean)
rownames(cell.ref.m) <- cell.ref.m$Group.1
cell.ref.m <- cell.ref.m[, -1]
cell.ref.m <- t(cell.ref.m)
cell.ref.sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts=cell.ref.m))
colData(cell.ref.sce)$Type <- colnames(cell.ref.sce)
load("EC.seurat.RData")
EC.seurat.data <- GetAssayData(EC.seurat, slot = "data")
EC.seurat.pred <- SingleR(test = EC.seurat.data, ref = cell.EC.sce, labels = cell.EC.sce$Type)
save(EC.seurat.pred, file = "SingleR-EC-prediction.RData")
table(EC.seurat.pred$labels, EC.seurat@active.ident)
