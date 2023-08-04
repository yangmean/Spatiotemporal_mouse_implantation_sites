library(tidyr)
load("../20210201-integrate/UE-harmony-pc30-UMAP-subname.RData")
DimPlot(UE.harmony.pc30)
#for single cell: UE.harmony.pc30
qc.count.gene.df <- data.frame(nCount = UE.harmony.pc30$nCount_RNA, 
                               nFeature = UE.harmony.pc30$nFeature_RNA, 
                               sample = "SC")
E65.1.data <- read.table("/Data3/yangmin/Jennie/20210728-spatial/R_bin_remove/E6.5-1.bin1.gem", sep = "\t", header = T, check.names = F)
E65.2.data <- read.table("/Data3/yangmin/Jennie/20210728-spatial/R_bin_remove/E6.5-2.bin1.gem", sep = "\t", header = T, check.names = F)
E75.1.data <- read.table("/Data3/yangmin/Jennie/20210728-spatial/R_bin_remove/E7.5-1.bin1.gem", sep = "\t", header = T, check.names = F)
E75.2.data <- read.table("/Data3/yangmin/Jennie/20210728-spatial/R_bin_remove/E7.5-2.bin1.gem", sep = "\t", header = T, check.names = F)
B6.data <- read.table("/Data3/yangmin/Jennie/20210728-spatial/R_bin_remove/FP200000337BR_B6.bin1.Lasso.gem", sep = "\t", header = T, check.names = F)
C1.data <- read.table("/Data3/yangmin/Jennie/20210728-spatial/R_bin_remove/FP200000335BR_C1.bin1.Lasso.gem", sep = "\t", header = T, check.names = F)
E95.data <- read.table("/Data3/yangmin/Jennie/20210728-spatial/R_bin_remove/E9.5-bin10.final.gem", sep = "\t", header = T, check.names = F)

#function
CreateStereoObject <- function(data, bin){
  temp <- data
  temp$x <- trunc(temp$x/bin) * bin
  temp$y <- trunc(temp$y/bin) * bin
  temp$cellID <- paste(temp$x, "_", temp$y, sep = "")
  temp <- aggregate(temp$MIDCounts, by = list(temp$cellID, temp$geneID), sum)
  colnames(temp) <- c("cellID", "geneID", "MIDCounts")
  temp$cellInx <- match(temp$cellID, unique(temp$cellID))
  temp$geneInx <- match(temp$geneID, unique(temp$geneID))
  mat <- sparseMatrix(i = temp$geneInx, j = temp$cellInx, x = temp$MIDCounts, 
                      dimnames = list(unique(temp$geneID), unique(temp$cellID)))
  temp.coord.df <- data.frame(cellname = colnames(mat))
  rownames(temp.coord.df) <- temp.coord.df$cellname
  temp.coord.df <- separate(temp.coord.df, col = cellname, sep = "_", into = c("x", "y"))
  temp <- CreateSeuratObject(mat, project = "temp", assay = "Spatial")
  temp$slice <- 1
  temp$region <- "temp"
  colnames(temp.coord.df) <- c("imagerow", "imagecol")
  temp.coord.df$imagerow <- as.numeric(temp.coord.df$imagerow)
  temp.coord.df$imagecol <- as.numeric(temp.coord.df$imagecol)
  temp@images$spatial <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = temp.coord.df)
  return(temp)
}
#bin10
E65.1.bin10 <- CreateStereoObject(E65.1.data, bin = 10)
E65.2.bin10 <- CreateStereoObject(E65.2.data, bin = 10)
E75.1.bin10 <- CreateStereoObject(E75.1.data, bin = 10)
E75.2.bin10 <- CreateStereoObject(E75.2.data, bin = 10)
B6.bin10 <- CreateStereoObject(B6.data, bin = 10)
C1.bin10 <- CreateStereoObject(C1.data, bin = 10)
E95.bin10 <- CreateStereoObject(E95.data, bin = 10)
save(E65.1.bin10, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin10.RData")
save(E65.2.bin10, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.2.bin10.RData")
save(E75.1.bin10, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.1.bin10.RData")
save(E75.2.bin10, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.2.bin10.RData")
save(B6.bin10, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/B6.bin10.RData")
save(C1.bin10, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/C1.bin10.RData")
save(E95.bin10, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95.bin10.RData")

load("20210818-spatial/UE-harmony-pc30-UMAP-subname-SCT.RData")
UE.harmony.SCT <- subset(UE.harmony.SCT, subset = subname != "cluster9-1")
IntegrationSC <- function(ST.data, sc.data){
  ST.data <- SCTransform(ST.data, assay = "Spatial", verbose = FALSE, ncells = 5000)
  ST.data <- RunPCA(ST.data, assay = "SCT", features = VariableFeatures(ST.data), verbose = F)
  ST.data <- FindNeighbors(ST.data, reduction = "pca", dims = 1:30, verbose = F)
  ST.data <- FindClusters(ST.data, resolution = 0.3, verbose = F)
  ST.data <- RunUMAP(ST.data, reduction = "pca", dims = 1:30, verbose = F)
  anchors <- FindTransferAnchors(reference = sc.data, query = ST.data, normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, refdata = sc.data$subname, 
                                    prediction.assay = TRUE, weight.reduction = ST.data[["pca"]], dims = 1:30)
  ST.data[["predictions"]] <- predictions.assay
  DefaultAssay(ST.data) <- "predictions"
  return(ST.data)
}

load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin10.RData")
E65.1.bin10.anno <- IntegrationSC(E65.1.bin10, UE.harmony.SCT)
save(E65.1.bin10.anno, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin10.map.RData")
rm(E65.1.bin10)
rm(E65.1.bin10.anno)
gc()

load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.2.bin10.RData")
E65.2.bin10.anno <- IntegrationSC(E65.2.bin10, UE.harmony.SCT)
save(E65.2.bin10.anno, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.2.bin10.map.RData")
rm(E65.2.bin10)
rm(E65.2.bin10.anno)
gc()

load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.1.bin10.RData")
E75.1.bin10.anno <- IntegrationSC(E75.1.bin10, UE.harmony.SCT)
save(E75.1.bin10.anno, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.1.bin10.map.RData")
rm(E75.1.bin10)
rm(E75.1.bin10.anno)
gc()

load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.2.bin10.RData")
E75.2.bin10.anno <- IntegrationSC(E75.2.bin10, UE.harmony.SCT)
save(E75.2.bin10.anno, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.2.bin10.map.RData")
rm(E75.2.bin10)
rm(E75.2.bin10.anno)
gc()

load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/B6.bin10.RData")
B6.bin10.anno <- IntegrationSC(B6.bin10, UE.harmony.SCT)
save(B6.bin10.anno, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/B6.bin10.map.RData")
rm(B6.bin10)
rm(B6.bin10.anno)
gc()

load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/C1.bin10.RData")
C1.bin10.anno <- IntegrationSC(C1.bin10, UE.harmony.SCT)
save(C1.bin10.anno, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/C1.bin10.map.RData")
rm(C1.bin10)
rm(C1.bin10.anno)
gc()

load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95.bin10.RData")
E95.bin10.anno <- IntegrationSC(E95.bin10, UE.harmony.SCT)
save(E95.bin10.anno, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95.bin10.map.RData")
rm(E95.bin10)
rm(E95.bin10.anno)
gc()

E65.1.bin100.anno <- IntegrationSC(E65.1.bin100, UE.harmony.SCT)
save(E65.1.bin100.anno, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin100.map.RData")
rm(E65.1.bin100)
rm(E65.1.bin100.anno)
gc()

#bin20
E65.1.bin20 <- CreateStereoObject(E65.1.data, bin = 20)
save(E65.1.bin20, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin20.RData")
E65.2.bin20 <- CreateStereoObject(E65.2.data, bin = 20)
save(E65.2.bin20, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.2.bin20.RData")
E75.1.bin20 <- CreateStereoObject(E75.1.data, bin = 20)
save(E75.1.bin20, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.1.bin20.RData")
E75.2.bin20 <- CreateStereoObject(E75.2.data, bin = 20)
save(E75.2.bin20, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.2.bin20.RData")
B6.bin20 <- CreateStereoObject(B6.data, bin = 20)
save(B6.bin20, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/B6.bin20.RData")
C1.bin20 <- CreateStereoObject(C1.data, bin = 20)
save(C1.bin20, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/C1.bin20.RData")
E95.bin20 <- CreateStereoObject(E95.data, bin = 20)
save(E95.bin20, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95.bin20.RData")

#bin100
E65.1.bin100 <- CreateStereoObject(E65.1.data, bin = 100)
save(E65.1.bin100, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin100.RData")
E65.2.bin100 <- CreateStereoObject(E65.2.data, bin = 100)
save(E65.2.bin100, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.2.bin100.RData")
E75.1.bin100 <- CreateStereoObject(E75.1.data, bin = 100)
save(E75.1.bin100, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.1.bin100.RData")
E75.2.bin100 <- CreateStereoObject(E75.2.data, bin = 100)
save(E75.2.bin100, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.2.bin100.RData")
B6.bin100 <- CreateStereoObject(B6.data, bin = 100)
save(B6.bin100, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/B6.bin100.RData")
C1.bin100 <- CreateStereoObject(C1.data, bin = 100)
save(C1.bin100, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/C1.bin100.RData")
E95.bin100 <- CreateStereoObject(E95.data, bin = 100)
save(E95.bin100, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95.bin100.RData")

#nCount/nGene
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65-1-bin50-spatial-phago-map-remove.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65-2-bin50-spatial-phago-map-remove.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75-1-bin50-spatial-phago-map-remove.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75-2-bin50-spatial-phago-map-remove.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/B6-bin50-spatial-phago-map.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/C1-bin50-spatial-phago-map.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95-bin50-spatial-phago-map-remove.RData")
Metadata.fun <- function(object){
  DefaultAssay(object) <-  "predictions"
  object.celltype <- GetAssayData(object, slot = "data")
  object.celltype <- apply(object.celltype, 2, function(x) rownames(object.celltype)[which.max(x)])
  object.celltype <- data.frame(cell_type = object.celltype, 
                                cell_name = names(object.celltype),
                                nFeature_RNA = object$nFeature_Spatial, 
                                nCount_RNA = object$nCount_Spatial)
}

UE.harmony.SCT.celltype <- data.frame(cell_type = UE.harmony.SCT$subname, 
                                      cell_name = colnames(UE.harmony.SCT),
                                      nFeature_RNA = UE.harmony.SCT$nFeature_RNA, 
                                      nCount_RNA = UE.harmony.SCT$nCount_RNA)
UE.harmony.SCT.celltype$sample <- "sc_RNA"
E65.1.bin50.anno.celltype <- Metadata.fun(E65.1.bin50)
E65.1.bin50.anno.celltype$sample <- "E65.1"
E65.2.bin50.anno.celltype <- Metadata.fun(E65.2.bin50)
E65.2.bin50.anno.celltype$sample <- "E65.2"
E75.1.bin50.anno.celltype <- Metadata.fun(E75.1.bin50)
E75.1.bin50.anno.celltype$sample <- "E75.1"
E75.2.bin50.anno.celltype <- Metadata.fun(E75.2.bin50)
E75.2.bin50.anno.celltype$sample <- "E75.2"
B6.bin50.anno.celltype <- Metadata.fun(B6.bin50)
B6.bin50.anno.celltype$sample <- "B6"
C1.bin50.anno.celltype <- Metadata.fun(C1.bin50)
C1.bin50.anno.celltype$sample <- "C1"
E95.bin50.anno.celltype <- Metadata.fun(E95.bin50)
E95.bin50.anno.celltype$sample <- "E95"
bin50.celltype.QC <- rbind(UE.harmony.SCT.celltype, E65.1.bin50.anno.celltype, E65.2.bin50.anno.celltype, E75.1.bin50.anno.celltype, 
                           E75.2.bin50.anno.celltype, B6.bin50.anno.celltype, C1.bin50.anno.celltype, E95.bin50.anno.celltype)
bin50.celltype.QC$Count_Feature <- bin50.celltype.QC$nCount_RNA/bin50.celltype.QC$nFeature_RNA
save(bin50.celltype.QC, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/bin50.celltype.QC.RData")
bin50.celltype.QC$sample <- factor(bin50.celltype.QC$sample, levels = c("sc_RNA", "E65.1", "E65.2", "E75.1", "E75.2", "B6", "C1", "E95"))
pdf("20220622-Comments/06.ST_QC/01.scRNA_ST_Count_Feature_celltype_violin.pdf", width = 20, height = 10)
ggplot() + 
  geom_violin(bin50.celltype.QC, mapping = aes(x = sample, y = Count_Feature, fill = sample), scale = "width") + 
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_fill_manual(values = c("darkred", "#AE6156", "#AE6156", "#D6987F", "#D6987F", "#8DB1C0", "#8DB1C0", "#5882A0", "#5882A0")) + 
  theme_bw()
dev.off()

#score_bin10
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin10.map.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.2.bin10.map.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.1.bin10.map.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.2.bin10.map.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/B6.bin10.map.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/C1.bin10.map.RData")
load("/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95.bin10.map.RData")
all.slide.score.2 <- merge(E65.1.bin10.anno, list(E65.2.bin10.anno, E75.1.bin10.anno, E75.2.bin10.anno, B6.bin10.anno, C1.bin10.anno, E95.bin10.anno), 
                           add.cell.ids = c("E65.1.bin10", "E65.2.bin10", "E75.1.bin10", "E75.2.bin10", "E85.B6.bin10", "E85.C1.bin10", "E95.bin50"))
all.slide.score.2.levels <- c(decidual.levels[c(1:7)], immune.levels[10], immune.levels[1:9], EC.levels, major.levels[c(4,5,7,11)], tropho.levels[c(0,3,2,5,1,4)])
all.slide.score.2.celltype.score <- GetAssayData(all.slide.score.2, slot = "data")
all.slide.score.2.celltype.score <- all.slide.score.2.celltype.score[-27, ]
all.slide.score.2.celltype.max <- apply(all.slide.score.2.celltype.score, 2, 
                                        function(x) rownames(all.slide.score.2.celltype.score)[which.max(x)])
all.slide.score.2.celltype.sort <- all.slide.score.2.celltype.score[, names(sort(all.slide.score.2.celltype.max))]
col_anno <- data.frame(cell_type = all.slide.score.2.celltype.max)
rownames(col_anno) <- names(all.slide.score.2.celltype.max)
col_anno$cell_type <- factor(col_anno$cell_type, levels = all.slide.score.2.levels)
col_anno$cell_name <- rownames(col_anno)
col_anno.sort <- data.frame()
for (n in all.slide.score.2.levels) {
  col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
}
col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type, 
                          time = matrix(unlist(strsplit(col_anno.sort$cell_name, "_")), ncol = 3, byrow = T)[, 1])
rownames(col_anno.df) <- rownames(col_anno.sort)
col_anno.df$time <- factor(col_anno.df$time, levels = c("E65.2.bin50", "E65.1.bin50", "E75.2.bin50", "E75.1.bin50", "E85.B6.bin50", "E85.C1.bin50", "E95.bin50"))
color <- brewer.pal(n =11, name = "RdBu")
color <- color[c(1:4, 7:11)]
color = colorRampPalette(rev(color))(100)
all.slide.score.2.celltype.sort.s <- all.slide.score.2.celltype.sort[all.slide.score.2.levels, rownames(col_anno.df)]
pdf("20210827-figures/figure2-spatial/Figure2-all-2-seurat-score-heatmap.pdf", width = 20, height = 10)
pheatmap(all.slide.score.2.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
         scale = "column", show_colnames = F, annotation_col = col_anno.df, 
         color = color)
dev.off()

########################################20221017 calculate cell type proportion in ST###########################################
E65.1.bin50.celltype.df <- read.table("20210818-spatial/E65-1-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E65.1.bin50.celltype.df <- E65.1.bin50.celltype.df[E65.1.bin50.celltype.df$cell_type != "Erythroid", ]
E65.1.bin50.celltype.df$sample <- "E65.1"
E65.2.bin50.celltype.df <- read.table("20210818-spatial/E65-2-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E65.2.bin50.celltype.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type != "Erythroid", ]
E65.2.bin50.celltype.df$sample <- "E65.2"
E75.1.bin50.celltype.df <- read.table("20210818-spatial/E75-1-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E75.1.bin50.celltype.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type != "Erythroid", ]
E75.1.bin50.celltype.df$sample <- "E75.1"
E75.2.bin50.celltype.df <- read.table("20210818-spatial/E75-2-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E75.2.bin50.celltype.df <- E75.2.bin50.celltype.df[E75.2.bin50.celltype.df$cell_type != "Erythroid", ]
E75.2.bin50.celltype.df$sample <- "E75.2"
B6.bin50.celltype.df <- read.table("20210818-spatial/E85-B6-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
B6.bin50.celltype.df <- B6.bin50.celltype.df[B6.bin50.celltype.df$cell_type != "Erythroid", ]
B6.bin50.celltype.df$sample <- "E85.B6"
C1.bin50.celltype.df <- read.table("20210818-spatial/E85-C1-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
C1.bin50.celltype.df <- C1.bin50.celltype.df[C1.bin50.celltype.df$cell_type != "Erythroid", ]
C1.bin50.celltype.df$sample <- "E85.C1"
E95.bin50.celltype.df <- read.table("20210818-spatial/E95-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E95.bin50.celltype.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type != "Erythroid", ]
E95.bin50.celltype.df$sample <- "E95"

All.ST.df <- rbind(E65.1.bin50.celltype.df, E65.2.bin50.celltype.df, E75.1.bin50.celltype.df, E75.2.bin50.celltype.df, 
                   B6.bin50.celltype.df, C1.bin50.celltype.df, E95.bin50.celltype.df)
save(All.ST.df, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/All-ST-dataframe.RData")
All.ST.cell.df <- data.frame(table(All.ST.df$cell_type, All.ST.df$sample), check.names = F)
All.ST.cell.df$celltype <- rownames(All.ST.cell.df)
cell_type.levels <- c(decidual.levels[c(1:7)], immune.levels[c(1:10)], EC.levels, major.levels[c(4:5,7,11)], tropho.levels[c(1:5)], "embryo")
color <- c(decidual.color[c(1:7)], immune.color[c(1:10)], EC.color, major.color[c(4:5,7,11)], tropho.color[c(1:5)], "lightgray")
All.ST.cell.df$Var1 <- factor(All.ST.cell.df$Var1, levels = cell_type.levels)
pdf("20220622-Comments/06.ST_QC/02.ST_cell_type_percentage.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(All.ST.cell.df, mapping = aes(x = Var2, y = Freq, fill = Var1), stat = "identity", position = "fill") + 
  scale_fill_manual(values = color) + 
  theme(panel.background = element_rect(color = "black", fill = "white"), 
        panel.grid = element_blank())
dev.off()

#bin3/14/140 UMI boxplot compare to other ST methods
E65.1.bin3 <- CreateStereoObject(E65.1.data, bin = 3)
E65.2.bin3 <- CreateStereoObject(E65.2.data, bin = 3)
E75.1.bin3 <- CreateStereoObject(E75.1.data, bin = 3)
E75.2.bin3 <- CreateStereoObject(E75.2.data, bin = 3)
B6.bin3 <- CreateStereoObject(B6.data, bin = 3)
C1.bin3 <- CreateStereoObject(C1.data, bin = 3)
E95.bin3 <- CreateStereoObject(E95.data, bin = 3)

E65.1.bin14 <- CreateStereoObject(E65.1.data, bin = 14)
E65.2.bin14 <- CreateStereoObject(E65.2.data, bin = 14)
E75.1.bin14 <- CreateStereoObject(E75.1.data, bin = 14)
E75.2.bin14 <- CreateStereoObject(E75.2.data, bin = 14)
B6.bin14 <- CreateStereoObject(B6.data, bin = 14)
C1.bin14 <- CreateStereoObject(C1.data, bin = 14)
E95.bin14 <- CreateStereoObject(E95.data, bin = 14)

E65.1.bin140 <- CreateStereoObject(E65.1.data, bin = 140)
E65.2.bin140 <- CreateStereoObject(E65.2.data, bin = 140)
E75.1.bin140 <- CreateStereoObject(E75.1.data, bin = 140)
E75.2.bin140 <- CreateStereoObject(E75.2.data, bin = 140)
B6.bin140 <- CreateStereoObject(B6.data, bin = 140)
C1.bin140 <- CreateStereoObject(C1.data, bin = 140)
E95.bin140 <- CreateStereoObject(E95.data, bin = 140)

save(E65.1.bin3, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin3.RData")
save(E65.2.bin3, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.2.bin3.RData")
save(E75.1.bin3, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.1.bin3.RData")
save(E75.2.bin3, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.2.bin3.RData")
save(B6.bin3, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E85.B6.bin3.RData")
save(C1.bin3, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E85.C1.bin3.RData")
save(E95.bin3, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95.bin3.RData")

save(E65.1.bin14, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin14.RData")
save(E65.2.bin14, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.2.bin14.RData")
save(E75.1.bin14, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.1.bin14.RData")
save(E75.2.bin14, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.2.bin14.RData")
save(B6.bin14, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E85.B6.bin14.RData")
save(C1.bin14, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E85.C1.bin14.RData")
save(E95.bin14, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95.bin14.RData")

save(E65.1.bin140, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.1.bin140.RData")
save(E65.2.bin140, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E65.2.bin140.RData")
save(E75.1.bin140, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.1.bin140.RData")
save(E75.2.bin140, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E75.2.bin140.RData")
save(B6.bin140, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E85.B6.bin140.RData")
save(C1.bin140, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E85.C1.bin140.RData")
save(E95.bin140, file = "/Data2/Users/yangmin/Jennie/20200801-SC/R/20220622-Comments/06.ST_QC/E95.bin140.RData")

#get UMI data
UMI.bin3.df <- data.frame(UMI = c(E65.1.bin3$nCount_Spatial, E65.2.bin3$nCount_Spatial, E75.1.bin3$nCount_Spatial, E75.2.bin3$nCount_Spatial, 
                             B6.bin3$nCount_Spatial, C1.bin3$nCount_Spatial, E95.bin10$nCount_Spatial/11), 
                     sample = c(rep("E65.1", dim(E65.1.bin3)[2]), rep("E65.2", dim(E65.2.bin3)[2]), rep("E75.1", dim(E75.1.bin3)[2]), rep("E75.2", dim(E75.2.bin3)[2]), 
                                rep("E85.B6", dim(B6.bin3)[2]), rep("E85.C1", dim(C1.bin3)[2]), rep("E95", dim(E95.bin10)[2])))
UMI.bin3.df$bin <- "bin3"
UMI.bin3.df$sample <- factor(UMI.bin3.df$sample, levels = c("E65.2","E65.1","E75.1","E75.2","E85.B6","E85.C1","E95"))

UMI.bin14.df <- data.frame(UMI = c(E65.1.bin14$nCount_Spatial, E65.2.bin14$nCount_Spatial, E75.1.bin14$nCount_Spatial, E75.2.bin14$nCount_Spatial, 
                                  B6.bin14$nCount_Spatial, C1.bin14$nCount_Spatial, E95.bin10$nCount_Spatial * 1.96), 
                                  sample = c(rep("E65.1", dim(E65.1.bin14)[2]), rep("E65.2", dim(E65.2.bin14)[2]), rep("E75.1", dim(E75.1.bin14)[2]), rep("E75.2", dim(E75.2.bin14)[2]), 
                                             rep("E85.B6", dim(B6.bin14)[2]), rep("E85.C1", dim(C1.bin14)[2]), rep("E95", dim(E95.bin10)[2])))
UMI.bin14.df$bin <- "bin14"
UMI.bin14.df$sample <- factor(UMI.bin14.df$sample, levels = c("E65.2","E65.1","E75.1","E75.2","E85.B6","E85.C1","E95"))

UMI.bin140.df <- data.frame(UMI = c(E65.1.bin140$nCount_Spatial, E65.2.bin140$nCount_Spatial, E75.1.bin140$nCount_Spatial, E75.2.bin140$nCount_Spatial, 
                                  B6.bin140$nCount_Spatial, C1.bin140$nCount_Spatial, E95.bin140$nCount_Spatial), 
                                  sample = c(rep("E65.1", dim(E65.1.bin140)[2]), rep("E65.2", dim(E65.2.bin140)[2]), rep("E75.1", dim(E75.1.bin140)[2]), rep("E75.2", dim(E75.2.bin140)[2]), 
                                             rep("E85.B6", dim(B6.bin140)[2]), rep("E85.C1", dim(C1.bin140)[2]), rep("E95", dim(E95.bin140)[2])))
UMI.bin140.df$bin <- "bin140"
UMI.bin140.df$sample <- factor(UMI.bin140.df$sample, levels = c("E65.2","E65.1","E75.1","E75.2","E85.B6","E85.C1","E95"))

UMI.df <- rbind(UMI.bin3.df, UMI.bin14.df, UMI.bin140.df)
save(UMI.df, file = "20220622-Comments/06.ST_QC/All-ST-bin3-14-140-UMI-count.RData")

pdf("20220622-Comments/06.ST_QC/03.ST-bin3-UMI-boxplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_boxplot(UMI.bin3.df, mapping = aes(x = sample, y = UMI, fill = sample), outlier.shape = NA) + 
  ylim(0, 150) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

pdf("20220622-Comments/06.ST_QC/03.ST-bin14-UMI-boxplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_boxplot(UMI.bin14.df, mapping = aes(x = sample, y = UMI, fill = sample), outlier.shape = NA) + 
  ylim(0, 2000) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

pdf("20220622-Comments/06.ST_QC/03.ST-bin140-UMI-boxplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_boxplot(UMI.bin140.df, mapping = aes(x = sample, y = UMI, fill = sample), outlier.shape = NA) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

load("20220622-Comments/05.DBA/DBA-E3-data.Rdata")
load("20220622-Comments/05.DBA/DBA-E6-data.Rdata")
load("20220622-Comments/05.DBA/DBA-F1-data.Rdata")
load("20220622-Comments/05.DBA/DBA-E3-data-bin20.Rdata")
load("20220622-Comments/05.DBA/DBA-E6-data-bin20.Rdata")
load("20220622-Comments/05.DBA/DBA-F1-data-bin20.Rdata")

DBA.E3.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-celltype.csv", sep = ",", header = T, check.names = F)
DBA.E6.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-celltype.csv", sep = ",", header = T, check.names = F)
DBA.F1.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin50-celltype.csv", sep = ",", header = T, check.names = F)
View(sort(DBA.E3.bin50.celltype %>% group_by(cell_type) %>% summarize(median_count = median(nCount_RNA), median_feature = median(nFeature_RNA))))