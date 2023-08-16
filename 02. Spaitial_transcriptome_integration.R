#Author: Min Yang

#1.  scRNA-seq data SCT conversion 
load("seruat.harmony.object.RData")
seruat.harmony.object.SCT <- seruat.harmony.object
seruat.harmony.object.SCT <- SCTransform(seruat.harmony.object.SCT, ncells = 3000, verbose = FALSE)
seruat.harmony.object.SCT <- RunPCA(seruat.harmony.object.SCT)
seruat.harmony.object.SCT <- RunUMAP(seruat.harmony.object.SCT, dims = 1:30)
save(seruat.harmony.object.SCT, file = "seruat.harmony.object.SCT.RData")

ST.data <- read.table("ST.bin1.Lasso.gem", sep = "\t", header = T, check.names = F)# read data matrix
#ST.data includes 4 columns: c1_geneID IDSC, c2_x, c3_y, c4_MIDCounts
Creat.ST.object <- function(ST.data, bin = 50){
  ST.data.bin50 <- ST.data
  ST.data.bin50$x <- trunc(ST.data.bin50$x/bin) * bin
  ST.data.bin50$y <- trunc(ST.data.bin50$y/bin) * bin
  ST.data.bin50$cellID <- paste(ST.data.bin50$x, "_", ST.data.bin50$y, sep = "")
  ST.data.bin50 <- aggregate(ST.data.bin50$MIDCounts, by = list(ST.data.bin50$cellID, ST.data.bin50$geneID), sum)
  colnames(ST.data.bin50) <- c("cellID", "geneID", "MIDCounts")
  ST.data.bin50$cellInx <- match(ST.data.bin50$cellID, unique(ST.data.bin50$cellID))
  ST.data.bin50$cellInx <- match(ST.data.bin50$cellID, unique(ST.data.bin50$cellID))
  ST.data.bin50$geneInx <- match(ST.data.bin50$geneID, unique(ST.data.bin50$geneID))
  mat <- sparseMatrix(i = ST.data.bin50$geneInx, j = ST.data.bin50$cellInx, x = ST.data.bin50$MIDCounts, 
                      dimnames = list(unique(ST.data.bin50$geneID), unique(ST.data.bin50$cellID)))
  ST.bin50.coord.df <- data.frame(cellname = colnames(mat))
  rownames(ST.bin50.coord.df) <- ST.bin50.coord.df$cellname
  ST.bin50.coord.df <- separate(ST.bin50.coord.df, col = cellname, sep = "_", into = c("x", "y"))
  ST.bin50 <- CreateSeuratObject(mat, project = "ST.bin50", assay = "Spatial")
  ST.bin50$slice <- 1
  ST.bin50$region <- "ST.bin50"
  colnames(ST.bin50.coord.df) <- c("imagerow", "imagecol")
  ST.bin50.coord.df$imagerow <- as.numeric(ST.bin50.coord.df$imagerow)
  ST.bin50.coord.df$imagecol <- as.numeric(ST.bin50.coord.df$imagecol)
  ST.bin50@images$ST_bin50 <- new(Class = "SlideSeq", assay = "spatial", key = "image_", coordinates = ST.bin50.coord.df)
  #SCT transform
  ST.bin50 <- SCTransform(ST.bin50, assay = "Spatial", verbose = FALSE)
  ST.bin50 <- RunPCA(ST.bin50, assay = "SCT", verbose = F)
  ElbowPlot(ST.bin50, ndims = 50)
  ST.bin50 <- FindNeighbors(ST.bin50, reduction = "pca", dims = 1:20)
  ST.bin50 <- FindClusters(ST.bin50, verbose = T, resolution = 0.3)
  ST.bin50 <- RunUMAP(ST.bin50, reduction = "pca", dims = 1:20)
  return(ST.bin50)
}

#2. Creat the seurat object for each data on each slide
ST.seurat.object <- Creat.ST.object(ST.data, bin)#bin %in% c(10, 20, 50, 100), ST.data %in% c(E65_1, E65_2, E75_1, E75_2, E85_1, E85_2, E95, CBA_E85, CBA_E95a, CBA_E95b)

#integrate scRNAs-seq and ST.seurat.object
options(future.globals.maxSize = 2 * 1024^3)#2G
plan("multiprocess", workers = 10)

anchors <- FindTransferAnchors(reference = seruat.harmony.object.SCT, query = ST.seurat.object, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = seruat.harmony.object.SCT$subname, 
                                  prediction.assay = TRUE, weight.reduction = ST.seurat.object[["pca"]], dims = 1:30)
ST.seurat.object[["predictions"]] <- predictions.assay
DefaultAssay(ST.seurat.object) <- "predictions"
save(ST.seurat.object, file = "ST.seurat.object_map_to_scRNA-seq.RData")

#3. Assign the cell type for each spatial spot
DefaultAssay(ST.seurat.object) <-  "predictions"
ST.seurat.object.celltype.score <- GetAssayData(ST.seurat.object, slot = "data")
ST.seurat.object.celltype.max <- apply(ST.seurat.object.celltype.score, 2, function(x) rownames(ST.seurat.object.celltype.score)[which.max(x)])
ST.seurat.object.celltype <- data.frame(cell_type = ST.seurat.object.celltype.max, 
                                        cell_name = names(ST.seurat.object.celltype.max), 
                                        GetTissueCoordinates(ST.seurat.object), 
                                        spatial_clusters = ST.seurat.object$seurat_clusters,
                                        nFeature_RNA = ST.seurat.object$nFeature_Spatial, 
                                        nCount_RNA = ST.seurat.object$nCount_Spatial)
write.table(ST.seurat.object.celltype, file = "ST.seurat.object.celltype.csv", sep = ",", quote = F)

#3.1 For E9.5 slide, we run integration twice to exactly annotate the cells from maternal and fetal.
#3.1.1 After first round integration, we extract the embryo regions as follows:
E95.bin50.embryo.cell <- c("Mesenchymal-epithelial transition", "Blood P", "Erythroid", "Mesenchymal stem cell", "Allantois mesodermal", "Fetal EC")
E95.embryo <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% E95.bin50.embryo.cell, ]
E95.embryo.coord.1 <- E95.embryo[E95.embryo$x >= 18000 & E95.embryo$x <= 21500 & E95.embryo$y >= 14500 & E95.embryo$y <= 24500, ]
E95.embryo.coord.2 <- E95.embryo[E95.embryo$x >= 21500 & E95.embryo$x <= 23200 & E95.embryo$y >= 15000 & E95.embryo$y <= 21000, ]
E95.embryo.coord <- rbind(E95.embryo.coord.1, E95.embryo.coord.2)
E95.embryo.coord$cell_type <- "embryo" #label all the cells in this region as embryo
#3.1.2 Exclude the spatial spots labeled as embryo and create a subset dataset (maternal) and integrate to the maternal-derived scRNA-seq again.
E95.bin50.seurat.maternal.object <- subset(E95.bin50.seurat.object, cells = setdiff(colnames(E95.bin50.seurat.object), rownames(E95.embryo.coord)))
subname <- unique(UE.harmony.SCT$subname)
seruat.harmony.object.SCT.maternal <- subset(seruat.harmony.object.SCT, subset = subname %in% setdiff(subname, E95.bin50.embryo.cell))
anchors <- FindTransferAnchors(reference = seruat.harmony.object.SCT.maternal, query = E95.bin50.seurat.maternal.object, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = seruat.harmony.object.SCT.maternal$subname, 
                                  prediction.assay = TRUE, weight.reduction = E95.bin50.seurat.maternal.object[["pca"]], dims = 1:30)
E95.bin50.seurat.maternal.object[["predictions"]] <- predictions.assay
DefaultAssay(E95.bin50.seurat.maternal.object) <- "predictions"

E95.bin50.seurat.object.celltype.maternal.score <- GetAssayData(E95.bin50.seurat.maternal.object, slot = "data")
E95.bin50.seurat.object.celltype.maternal.max <- apply(E95.bin50.seurat.object.celltype.score, 2, function(x) rownames(E95.bin50.seurat.object.celltype.score)[which.max(x)])
E95.bin50.seurat.object.celltype.maternal <- data.frame(cell_type = E95.bin50.seurat.object.celltype.maternal.max, 
                                                        cell_name = names(E95.bin50.seurat.object.celltype.maternal.max), 
                                                        GetTissueCoordinates(E95.bin50.seurat.maternal.object), 
                                                        spatial_clusters = E95.bin50.seurat.maternal.object$seurat_clusters,
                                                        nFeature_RNA = E95.bin50.seurat.maternal.object$nFeature_Spatial, 
                                                        nCount_RNA = E95.bin50.seurat.maternal.object$nCount_Spatial)
E95.bin50.celltype.df <- rbind(E95.bin50.seurat.object.celltype.maternal, E95.embryo.coord) #This is the final integration results for E9.5

#4. CBA samples follow the same as the normal samples. E95b sample follows the same strategy as normal E9.5 sample



