#Author: Min Yang
#13 Major cell type identification and subclusters of DSCs/Immune cells/ECs/Trophoblast cells/iDSCs

#import all the libraries
library(ggplot2)
library(Seurat)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(patchwork)
library(stringr)

#We have 10 normal and 1 CBA/J samples in total, as listed in Table S1 and Figure S1B

###1. Create Seruat object and filter cells/genes follow the rules in Table S1 
object.data <- Read10X("filtered_feature_bc_matrix")
seurat.object <- CreateSeuratObject(object.data, min.cells = 3, min.features = 200)
seurat.object$percent.mt <- PercentageFeatureSet(seurat.object, pattern = "^mt-")
seurat.object <- subset(seurat.object, subset = percent.mt < 5 & 
                          nFeature_RNA > 500 & 
                          nFeature_RNA<quantile(U65$nFeature_RNA, 0.99))
save(seurat.object, file = "seurat.object-filter.Rdata")

###2. Merge all the normal seurat objects and integrate according to the sample time(T5.5, T6.5, T7.5, T8.5, T9.5, T10.5)
#Add the time metadata for each object
object.list <- list(object.seurat1, object.seurat2, ...)
names(object.list) <- c("object.seurat1", "object.seurat2", ...)
merged.seurat.object <- merge(object.list[[1]], y = object.list[2:10], 
                              add.cell.ids = names(object.list), 
                              project = "UE.merge.filter")
#Using harmony to integrate
seruat.harmony.object <- RunHarmony(merged.seurat.object, "time", plot_convergence = TRUE)
DimPlot(object = seruat.harmony.object, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
seruat.harmony.object <- RunUMAP(seruat.harmony.object, reduction = "harmony", dims = 1:30)
save(seruat.harmony.object, file = "seruat.harmony.object-pc30-UMAP.RData")
seruat.harmony.object <- FindNeighbors(object = seruat.harmony.object, dims = 1:30, reduction = "harmony")
seruat.harmony.object <- FindClusters(seruat.harmony.object, resolution = 0.8)
DimPlot(seruat.harmony.object, pt.size = 0.15, label = T)
save(seruat.harmony.object, file = "seruat.harmony.object.RData")

###3. Rename the cell clusters into 13 major cell types under resolution=0.8

###4. subset DSCs/Immune cells/ECs/Trophoblast cells to reclustering
#4.1 subclustering DSCs
DSC.seurat.object <- susbet(seruat.harmony.object, subset = cellname %in% c("DSC"))
DSC.seurat.object <- NormalizeData(DSC.seurat.object) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = VariableFeatures(DSC.seurat.object))
DSC.seurat.object <- RunPCA(DSC.seurat.object, features = VariableFeatures(DSC.seurat.object))
DSC.seurat.object <- RunHarmony(DSC.seurat.object, "time", plot_convergence = TRUE, reduction = "pca")
DimPlot(object = DSC.seurat.object, reduction = "harmony", pt.size = .1, group.by = "time")
#renaming the DSC subclusters as: D1_eSF, D2_pre-DSC, D3_proliferating DSC, D4_Biomineral-regulatory DSC, D5_Nourishing DSC, D6_Angiogenic DSC, D7_Postmature DSC
#removing the low quality and contaminated subclusters, which both expressed trophoblast and DSC featured genes.

#4.2 subclustering immune cells
Immune.seurat.object <- susbet(seruat.harmony.object, subset = cellname %in% c("Immune cell"))
Immune.seurat.object <- NormalizeData(Immune.seurat.object) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = VariableFeatures(Immune.seurat.object))
Immune.seurat.object <- RunPCA(Immune.seurat.object, features = VariableFeatures(Immune.seurat.object))
Immune.seurat.object <- RunHarmony(Immune.seurat.object, "time", plot_convergence = TRUE, reduction = "pca")
DimPlot(object = Immune.seurat.object, reduction = "harmony", pt.size = .1, group.by = "time")
#renaming the Immune cell clusters as: Monocyte, Macrophage, NK, Proliferating Mac, DC-1, DC-2, iDSC, T cell, B cell, Neutrophil

#4.3 subclustering ECs
EC.seurat.object <- susbet(seruat.harmony.object, subset = cellname %in% c("EC"))
EC.seurat.object <- NormalizeData(EC.seurat.object) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(features = VariableFeatures(EC.seurat.object))
EC.seurat.object <- RunPCA(EC.seurat.object, features = VariableFeatures(EC.seurat.object))
EC.seurat.object <- RunHarmony(EC.seurat.object, "time", plot_convergence = TRUE, reduction = "pca")
DimPlot(object = EC.seurat.object, reduction = "harmony", pt.size = .1, group.by = "time")
#renaming the EC cell clusters as: Angiogenic EC, Proliferating EC, Venous EC1, Venous EC2, Arterial EC, Lymphatic EC

#4.4 subclustering Trophoblast cells
Trophoblast.seurat.object <- susbet(seruat.harmony.object, subset = cellname %in% c("Trophoblast"))
Trophoblast.seurat.object <- FindVariableFeatures(Trophoblast.seurat.object, selection.method = "vst", nfeatures = 2000)
Trophoblast.seurat.object <- ScaleData(Trophoblast.seurat.object, features = VariableFeatures(Trophoblast.seurat.object))
Trophoblast.seurat.object <- RunPCA(Trophoblast.seurat.object, features = VariableFeatures(Trophoblast.seurat.object))
Trophoblast.seurat.object <- RunHarmony(Trophoblast.seurat.object, "time", plot_convergence = TRUE, reduction = "pca")
DimPlot(object = Trophoblast.seurat.object, reduction = "harmony", pt.size = .1, group.by = "time")
#renaming the EC cell clusters as: TBPC1, TBPC2, EPC1, EPC2, TGC, SpT

#4.5 subclustering iDSCs
iDSC <- subset(Immune.seurat.object, subset = cellname %in% c("iDSC"))
iDSC <- SCTransform(iDSC, ncells = 3000, verbose = FALSE, return.only.var.genes = FALSE)
iDSC <- RunPCA(iDSC)
iDSC <- RunUMAP(iDSC, dims = 1:30)
iDSC <- FindNeighbors(iDSC, dims = 1:20)
iDSC <- FindClusters(iDSC, resolution = 0.1)
save(iDSC, file = "iDSC.seurat.object.RData")

#4.6 After reclustering, then, saving the names of subclusters in subname metadta

#CBA/J sample clustering
CBA.data <- Read10X("filtered_feature_bc_matrix/")
CBA <- CreateSeuratObject(CBA.data, min.cells = 3, min.features = 200, project = "CBA")
CBA$percent.mt <- PercentageFeatureSet(CBA, pattern = "^mt-")
CBA <- subset(CBA, subset = percent.mt < 10 & 
                  nFeature_RNA > 500 & 
                  nFeature_RNA<quantile(CBA$nFeature_RNA, 0.99))
CBA <- NormalizeData(CBA, normalization.method = "LogNormalize", scale.factor = 10000)
CBA <- FindVariableFeatures(CBA, selection.method = "vst", nfeatures = 2000)
CBA <- ScaleData(CBA, features = VariableFeatures(CBA))
CBA <- RunPCA(CBA, features = VariableFeatures(CBA), npcs = 70)
ElbowPlot(CBA, ndims = 70)
CBA <- FindNeighbors(object = CBA, dims = 1:30)
CBA <- RunUMAP(CBA, dims = 1:30, seed.use = 0912)
CBA <- FindClusters(CBA, resolution = 0.3)
DimPlot(CBA, label = T)
save(CBA, file = "CBA.seurat.object.Rdata")
#Rename the cell clusters















