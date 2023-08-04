library("FactoMineR")
library("factoextra")
library("pheatmap")
library("limma")
library("Seurat")
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library('harmony')
library("reshape2")
DBA.1 <- Read10X("../../20220909-DBA-SC/CBA0825/outs/filtered_feature_bc_matrix/")
DBA.1 <- CreateSeuratObject(DBA.1, min.cells = 3, min.features = 200, project = "DBA.1")
DBA.1$percent.mt <- PercentageFeatureSet(DBA.1, pattern = "^mt-")
VlnPlot(DBA.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
DBA.1 <- subset(DBA.1, subset = percent.mt < 10 & 
                  nFeature_RNA > 500 & 
                  nFeature_RNA<quantile(DBA.1$nFeature_RNA, 0.99))
DBA.1 <- NormalizeData(DBA.1, normalization.method = "LogNormalize", scale.factor = 10000)
DBA.1 <- FindVariableFeatures(DBA.1, selection.method = "vst", nfeatures = 2000)
DBA.1 <- ScaleData(DBA.1, features = VariableFeatures(DBA.1))
DBA.1 <- RunPCA(DBA.1, features = VariableFeatures(DBA.1), npcs = 70)
ElbowPlot(DBA.1, ndims = 70)
DBA.1 <- FindNeighbors(object = DBA.1, dims = 1:30)
DBA.1 <- RunUMAP(DBA.1, dims = 1:30, seed.use = 0912)
DBA.1 <- FindClusters(DBA.1, resolution = 0.3)
DimPlot(DBA.1, label = T)
FeaturePlot(DBA.1, features = c("Krt8", "H19"))
DBA.1.markers <- FindAllMarkers(DBA.1, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.1.markers.top5 <- DBA.1.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(DBA.1, features = DBA.1.markers.top5$gene)
save(DBA.1, file = "20210827-figures/files/DBA.1-20220912.Rdata")
pdf("20220622-Comments/05.DBA/00.figures/07.DBA-major-genes-featureplot.pdf", width = 15, height = 15)
FeaturePlot(DBA.1, features = c("Ptprc", "Gzma", "Hand2", "Pgr", "Kdr", "Rgs5", "Epcam"), pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()
DBA.1 <- RenameIdents(DBA.1, "1" = "Decidual", "4" = "Decidual", "5" = "Decidual", "6" = "Decidual", "8" = "Decidual", "13" = "Decidual", 
                      "0" = "Immune", "2" = "Immune", "12" = "Immune", "11" = "Immune", "9" = "Immune", "3" = "EC", "7" = "SMC", "10" = "Epithelial")
pdf("20220622-Comments/05.DBA/00.figures/08.DBA1-UMAP-color-major-cluster.pdf", width = 10, height = 10)
DimPlot(DBA.1, pt.size = 1, cols = major.color[c(1:5)])
dev.off()

#subset DSCs to integrate with normal samples
DBA.1.DSC <- subset(DBA.1, subset = seurat_clusters %in% c(1,4,6,8,13))
DBA.1.DSC$time <- "DBA.1"
DBA.1.DSC$cellname <- "DBA.decidua"
load("FB-immune-DSC-harmony-20211214.RData")
UE.DBA.decidual <- merge(UE.FB, DBA.1.DSC)
UE.DBA.decidual <- NormalizeData(UE.DBA.decidual, normalization.method = "LogNormalize", scale.factor = 10000)
UE.DBA.decidual <- FindVariableFeatures(UE.DBA.decidual, selection.method = "vst", nfeatures = 2000)
UE.DBA.decidual <- ScaleData(UE.DBA.decidual, features = VariableFeatures(UE.DBA.decidual))
UE.DBA.decidual <- RunPCA(UE.DBA.decidual, features = VariableFeatures(UE.DBA.decidual))
UE.DBA.decidual <- RunHarmony(UE.DBA.decidual, "time", plot_convergence = TRUE, reduction = "pca")
DimPlot(object = UE.DBA.decidual, reduction = "harmony", pt.size = .1, group.by = "time")
ElbowPlot(UE.DBA.decidual, ndims = 80, reduction = "harmony")
UE.DBA.decidual <- RunUMAP(UE.DBA.decidual, reduction = "harmony", dims = 1:30)
UE.DBA.decidual <- FindNeighbors(UE.DBA.decidual, reduction = "harmony", dims = 1:30)
UE.DBA.decidual <- FindClusters(UE.DBA.decidual, resolution = 0.2)
UE.DBA.decidual <- RenameIdents(UE.DBA.decidual, "6" = "D1_Lum", "4" = "D2_Sfrp4", "5" = "D3_Top2a", "2" = "D4_Ptn", "0" = "D5_Gatm", "1" = "D6_S100a8", "3" = "D3_Top2a", 
                                "8" = "D7_Prl8a2", "7" = "iDSC_Ptprc", "9" = "iDSC_Ptprc")
save(UE.DBA.decidual, file = "20210827-figures/files/DBA.1-DSC-20220912.Rdata")
pdf("20220622-Comments/05.DBA/00.figures/02.DBA1-deciuda-integrate-UMAP-color-cluster.pdf", width = 10, height = 10)
DimPlot(UE.DBA.decidual, cols = c(decidual.color[1:7], immune.color[10]), pt.size = 0.5)
dev.off()
DimPlot(UE.DBA.decidual, group.by = "cellname", label = T)
pdf("20220622-Comments/05.DBA/00.figures/01.DBA1-deciuda-integrate-UMAP-color-sample.pdf", width = 10, height = 10)
DimPlot(UE.DBA.decidual, group.by = "time", cols = c("#00BFFF", rep("#FA8072", 6)), pt.size = 0.5)
dev.off()
UE.DBA.decidual.markers <- FindAllMarkers(UE.DBA.decidual, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
UE.DBA.decidual.markers.top5 <- UE.DBA.decidual.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(UE.DBA.decidual, features = UE.DBA.decidual.markers.top5$gene)
#subset DSCs to calculate percentage
UE.DBA.decidual$cellcluster <- UE.DBA.decidual@active.ident
UE.DBA.DSC <- subset(UE.DBA.decidual, subset = cellcluster != 'iDSC_Ptprc')
UE.DBA.DSC.proportion <- data.frame(table(UE.DBA.DSC$time, UE.DBA.DSC$cellcluster), check.names = F)
UE.DBA.DSC.proportion$Var2 <- factor(UE.DBA.DSC.proportion$Var2, 
                                     levels = c("D1_Lum", "D2_Sfrp4", "D3_Top2a", "D4_Ptn", "D5_Gatm", "D6_S100a8", "D7_Prl8a2"))
UE.DBA.DSC.proportion$Var1 <- factor(UE.DBA.DSC.proportion$Var1, levels = c("T55", "T65", "T75", "T85", "T95", "T105", "DBA.1"))
decidual.color <- c("#6A5ACD","#9ACD32","#FB8072","#8DD3C7","#08519C","#FDB462","#228B22","#FCCDE5","#D9D9D9")
pdf("20220622-Comments/05.DBA/00.figures/03.DBA-decidua-cluster-proportion.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(UE.DBA.DSC.proportion, mapping = aes(x = Var1, y = Freq, fill = Var2), stat = "identity", position = "fill") + 
  scale_fill_manual(values = decidual.color) + 
  theme_classic()
dev.off()
#subset iDSCs to calculate percentage
UE.DBA.iDSC <- subset(UE.DBA.decidual, subset = cellcluster == 'iDSC_Ptprc')
UE.DBA.iDSC <- SCTransform(UE.DBA.iDSC, ncells = 3000, verbose = FALSE, return.only.var.genes = FALSE)
UE.DBA.iDSC <- RunPCA(UE.DBA.iDSC)
UE.DBA.iDSC <- RunUMAP(UE.DBA.iDSC, dims = 1:30)
UE.DBA.iDSC <- FindNeighbors(UE.DBA.iDSC, dims = 1:20)
UE.DBA.iDSC <- FindClusters(UE.DBA.iDSC, resolution = 0.2)
DoHeatmap(UE.DBA.iDSC, features = UE.DBA.iDSC.markers.top10$gene) + NoLegend()
UE.DBA.iDSC <- RenameIdents(UE.DBA.iDSC, "0" = "iDSC_0", "2" = "iDSC_1", "3" = "iDSC_2", "1" = "iDSC_DBA")
UE.DBA.iDSC$cellcluster <- UE.DBA.iDSC@active.ident
save(UE.DBA.iDSC, file = "20210827-figures/files/DBA.1-iDSC-20220912.Rdata")
pdf("20220622-Comments/05.DBA/00.figures/05.DBA-deciuda-iDSC-UMAP-color-clusters.pdf", width = 10, height = 10)
DimPlot(UE.DBA.iDSC, pt.size = 1, cols = c("#ED736E", "#BD7FF8", "#55BB3C", "#5B9DF8"))
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/04.DBA-deciuda-iDSC-UMAP-color-sample.pdf", width = 10, height = 10)
DimPlot(UE.DBA.iDSC, pt.size = 1, group.by = "time", cols = c("#00BFFF", rep("#FA8072", 6)))
dev.off()
UE.DBA.iDSC.markers <- FindAllMarkers(UE.DBA.iDSC, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
UE.DBA.iDSC.markers.top20 <- UE.DBA.iDSC.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
pdf("20220622-Comments/05.DBA/00.figures/05.DBA-deciuda-iDSC-DEGs-heatmap.pdf", width = 10, height = 10)
DoHeatmap(UE.DBA.iDSC, features = UE.DBA.iDSC.markers.top20$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()
UE.DBA.iDSC.proportion <- data.frame(table(UE.DBA.iDSC$time, UE.DBA.iDSC$cellcluster), check.names = F)

library("clusterProfiler")
enrich.GO.function <- function(geneset){
  enrich.GO <- enrichGO(gene = geneset,
                        OrgDb = "org.Mm.eg.db",
                        keyType = "SYMBOL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)
  enrich.GO <- simplify(enrich.GO, cutoff = 0.7, by = "p.adjust", select_fun = min)
  enrich.GO <- enrich.GO@result[enrich.GO@result$pvalue<0.01, ]
  return(enrich.GO)
}
iDSC.GO.df <- data.frame()
for (c in unique(UE.DBA.iDSC.markers$cluster)) {
  print(c)
  geneset <- UE.DBA.iDSC.markers[UE.DBA.iDSC.markers$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  iDSC.GO.df <- rbind(iDSC.GO.df, tmp)
}
iDSC.GO.df$GO_Term <- factor(iDSC.GO.df$GO_Term, levels = rev(unique(iDSC.GO.df$GO_Term)))
iDSC.GO.df$cluster <- factor(iDSC.GO.df$cluster, levels = c("iDSC_0", "iDSC_1", "iDSC_2", "iDSC_DBA"))
iDSC.GO.df.top5 <- iDSC.GO.df %>% group_by(cluster) %>% top_n(n = 5, -log(pvalue))
pdf("20220622-Comments/05.DBA/00.figures/06.DBA-deciuda-iDSC-DEGs-GO-terms-barplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(iDSC.GO.df.top5, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()

#EC integration
DBA.1.EC <- subset(DBA.1, subset = seurat_clusters %in% c(3))
DBA.1.EC$time <- "DBA.1"
DBA.1.EC$cellname <- "DBA.EC"
load("../20210201-integrate/UE-EC-harmony-pc30-rename.RData")
UE.EC <- RenameIdents(UE.EC, "EC-2" = "Cycling EC", "EC-0" = "Angiogeneic EC", "EC-4" = "Arterial EC", "EC-1" = "Venous EC1", "EC-3" = "Venous EC2", "EC-5" = "Lympatic EC")
DimPlot(UE.EC)
UE.EC$cellname <- UE.EC@active.ident
UE.DBA1.EC <- merge(UE.EC, DBA.1.EC)
UE.DBA1.EC <- FindVariableFeatures(UE.DBA1.EC, selection.method = "vst", nfeatures = 2000)
UE.DBA1.EC$CC.Score <- UE.DBA1.EC$S.Score - UE.DBA1.EC$G2M.Score
UE.DBA1.EC <- ScaleData(UE.DBA1.EC, features = VariableFeatures(UE.DBA1.EC))
UE.DBA1.EC <- RunPCA(UE.DBA1.EC, features = VariableFeatures(UE.DBA1.EC))
UE.DBA1.EC <- RunHarmony(UE.DBA1.EC, "time", plot_convergence = TRUE, reduction = "pca")
DimPlot(object = UE.DBA1.EC, reduction = "harmony", pt.size = .1, group.by = "time")
ElbowPlot(UE.DBA1.EC, ndims = 80, reduction = "harmony")
UE.DBA1.EC <- RunUMAP(UE.DBA1.EC, reduction = "harmony", dims = 1:30)
UE.DBA1.EC <- FindNeighbors(UE.DBA1.EC, reduction = "harmony", dims = 1:30)
UE.DBA1.EC <- FindClusters(UE.DBA1.EC, resolution = 0.1)
pdf("20220622-Comments/05.DBA/00.figures/09.DBA1-EC-integrate-UMAP-color-sample.pdf", width = 10, height = 10)
DimPlot(UE.DBA1.EC, group.by = "time", cols = c("#00BFFF", rep("#FA8072", 6)), pt.size = 1)
dev.off()

DimPlot(UE.DBA1.EC, group.by = "cellname", label = T) 
UE.DBA1.EC <- RenameIdents(UE.DBA1.EC, "0" = "EC0_Angiogenic", "1" = "EC1_Venous EC1", "2" = "EC2_Cycling", "3" = "EC3_Venouse EC2", 
                           "4" = "EC4_Artery", "5" = "EC3_Venouse EC2", "6" = "EC5_Lymphotic")
UE.DBA1.EC$sample <- UE.DBA1.EC$time
UE.DBA1.EC$sample <- gsub("T.*", "Normal", UE.DBA1.EC$sample)
UE.DBA1.EC$cellcluster <- UE.DBA1.EC@active.ident
save(UE.DBA1.EC, file = "20220622-Comments/05.DBA/UE-DBA1-EC-integrate.RData")

pdf("20220622-Comments/05.DBA/00.figures/09.DBA1-EC-integrate-UMAP-color-cluster.pdf", width = 10, height = 10)
DimPlot(UE.DBA1.EC, cols = EC.color, pt.size = 1)
dev.off()
#proportion
UE.DBA1.EC.df <- data.frame(table(UE.DBA1.EC$time, UE.DBA1.EC@active.ident), check.names = F)
UE.DBA1.EC.df$Var1 <- factor(UE.DBA1.EC.df$Var1, levels = c("T55", "T65", "T75", "T85", "T95", "T105", "DBA.1"))
pdf("20220622-Comments/05.DBA/00.figures/10.DBA1-EC-integrate-proportion-barplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(UE.DBA1.EC.df, mapping = aes(x = Var1, y = Freq, fill = Var2), stat = "identity", position = "fill") + 
  scale_fill_manual(values = EC.color) + 
  theme_classic()
dev.off()

#############################20221101 DBA2 loss pregancy#######################
DBA.2 <- Read10X("../../20220909-DBA-SC/CBA0901/outs/filtered_feature_bc_matrix/")
DBA.2 <- CreateSeuratObject(DBA.2, min.cells = 3, min.features = 200, project = "DBA.2")
DBA.2$percent.mt <- PercentageFeatureSet(DBA.2, pattern = "^mt-")
VlnPlot(DBA.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
DBA.2 <- subset(DBA.2, subset = percent.mt < 10 & 
                  nFeature_RNA > 500 & 
                  nCount_RNA > 2000 &
                  nFeature_RNA<quantile(DBA.2$nFeature_RNA, 0.99))
DBA.2 <- NormalizeData(DBA.2, normalization.method = "LogNormalize", scale.factor = 10000)
DBA.2 <- FindVariableFeatures(DBA.2, selection.method = "vst", nfeatures = 2000)
DBA.2 <- ScaleData(DBA.2, features = rownames(DBA.2))
DBA.2 <- RunPCA(DBA.2, features = VariableFeatures(DBA.2), npcs = 70)
ElbowPlot(DBA.2, ndims = 70)
DBA.2 <- FindNeighbors(object = DBA.2, dims = 1:30)
DBA.2 <- FindClusters(DBA.2, resolution = 0.3)
DBA.2 <- RunUMAP(DBA.2, dims = 1:30, seed.use = 0912)
DimPlot(DBA.2, label = T)
DBA.2.markers <- FindAllMarkers(DBA.2, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.2.markers.top5 <- DBA.2.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(DBA.2, features = DBA.2.markers.top5$gene)
save(DBA.2, file = "20210827-figures/files/DBA.2-20220912.Rdata")
FeaturePlot(DBA.2, features = c("Ptprc", "Gzma", "Hand2", "Pgr", "Kdr", "Rgs5", "Epcam"), pt.size = 0.5, ncol = 3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

####comparing the DSC/iDSC clusters between normal and DBA model
UE.DBA.decidual$sample <- UE.DBA.decidual$orig.ident
UE.DBA.decidual$sample[grep("[U,E].*", UE.DBA.decidual$sample)] <- "Normal"
UE.DBA.decidual.data <- AverageExpression(UE.DBA.decidual, add.ident = c("sample"), assays = "RNA")
UE.DBA.decidual.data <- UE.DBA.decidual.data$RNA
UE.DBA.decidual.data.cor <- cor(UE.DBA.decidual.data)
UE.DBA.decidual.data.cor <- UE.DBA.decidual.data.cor[c(16,15,11,9,14,10,13,12), c(5,3,1,4,2,7,6,8)]
pheatmap(UE.DBA.decidual.data.cor, scale = "row", cluster_rows = F, cluster_cols = F)

UE.DBA.decidual.data.nor <- data.frame(t(apply(UE.DBA.decidual.data, 1, function(x){x/colSums(UE.DBA.decidual.data) * 1000000})), check.names = F)
UE.DBA.decidual.data.nor.D1 <- (UE.DBA.decidual.data.nor$D1_Lum_Normal + 1)/(UE.DBA.decidual.data.nor$D1_Lum_DBA.1 + 1)
UE.DBA.decidual.data.nor.D2 <- (UE.DBA.decidual.data.nor$D2_Sfrp4_Normal + 1)/(UE.DBA.decidual.data.nor$D2_Sfrp4_DBA.1 + 1)
UE.DBA.decidual.data.nor.D3 <- (UE.DBA.decidual.data.nor$D3_Top2a_Normal + 1)/(UE.DBA.decidual.data.nor$D3_Top2a_DBA.1 + 1)
UE.DBA.decidual.data.nor.D4 <- (UE.DBA.decidual.data.nor$D4_Ptn_Normal + 1)/(UE.DBA.decidual.data.nor$D4_Ptn_DBA.1 + 1)
UE.DBA.decidual.data.nor.D5 <- (UE.DBA.decidual.data.nor$D5_Gatm_Normal + 1)/(UE.DBA.decidual.data.nor$D5_Gatm_DBA.1 + 1)
UE.DBA.decidual.data.nor.D6 <- (UE.DBA.decidual.data.nor$D6_S100a8_Normal + 1)/(UE.DBA.decidual.data.nor$D6_S100a8_DBA.1 + 1)
UE.DBA.decidual.data.nor.D7 <- (UE.DBA.decidual.data.nor$D7_Prl8a2_Normal + 1)/(UE.DBA.decidual.data.nor$D7_Prl8a2_DBA.1 + 1)
UE.DBA.decidual.data.nor.iDSC <- (UE.DBA.decidual.data.nor$iDSC_Ptprc_Normal + 1)/(UE.DBA.decidual.data.nor$iDSC_Ptprc_DBA.1 + 1)
ggplot() + 
  geom_point(UE.DBA.decidual.data.nor, mapping = aes(x = log2(D1_Lum_Normal+1), y = log2(D1_Lum_DBA.1+1)), size = 0.1)
ggplot() + 
  geom_point(UE.DBA.decidual.data.nor, mapping = aes(x = log2(D2_Sfrp4_Normal+1), y = log2(D2_Sfrp4_DBA.1+1)), size = 0.1)
ggplot() + 
  geom_point(UE.DBA.decidual.data.nor, mapping = aes(x = log2(D3_Top2a_Normal+1), y = log2(D3_Top2a_DBA.1+1)), size = 0.1) + 
  ylim(0, 15) + 
  xlim(0, 15)
ggplot() + 
  geom_point(UE.DBA.decidual.data.nor, mapping = aes(x = log2(D4_Ptn_Normal+1), y = log2(D4_Ptn_DBA.1+1)), size = 0.1)+ 
  ylim(0, 15) + 
  xlim(0, 15)
ggplot() + 
  geom_point(UE.DBA.decidual.data.nor, mapping = aes(x = log2(D5_Gatm_Normal+1), y = log2(D5_Gatm_DBA.1+1)), size = 0.1)+ 
  ylim(0, 15) + 
  xlim(0, 15)
ggplot() + 
  geom_point(UE.DBA.decidual.data.nor, mapping = aes(x = log2(D6_S100a8_Normal+1), y = log2(D6_S100a8_DBA.1+1)), size = 0.1)+ 
  ylim(0, 15) + 
  xlim(0, 15)
ggplot() + 
  geom_point(UE.DBA.decidual.data.nor, mapping = aes(x = log2(D7_Prl8a2_Normal+1), y = log2(D7_Prl8a2_DBA.1+1)), size = 0.1)+ 
  ylim(0, 15) + 
  xlim(0, 15)
ggplot() + 
  geom_point(UE.DBA.decidual.data.nor, mapping = aes(x = log2(iDSC_Ptprc_Normal+1), y = log2(iDSC_Ptprc_DBA.1+1)), size = 0.1)

D1.lm <- lm(log2(D1_Lum_Normal+1) ~ log2(D1_Lum_DBA.1+1), data = UE.DBA.decidual.data.nor)
summary(D1.lm) # R^2=0.50   0.87
D2.lm <- lm(log2(D2_Sfrp4_Normal+1) ~ log2(D2_Sfrp4_DBA.1+1), data = UE.DBA.decidual.data.nor)
summary(D2.lm) # R^2=0.68   0.91   
D3.lm <- lm(log2(D3_Top2a_Normal+1) ~ log2(D3_Top2a_DBA.1+1), data = UE.DBA.decidual.data.nor)
summary(D3.lm) # R^2=0.62   0.96
D4.lm <- lm(log2(D4_Ptn_Normal+1) ~ log2(D4_Ptn_DBA.1+1), data = UE.DBA.decidual.data.nor)
summary(D4.lm) # R^2=0.67   0.97
D5.lm <- lm(log2(D5_Gatm_Normal+1) ~ log2(D5_Gatm_DBA.1+1), data = UE.DBA.decidual.data.nor)
summary(D5.lm) # R^2=0.64   0.97
D6.lm <- lm(log2(D6_S100a8_Normal+1) ~ log2(D6_S100a8_DBA.1+1), data = UE.DBA.decidual.data.nor)
summary(D6.lm) # R^2=0.77   0.97
D7.lm <- lm(log2(D7_Prl8a2_Normal+1) ~ log2(D7_Prl8a2_DBA.1+1), data = UE.DBA.decidual.data.nor)
summary(D7.lm) # R^2=0.51   0.81
iDSC.lm <- lm(log2(iDSC_Ptprc_Normal+1) ~ log2(iDSC_Ptprc_DBA.1+1), data = UE.DBA.decidual.data.nor)
summary(iDSC.lm) # R^2=0.83   0.95

#subcluster DBA-DSC
load("20210827-figures/files/DBA.1-20220912.Rdata")
DBA.1.DSC <- subset(DBA.1, subset = seurat_clusters %in% c(1,4,6,8,13))
DBA.1.DSC$time <- "DBA.1"
DBA.1.DSC$cellname <- "DBA.1"
DBA.1.DSC <- NormalizeData(DBA.1.DSC, normalization.method = "LogNormalize", scale.factor = 10000)
DBA.1.DSC <- FindVariableFeatures(DBA.1.DSC, selection.method = "vst", nfeatures = 2000)
DBA.1.DSC <- ScaleData(DBA.1.DSC, features = rownames(DBA.1.DSC))
DBA.1.DSC <- RunPCA(DBA.1.DSC, features = VariableFeatures(DBA.1.DSC), npcs = 70)
ElbowPlot(DBA.1.DSC, ndims = 70)
DBA.1.DSC <- FindNeighbors(object = DBA.1.DSC, dims = 1:30)
DBA.1.DSC <- FindClusters(DBA.1.DSC, resolution = 0.4)
DBA.1.DSC <- RunUMAP(DBA.1.DSC, dims = 1:30, seed.use = 0912)
DimPlot(DBA.1.DSC, label = T)
DBA.1.DSC.markers <- FindAllMarkers(DBA.1.DSC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.1.DSC.markers.top5 <- DBA.1.DSC.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(DBA.1.DSC, features = DBA.1.DSC.markers.top5$gene)
DBA.1.DSC <- RenameIdents(DBA.1.DSC, "7" = "D1", "6" = "D2", "3" = "D3", "4" = "D3", "1" = "D4",
                          "2" = "D5", "0" = "D6", "8" = "D8", "5" = "iDSC_DBA1", "9" = "iDSC_DBA1")
DBA.1.DSC$cellcluster <- DBA.1.DSC@active.ident
DBA.1.DSC$sample <- "DBA.1"
save(DBA.1.DSC, file = "20210827-figures/files/DBA.1-DSC-20221125.Rdata")
#gene expression 
D1.genes <- c("Lum", "Dcn", "Postn", "Col3a1", "Col1a2", "Col1a1", "Igfbp6")
FeaturePlot(DBA.1.DSC, features = D1.genes, ncol = 7, pt.size = 0.1)
VlnPlot(DBA.1.DSC, features = D1.genes, ncol = 7, pt.size = 0)
FeaturePlot(UE.decidual, features = D1.genes, ncol = 7, pt.size = 0.1)

D2.genes <- c("Sfrp4", "Igfbp3", "Igfbp5", "Ifitm1", "Igf1", "Hsd11b2", "Mmp2")
FeaturePlot(DBA.1.DSC, features = D2.genes, ncol = 7, pt.size = 0.1)
VlnPlot(DBA.1.DSC, features = D2.genes, ncol = 7, pt.size = 0)
FeaturePlot(UE.decidual, features = D2.genes, ncol = 7, pt.size = 0.1)

D3.genes <- c("Hist1h2ap", "Hist1h1b", "Hmgb2", "Top2a", "Ube2c")
FeaturePlot(DBA.1.DSC, features = D3.genes, ncol = 7, pt.size = 0.1)
VlnPlot(DBA.1.DSC, features = D3.genes, ncol = 7, pt.size = 0)
FeaturePlot(UE.decidual, features = D3.genes, ncol = 7, pt.size = 0.1)

D4.genes <- c("Ptn", "Bmp2", "Hmgn5", "Wnt6", "Wnt5a", "Fmo1", "Hspd1")
FeaturePlot(DBA.1.DSC, features = D4.genes, ncol = 7, pt.size = 0.1)
VlnPlot(DBA.1.DSC, features = D4.genes, ncol = 7, pt.size = 0)
FeaturePlot(UE.decidual, features = D4.genes, ncol = 7, pt.size = 0.1)

D5.genes <- c("Gatm", "Anxa8", "Cxcl14", "Ass1", "Ano1", "Tdo2", "Ebpl", "Nampt")
FeaturePlot(DBA.1.DSC, features = D5.genes, ncol = 8, pt.size = 0.1)
VlnPlot(DBA.1.DSC, features = D5.genes, ncol = 8, pt.size = 0)
FeaturePlot(UE.decidual, features = D5.genes, ncol = 8, pt.size = 0.1)

D6.genes <- c("S100a8", "Angpt2", "Lgals3", "Il1r2", "Tnfrsf11b", "Psca")
FeaturePlot(DBA.1.DSC, features = D6.genes, ncol = 6, pt.size = 0.1)
VlnPlot(DBA.1.DSC, features = D6.genes, ncol = 6, pt.size = 0)
FeaturePlot(UE.decidual, features = D6.genes, ncol = 6, pt.size = 0.1)

D7.genes <- c("Procr", "Ada", "Erv3", "Htra1", "Prl8a2", "Psap", "Smpd1")
FeaturePlot(DBA.1.DSC, features = D7.genes, ncol = 7, pt.size = 0.1)
VlnPlot(DBA.1.DSC, features = D7.genes, ncol = 7, pt.size = 0)
FeaturePlot(UE.decidual, features = D7.genes, ncol = 7, pt.size = 0.1)

#20221122
#DBA.2 DSC subclsuter analysis
load("20210827-figures/files/DBA.2-20220912.Rdata")
DBA.2.DSC <- subset(DBA.2, subset = seurat_clusters %in% c(0, 1, 13))
DBA.2.DSC$time <- "DBA.2"
DBA.2.DSC$cellname <- "DBA.2.decidua"
DBA.2.DSC <- NormalizeData(DBA.2.DSC, normalization.method = "LogNormalize", scale.factor = 10000)
DBA.2.DSC <- FindVariableFeatures(DBA.2.DSC, selection.method = "vst", nfeatures = 2000)
DBA.2.DSC <- ScaleData(DBA.2.DSC, features = rownames(DBA.2.DSC))
DBA.2.DSC <- RunPCA(DBA.2.DSC, features = VariableFeatures(DBA.2.DSC), npcs = 70)
ElbowPlot(DBA.2.DSC, ndims = 70)
DBA.2.DSC <- FindNeighbors(object = DBA.2.DSC, dims = 1:30)
DBA.2.DSC <- RunUMAP(DBA.2.DSC, dims = 1:30, seed.use = 0912)
DBA.2.DSC <- FindClusters(DBA.2.DSC, resolution = 0.7)
DimPlot(DBA.2.DSC, label = T)
DBA.2.DSC.markers <- FindAllMarkers(DBA.2.DSC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.2.DSC.markers.top5 <- DBA.2.DSC.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(DBA.2.DSC, features = DBA.2.DSC.markers.top5$gene)
FeaturePlot(DBA.2.DSC, features = D7.genes, ncol = 8, pt.size = 0.1)
save(DBA.2.DSC, file = "20210827-figures/files/DBA.2-DSC-20221125.Rdata")

#integrate with Normal DSC
DBA.2.DSC$DBA_cluster <- paste("cluster_", DBA.2.DSC$seurat_clusters, sep = "")
load("FB-immune-DSC-harmony-20211214.RData")
UE.DBA.2.decidual <- merge(UE.FB, DBA.2.DSC)
UE.DBA.2.decidual <- NormalizeData(UE.DBA.2.decidual, normalization.method = "LogNormalize", scale.factor = 10000)
UE.DBA.2.decidual <- FindVariableFeatures(UE.DBA.2.decidual, selection.method = "vst", nfeatures = 2000)
UE.DBA.2.decidual <- ScaleData(UE.DBA.2.decidual, features = VariableFeatures(UE.DBA.2.decidual))
UE.DBA.2.decidual <- RunPCA(UE.DBA.2.decidual, features = VariableFeatures(UE.DBA.2.decidual))
UE.DBA.2.decidual <- RunHarmony(UE.DBA.2.decidual, "time", plot_convergence = TRUE, reduction = "pca")
DimPlot(object = UE.DBA.2.decidual, reduction = "harmony", pt.size = .1, group.by = "time")
ElbowPlot(UE.DBA.2.decidual, ndims = 80, reduction = "harmony")
UE.DBA.2.decidual <- RunUMAP(UE.DBA.2.decidual, reduction = "harmony", dims = 1:30)
UE.DBA.2.decidual <- FindNeighbors(UE.DBA.2.decidual, reduction = "harmony", dims = 1:30)
UE.DBA.2.decidual <- FindClusters(UE.DBA.2.decidual, resolution = 0.2)
DimPlot(UE.DBA.2.decidual, group.by = "DBA_cluster", label = T)
DimPlot(UE.DBA.2.decidual, group.by = "cellname", label = T)

#20221125 查看DSCs亚群、iDSC类群高表达的基因在Normal和DBA1中的表达情况
load("20210827-figures/files/DBA.1-DSC-20220912.Rdata")
DimPlot(UE.DBA.decidual)
#UE.DBA.decidual$sample <- gsub(".*-.*", "Normal", UE.DBA.decidual$cellname)
UE.DBA.decidual$sample_cluster <- paste(gsub(".*-.*", "Normal", UE.DBA.decidual$cellname), UE.DBA.decidual@active.ident, sep = "_")
UE.DBA.decidual.sub <- subset(UE.DBA.decidual, subset = sample_cluster %in% 
                                 c("Normal_D3_Top2a","Normal_D5_Gatm","Normal_D2_Sfrp4","Normal_D1_Lum","Normal_D4_Ptn","Normal_D7_Prl8a2",        
                                   "Normal_D6_S100a8","Normal_iDSC_Ptprc","DBA.decidua_D4_Ptn","DBA.decidua_D6_S100a8","DBA.decidua_D3_Top2a","DBA.decidua_iDSC_Ptprc",
                                   "DBA.decidua_D7_Prl8a2","DBA.decidua_D5_Gatm"))
UE.DBA.decidual.sub@active.ident <- factor(UE.DBA.decidual.sub$sample_cluster, 
                                           levels = c("Normal_D1_Lum","Normal_D2_Sfrp4","Normal_D3_Top2a","DBA.decidua_D3_Top2a",
                                                      "Normal_D4_Ptn","DBA.decidua_D4_Ptn","Normal_D5_Gatm","DBA.decidua_D5_Gatm",
                                                      "Normal_D6_S100a8","DBA.decidua_D6_S100a8","Normal_D7_Prl8a2","DBA.decidua_D7_Prl8a2",     
                                                      "Normal_iDSC_Ptprc","DBA.decidua_iDSC_Ptprc"))
UE.DBA.decidual.sub <- NormalizeData(UE.DBA.decidual.sub)
UE.DBA.decidual.sub <- ScaleData(UE.DBA.decidual.sub, features = rownames(UE.DBA.decidual.sub))
DSC.markers <- read.table("20210827-figures/files/UE-harmony-decidual-markers.csv", sep = ",", header = T, check.names = F, stringsAsFactors = F)
iDSC.markers <- read.table("20210827-figures/files/iDSC-compare-DSC-DEGs.csv", sep = ",", check.names = F, header = T, stringsAsFactors = F)
iDSC.markers.immune <- iDSC.markers[iDSC.markers$cluster == "immune", ]
DSC.iDSC.markers <- rbind(DSC.markers, iDSC.markers.immune)
DSC.iDSC.markers$cluster <- factor(DSC.iDSC.markers$cluster, levels = c("Lum-D", "Sfrp4-D", "Top2a-D", "Ptn-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "immune"))
DSC.iDSC.markers.top10 <- DSC.iDSC.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
pdf("20220622-Comments/05.DBA/00.figures/11.DSC-iDSC-markers-heatmap.pdf", width = 30, height = 20)
DoHeatmap(UE.DBA.decidual.sub, features = DSC.iDSC.markers.top10$gene[1:80]) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()
#add module
UE.DBA.decidual.sub <- AddModuleScore(UE.DBA.decidual.sub, features = list(DSC.iDSC.markers[DSC.iDSC.markers$cluster == "Lum-D", ]$gene), name = "D1.markers")
UE.DBA.decidual.sub <- AddModuleScore(UE.DBA.decidual.sub, features = list(DSC.iDSC.markers[DSC.iDSC.markers$cluster == "Sfrp4-D", ]$gene), name = "D2.markers")
UE.DBA.decidual.sub <- AddModuleScore(UE.DBA.decidual.sub, features = list(DSC.iDSC.markers[DSC.iDSC.markers$cluster == "Top2a-D", ]$gene), name = "D3.markers")
UE.DBA.decidual.sub <- AddModuleScore(UE.DBA.decidual.sub, features = list(DSC.iDSC.markers[DSC.iDSC.markers$cluster == "Ptn-D", ]$gene), name = "D4.markers")
UE.DBA.decidual.sub <- AddModuleScore(UE.DBA.decidual.sub, features = list(DSC.iDSC.markers[DSC.iDSC.markers$cluster == "Gatm-D", ]$gene), name = "D5.markers")
UE.DBA.decidual.sub <- AddModuleScore(UE.DBA.decidual.sub, features = list(DSC.iDSC.markers[DSC.iDSC.markers$cluster == "S100a8-D", ]$gene), name = "D6.markers")
UE.DBA.decidual.sub <- AddModuleScore(UE.DBA.decidual.sub, features = list(DSC.iDSC.markers[DSC.iDSC.markers$cluster == "Prl8a2-D", ]$gene), name = "D7.markers")
UE.DBA.decidual.sub <- AddModuleScore(UE.DBA.decidual.sub, features = list(DSC.iDSC.markers[DSC.iDSC.markers$cluster == "immune", ]$gene), name = "immune.markers")
FeaturePlot(UE.DBA.decidual.sub, features = "D1.markers1", split.by = "sample")
FeaturePlot(UE.DBA.decidual.sub, features = "D2.markers1", split.by = "sample")
FeaturePlot(UE.DBA.decidual.sub, features = "D3.markers1", split.by = "sample")
FeaturePlot(UE.DBA.decidual.sub, features = "D4.markers1", split.by = "sample")
FeaturePlot(UE.DBA.decidual.sub, features = "D5.markers1", split.by = "sample")
FeaturePlot(UE.DBA.decidual.sub, features = "D6.markers1", split.by = "sample")
FeaturePlot(UE.DBA.decidual.sub, features = "D7.markers1", split.by = "sample")
FeaturePlot(UE.DBA.decidual.sub, features = "immune.markers1", split.by = "sample")

Major.markers <- read.table("20210827-figures/files/UE-harmony-major-DEGs.csv", sep = ",", check.names = F, stringsAsFactors = F, header = T)
Major.markers.immune <- Major.markers[Major.markers$cluster %in% "Immune", ]
DSC.iDSC.markers.major <- rbind(DSC.markers, Major.markers.immune)
UE.DBA.decidual.sub <- AddModuleScore(UE.DBA.decidual.sub, features = list(DSC.iDSC.markers.major[DSC.iDSC.markers.major$cluster == "Immune", ]$gene), name = "major.immune.markers")
FeaturePlot(UE.DBA.decidual.sub, features = "major.immune.markers1", split.by = "sample")
VlnPlot(UE.DBA.decidual.sub, features = c("D1.markers1", "D2.markers1", "D3.markers1", "D4.markers1", 
                                          "D5.markers1", "D6.markers1", "D7.markers1", "immune.markers1"), group.by = "sample", pt.size = 0)

#直接在图中show出iDSC的三个亚群
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
#FB.immune.SCT
DSC.subcluster.meta <- data.frame(cellname = colnames(UE.DBA.decidual), 
                                   cell_subcluster = UE.DBA.decidual$cellname, stringsAsFactors = F)
DSC.subcluster.meta <- DSC.subcluster.meta[DSC.subcluster.meta$cell_subcluster != "FB-immune", ]
iDSC.subcluster.meta <- data.frame(cellname = colnames(FB.immune.SCT), 
                                   cell_subcluster = paste("iDSC", FB.immune.SCT@active.ident, sep = ""), stringsAsFactors = F)
DSC.iDSC.subcluster.meta <- rbind(DSC.subcluster.meta, iDSC.subcluster.meta)
rownames(DSC.iDSC.subcluster.meta) <- DSC.iDSC.subcluster.meta$cellname
DSC.iDSC.subcluster.meta$cellname <- NULL
UE.DBA.decidual <- AddMetaData(UE.DBA.decidual, metadata = DSC.iDSC.subcluster.meta)
save(UE.DBA.decidual, file = "20210827-figures/files/UE-DSC-iDSC_sub-DBA-20221126.Rdata")
UE.DBA.decidual$cell_subcluster <- factor(UE.DBA.decidual$cell_subcluster, levels = c("Lum-D", "Sfrp4-D", "Top2a-D", "Ptn-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "iDSC0", "iDSC1", "iDSC2", "DBA.decidua"))
highlight.cell <- subset(UE.DBA.decidual, subset = cell_subcluster == "DBA.decidua")
highlight.cell <- colnames(highlight.cell)
pdf("20220622-Comments/05.DBA/00.figures/12.DSC-iDSC-subclsuter-dimplot-DBAgrey.pdf", width = 11.5, height = 9.5)
DimPlot(UE.DBA.decidual, group.by = "cell_subcluster", cols = c(decidual.color[1:7], "#DB7E71", "#61B057", "#6E91C9", "grey50"), pt.size = 1)
dev.off()

#Add iDSC0/1/2三个subcluster的基因表达score
iDSC.sub.DEGs <- read.table("20210827-figures/files/iDSC-subclusters-markers.csv", sep = ",", check.names = F, header = T, stringsAsFactors = F)
UE.DBA.decidual <- AddModuleScore(UE.DBA.decidual, features = list(iDSC.sub.DEGs[iDSC.sub.DEGs$cluster == "0", ]$gene), name = "iDSC0")
UE.DBA.decidual <- AddModuleScore(UE.DBA.decidual, features = list(iDSC.sub.DEGs[iDSC.sub.DEGs$cluster == "1", ]$gene), name = "iDSC1")
UE.DBA.decidual <- AddModuleScore(UE.DBA.decidual, features = list(iDSC.sub.DEGs[iDSC.sub.DEGs$cluster == "2", ]$gene), name = "iDSC2")
VlnPlot(UE.DBA.decidual, features = c("iDSC01", "iDSC11", "iDSC21"), group.by = "sample", pt.size = 0)
UE.DBA.iDSC <- subset(UE.DBA.decidual, subset = cellcluster == 'iDSC_Ptprc')
UE.DBA.iDSC <- subset(UE.DBA.iDSC, subset = cell_subcluster %in% c("iDSC0", "iDSC1", "iDSC2", "DBA.decidua"))
UE.DBA.iDSC <- SCTransform(UE.DBA.iDSC, ncells = 3000, verbose = FALSE, return.only.var.genes = FALSE)
UE.DBA.iDSC <- RunPCA(UE.DBA.iDSC)
UE.DBA.iDSC <- RunHarmony(UE.DBA.iDSC, group.by.vars = "sample", plot_convergence = TRUE, reduction = "pca", assay.use = "SCT")
UE.DBA.iDSC <- RunUMAP(UE.DBA.iDSC, dims = 1:30, reduction = "harmony")
UE.DBA.iDSC <- FindNeighbors(UE.DBA.iDSC, dims = 1:30, reduction = "harmony")
UE.DBA.iDSC <- FindClusters(UE.DBA.iDSC, resolution = 0.2)
DimPlot(UE.DBA.iDSC, label = T)
DimPlot(UE.DBA.iDSC, group.by = "sample")
UE.DBA.iDSC.markers <- FindAllMarkers(UE.DBA.iDSC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
UE.DBA.iDSC.markers.top10 <- UE.DBA.iDSC.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(UE.DBA.iDSC, features = UE.DBA.iDSC.markers.top10$gene) + NoLegend()
UE.DBA.iDSC <- RenameIdents(UE.DBA.iDSC, "0" = "iDSC_0", "2" = "iDSC_1", "3" = "iDSC_2", "1" = "iDSC_DBA")
UE.DBA.iDSC$cellcluster <- UE.DBA.iDSC@active.ident
DimPlot(UE.DBA.iDSC, group.by = "cell_subcluster")

##20221127 immune cell subclusters
load("20210827-figures/files/DBA.1-20220912.Rdata")
DBA.1.immune <- subset(DBA.1, seurat_clusters %in% c(0,2,12,11,9))
DBA.1.immune <- NormalizeData(DBA.1.immune, normalization.method = "LogNormalize", scale.factor = 10000)
DBA.1.immune <- FindVariableFeatures(DBA.1.immune, selection.method = "vst", nfeatures = 2000)
DBA.1.immune <- ScaleData(DBA.1.immune, features = rownames(DBA.1.immune))
DBA.1.immune <- RunPCA(DBA.1.immune, features = VariableFeatures(DBA.1.immune), npcs = 70)
ElbowPlot(DBA.1.immune, ndims = 70)
DBA.1.immune <- FindNeighbors(object = DBA.1.immune, dims = 1:30)
DBA.1.immune <- FindClusters(DBA.1.immune, resolution = 0.3)
DBA.1.immune <- RunUMAP(DBA.1.immune, dims = 1:30, seed.use = 0912)
DimPlot(DBA.1.immune, label = T)
DBA.1.immune.markers <- FindAllMarkers(DBA.1.immune, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.1.immune.markers.top10 <- DBA.1.immune.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(DBA.1.immune, features = DBA.1.immune.markers.top10$gene)
DBA.1.immune <- RenameIdents(DBA.1.immune, "0" = "Monocyte", "2" = "Monocyte", "1" = "Mac", "3" = "Proliferating Mac", 
                            "5" = "NK", "8" = "T cell", "4" = "iDSC_DBA1", "7" = "Neutrophil", "6" = "DC-1")
DBA.1.immune$sample <- 'DBA.1'
DBA.1.immune$cellcluster <- DBA.1.immune@active.ident
save(DBA.1.immune, file = "20210827-figures/files/DBA.1-immune-20221130.RData")

load("20210827-figures/files/DBA.1-DSC-20221125.Rdata")
DBA.1.decidua.iDSC <- subset(DBA.1.DSC, subset = seurat_clusters %in% c(5))
DBA.1.decidua.iDSC$cellcluster <- "DSC.iDSC"
DBA.1.immune.m <- merge(DBA.1.immune, DBA.1.decidua.iDSC)
DBA.1.immune.m <- NormalizeData(DBA.1.immune.m, normalization.method = "LogNormalize", scale.factor = 10000)
DBA.1.immune.m <- ScaleData(DBA.1.immune.m, features = rownames(DBA.1.immune.m))
VariableFeatures(DBA.1.immune.m) <- VariableFeatures(DBA.1.immune)
DBA.1.immune.m <- RunPCA(DBA.1.immune.m, features = VariableFeatures(DBA.1.immune.m), npcs = 70)
ElbowPlot(DBA.1.immune.m, ndims = 70)
DBA.1.immune.m <- FindNeighbors(object = DBA.1.immune.m, dims = 1:20)
DBA.1.immune.m <- FindClusters(DBA.1.immune.m, resolution = 0.3)
DBA.1.immune.m <- RunUMAP(DBA.1.immune.m, dims = 1:20, seed.use = 0912)
DimPlot(DBA.1.immune.m, label = T)
DimPlot(DBA.1.immune.m, label = T, group.by = "cellcluster")
DBA.1.immune.m$cellcluster <- gsub("DSC.iDSC", "iDSC_DBA1", DBA.1.immune.m$cellcluster)
DBA.1.immune.m@active.ident <- factor(DBA.1.immune.m$cellcluster)
save(DBA.1.immune.m, file = "20210827-figures/files/DBA.1-immune-20221215.RData")

#integrate
load("20210827-figures/files/UE-immune-harmony-web.RData")
DimPlot(UE.immune.harmony)
UE.immune.harmony$sample <- "Normal"
UE.DBA.immune <- merge(UE.immune.harmony, DBA.1.immune)
UE.DBA.immune <- NormalizeData(UE.DBA.immune)
UE.DBA.immune <- ScaleData(UE.DBA.immune, features = rownames(UE.DBA.immune))
UE.DBA.immune <- FindVariableFeatures(UE.DBA.immune, nfeatures = 2000, selection.method = "vst")
UE.DBA.immune <- RunPCA(UE.DBA.immune, features = VariableFeatures(UE.DBA.immune))
UE.DBA.immune <- RunHarmony(UE.DBA.immune, "time", plot_convergence = TRUE, reduction = "pca")
DimPlot(object = UE.DBA.immune, reduction = "harmony", pt.size = .1, group.by = "time")
ElbowPlot(UE.DBA.immune, ndims = 80, reduction = "harmony")
UE.DBA.immune <- RunUMAP(UE.DBA.immune, reduction = "harmony", dims = 1:23)
UE.DBA.immune <- FindNeighbors(UE.DBA.immune, reduction = "harmony", dims = 1:23)
UE.DBA.immune <- FindClusters(UE.DBA.immune, resolution = 0.4)
DimPlot(UE.DBA.immune, label = T)
DimPlot(UE.DBA.immune, group.by = "sample", label = T)
DimPlot(UE.DBA.immune, group.by = "sub_clusters", label = T)
UE.DBA.immune.iDSC <- subset(UE.DBA.immune, subset = seurat_clusters %in% c("4"))
DimPlot(DBA.1, cells.highlight = colnames(UE.DBA.immune.iDSC), sizes.highlight = 0.1)

DBA.1.immune.iDSC <- subset(DBA.1.immune, subset = seurat_clusters == 4)
DBA.1.DSC.iDSC <- subset(DBA.1.DSC, subset = seurat_clusters %in% c(5, 9))
DBA.1.immune.iDSC$major <- "Immune"
DBA.1.DSC.iDSC$major <- "DSC"
DBA.1.iDSC <- merge(DBA.1.immune.iDSC, DBA.1.DSC.iDSC)
DBA.1.iDSC <- SCTransform(DBA.1.iDSC)
#DBA.1.iDSC <- FindVariableFeatures(DBA.1.iDSC, selection.method = "vst", nfeatures = 2000)
#DBA.1.iDSC <- ScaleData(DBA.1.iDSC, features = rownames(DBA.1.iDSC))
DBA.1.iDSC <- RunPCA(DBA.1.iDSC, features = VariableFeatures(DBA.1.iDSC), npcs = 70, assay = "SCT")
ElbowPlot(DBA.1.iDSC, ndims = 70)
DBA.1.iDSC <- FindNeighbors(object = DBA.1.iDSC, dims = 1:15, assay = "SCT")
DBA.1.iDSC <- RunUMAP(DBA.1.iDSC, dims = 1:15, seed.use = 0912, assay = "SCT")
DBA.1.iDSC <- FindClusters(DBA.1.iDSC, resolution = 0.5)
DimPlot(DBA.1.iDSC, label = T)
DimPlot(DBA.1.iDSC, label = T, group.by = "major")
DBA.1.iDSC.markers <- FindAllMarkers(DBA.1.iDSC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.1.iDSC.markers.top10 <- DBA.1.iDSC.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(DBA.1.iDSC, features = DBA.1.iDSC.markers.top10$gene)
FeaturePlot(DBA.1.iDSC, features = c("Ptprc", "Gzma", "Hand2", "Lyz2"))
DBA.1.iDSC <- RenameIdents(DBA.1.iDSC, "0" = "DBA.1.iDSC0", "1" = "DBA.1.iDSC1", 
                           "2" = "DBA.1.iDSC2", "3" = "DBA.1.iDSC3", "4" = "DBA.1.iDSC4")
save(DBA.1.iDSC, file = "20210827-figures/files/DBA.1-iDSC-subclusters-20221129.RData")
load("20210827-figures/files/UE-DSC-web.RData")
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221129.RData")
DBA.1.iDSC.decidual <- merge(UE.decidual.susbet, DBA.1.iDSC)
DBA.1.iDSC.decidual$sub_clusters[is.na(DBA.1.iDSC.decidual$sub_clusters)] <- "DBA.1.iDSC"
DBA.1.iDSC.decidual@active.ident <- factor(DBA.1.iDSC.decidual$sub_clusters)
DBA.1.iDSC.decidual.RNA <- AverageExpression(DBA.1.iDSC.decidual, assays = "RNA", slot = "data")
plclust(hclust(d = dist(x = 1-cor(DBA.1.iDSC.decidual.RNA$RNA))))

#merge with immune
load("20210827-figures/files/UE-immune-harmony-web.RData")
UE.immune.harmony.sub <- subset(UE.immune.harmony, sub_clusters != "iDSC")
DBA.1.iDSC.immune <- merge(UE.immune.harmony.sub, DBA.1.iDSC)
DBA.1.iDSC.immune.RNA <- AverageExpression(DBA.1.iDSC.immune, assays = "RNA", slot = "data")
DBA.1.iDSC.immune.RNA <- DBA.1.iDSC.immune.RNA$RNA
DBA.1.iDSC.immune.RNA.cor <- cor(DBA.1.iDSC.immune.RNA)
DBA.1.iDSC.immune.RNA.cor <- DBA.1.iDSC.immune.RNA.cor[c(3:7), c(1:2,8:14)]
pheatmap(DBA.1.iDSC.immune.RNA.cor, scale = "none", cluster_rows = F, cluster_cols = F)

load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
FB.immune.SCT$sub_clusters <- paste("Normal.iDSC", FB.immune.SCT$seurat_clusters, sep = "")
FB.immune.SCT@active.ident <- factor(FB.immune.SCT$sub_clusters)
DBA.1.iDSC.normal.iDSC <- merge(FB.immune.SCT, DBA.1.iDSC)
DBA.1.iDSC.normal.iDSC.RNA <- AverageExpression(DBA.1.iDSC.normal.iDSC, assays = "RNA", slot = "data")
DBA.1.iDSC.normal.iDSC.RNA <- DBA.1.iDSC.normal.iDSC.RNA$RNA
DBA.1.iDSC.normal.iDSC.RNA.cor <- cor(DBA.1.iDSC.normal.iDSC.RNA)
DBA.1.iDSC.normal.iDSC.RNA.cor <- DBA.1.iDSC.normal.iDSC.RNA.cor[c(1:5), c(6:8)]
pheatmap(DBA.1.iDSC.normal.iDSC.RNA.cor, scale = "none", cluster_rows = F, cluster_cols = F)

#20221204 iDSC DEG expression in DBA-iDSC 4 subclusters 
iDSC.normal.DEGs <- read.table("20210827-figures/files/iDSC-subclusters-markers.csv", sep = ",", header = T, check.names = F)
iDSC.normal.DEGs.top20 <- iDSC.normal.DEGs %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
DBA.1.iDSC <- NormalizeData(DBA.1.iDSC)
DBA.1.iDSC <- ScaleData(DBA.1.iDSC, features = rownames(DBA.1.iDSC))
DBA.1.iDSC$subcluster <- DBA.1.iDSC@active.ident
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221129.RData")
DoHeatmap(DBA.1.iDSC, features = as.character(iDSC.normal.DEGs.top20$gene))

load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
FB.immune.SCT$subcluster <- paste("Normal", FB.immune.SCT@active.ident, sep = "")
FB.immune.SCT@active.ident <- factor(FB.immune.SCT$subcluster)
DBA.Normal.iDSC <- merge(FB.immune.SCT, DBA.1.iDSC)
DBA.Normal.iDSC <- NormalizeData(DBA.Normal.iDSC)
DBA.Normal.iDSC <- ScaleData(DBA.Normal.iDSC, features = rownames(DBA.Normal.iDSC))
DBA.Normal.iDSC@active.ident <- factor(DBA.Normal.iDSC$subcluster)
DoHeatmap(DBA.Normal.iDSC, features = as.character(iDSC.normal.DEGs.top20$gene))

#integrate DBA.iDSC 4 subclusters and Normal.iDSC 3 subclusters
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221129.RData")
DimPlot(DBA.1.iDSC)
FB.immune.SCT$subcluster <- paste("Normal", FB.immune.SCT@active.ident, sep = "")
FB.immune.SCT$sample <- "Normal"
#CCA
features <- SelectIntegrationFeatures(object.list = list(FB.immune.SCT, DBA.1.iDSC), nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = list(FB.immune.SCT, DBA.1.iDSC), anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                  anchor.features = features)
iDSC.combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
iDSC.combined.sct <- RunPCA(iDSC.combined.sct, verbose = FALSE)
iDSC.combined.sct <- RunUMAP(iDSC.combined.sct, reduction = "pca", dims = 1:30)
iDSC.combined.sct <- FindNeighbors(iDSC.combined.sct, dims = 1:30)
iDSC.combined.sct <- FindClusters(iDSC.combined.sct, resolution = 0.2)
DimPlot(iDSC.combined.sct, label = T)
DimPlot(iDSC.combined.sct, group.by = "sample")
DimPlot(iDSC.combined.sct, group.by = "subcluster", split.by = "sample", label = T)

#20221212 
#subset DSC in DBA.1
load("20210827-figures/files/DBA.1-DSC-20221125.Rdata")
DBA.1.DSC@active.ident <- DBA.1.DSC$seurat_clusters
DBA.1.DSC <- subset(DBA.1.DSC, subset = seurat_clusters %in% c(0,1,2,3,4,6,8))#cluster7:Kdr+.  cluster9:Krt8+   cluster5:Ptprc+, iDSC
FeaturePlot(DBA.1.DSC, features = c("Prl8a2", "Krt8", "nFeature_RNA"))
DBA.1.DSC <- NormalizeData(DBA.1.DSC, normalization.method = "LogNormalize", scale.factor = 10000)
DBA.1.DSC <- FindVariableFeatures(DBA.1.DSC, selection.method = "vst", nfeatures = 2000)
DBA.1.DSC <- ScaleData(DBA.1.DSC, features = rownames(DBA.1.DSC))
DBA.1.DSC <- RunPCA(DBA.1.DSC, features = VariableFeatures(DBA.1.DSC), npcs = 70)
ElbowPlot(DBA.1.DSC, ndims = 70)
DBA.1.DSC <- FindNeighbors(object = DBA.1.DSC, dims = 1:20)
DBA.1.DSC <- FindClusters(DBA.1.DSC, resolution = 0.3)
DBA.1.DSC <- RunUMAP(DBA.1.DSC, dims = 1:20, seed.use = 0912)
DBA.1.DSC <- RenameIdents(DBA.1.DSC, "2" = "D3_Top2a", "1" = "D4_Ptn", "3" = "D5_Gatm", "4" = "D5_Gatm", 
                          "0" = "D6_S100a8", "5" = "D7_Prl8a2")
DimPlot(DBA.1.DSC, label = T)
DBA.1.DSC.markers <- FindAllMarkers(DBA.1.DSC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.1.DSC.markers.top5 <- DBA.1.DSC.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(DBA.1.DSC, features = DBA.1.DSC.markers.top5$gene)
DBA.1.DSC$cellcluster <- DBA.1.DSC@active.ident
save(DBA.1.DSC, file = "20210827-figures/files/DBA.1-DSC-20221212.Rdata")

#iDSC
DBA.1.immune.iDSC <- subset(DBA.1.immune, subset = seurat_clusters == 4)
DBA.1.DSC.iDSC <- subset(DBA.1.DSC, subset = seurat_clusters %in% c(5))
DBA.1.immune.iDSC$major <- "Immune"
DBA.1.DSC.iDSC$major <- "DSC"
DBA.1.iDSC <- merge(DBA.1.immune.iDSC, DBA.1.DSC.iDSC)
DBA.1.iDSC <- SCTransform(DBA.1.iDSC)
DBA.1.iDSC <- RunPCA(DBA.1.iDSC, npcs = 70, assay = "SCT")
ElbowPlot(DBA.1.iDSC, ndims = 70)
DBA.1.iDSC <- FindNeighbors(object = DBA.1.iDSC, dims = 1:15, assay = "SCT")
DBA.1.iDSC <- RunUMAP(DBA.1.iDSC, dims = 1:15, seed.use = 0912, assay = "SCT")
DBA.1.iDSC <- FindClusters(DBA.1.iDSC, resolution = 0.5)
DBA.1.iDSC.markers <- FindAllMarkers(DBA.1.iDSC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.1.iDSC.markers.top10 <- DBA.1.iDSC.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(DBA.1.iDSC, features = DBA.1.iDSC.markers.top10$gene)
FeaturePlot(DBA.1.iDSC, features = c("Ptprc", "Gzma", "Hand2", "Lyz2"))
DBA.1.iDSC <- RenameIdents(DBA.1.iDSC, "2" = "DBA.1.iDSC0-1", "1" = "DBA.1.iDSC0-2", 
                           "0" = "DBA.1.iDSC1", "3" = "DBA.1.iDSC1", "4" = "DBA.1.iDSC2")
DBA.1.iDSC$cellcluster <- DBA.1.iDSC@active.ident
save(DBA.1.iDSC, file = "20210827-figures/files/DBA.1-iDSC-subclusters-20221212.RData")
DimPlot(DBA.1.iDSC, label = T)
DimPlot(DBA.1.iDSC, label = T, group.by = "major")
DBA.1.iDSC$cellcluster <- paste("DBA.1.iDSC", DBA.1.iDSC$seurat_clusters, sep = "")

#计算DBA-iDSC cluster 与DSC的相似性
load("20210827-figures/files/UE-DSC-web.RData")
DBA.1.iDSC$sub_clusters <- "DBA.1.iDSC"
DBA.1.iDSC.decidual <- merge(UE.decidual.susbet, DBA.1.iDSC)
DBA.1.iDSC.decidual@active.ident <- factor(DBA.1.iDSC.decidual$sub_clusters)
DBA.1.iDSC.decidual.RNA <- AverageExpression(DBA.1.iDSC.decidual, assays = "RNA", slot = "data")
plclust(hclust(d = dist(x = 1-cor(DBA.1.iDSC.decidual.RNA$RNA))))

#计算DBA.1.iDSC的四个subclusters和D3-D6以及免疫亚群之间的相关性，考虑是否加上Normal的作为control
library(gggibbous)
load("20210827-figures/files/UE-immune-harmony-web.RData")
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221212.RData")
UE.immune.harmony.sub <- subset(UE.immune.harmony, subset = sub_clusters != "iDSC")
DBA.1.iDSC$sub_clusters <- DBA.1.iDSC$cellcluster
FB.immune.SCT$sub_clusters <- paste("Normal.iDSC", FB.immune.SCT@active.ident, sep = "")
UE.DSC.immune.DBA.iDSC.RNA <- merge(UE.decidual.susbet, list(UE.immune.harmony.sub, FB.immune.SCT, DBA.1.iDSC))
UE.DSC.immune.DBA.iDSC.RNA@active.ident <- factor(UE.DSC.immune.DBA.iDSC.RNA$sub_clusters)
UE.DSC.immune.DBA.iDSC.RNA <- AverageExpression(UE.DSC.immune.DBA.iDSC.RNA, assays = "RNA", slot = "data")
UE.DSC.immune.DBA.iDSC.cor <- data.frame(cor(UE.DSC.immune.DBA.iDSC.RNA$RNA), check.names = F)
save(UE.DSC.immune.DBA.iDSC.cor, file = "20220622-Comments/05.DBA/DSC-Normal_iDSC-DBA_iDSC-immune-correlation.RData")
UE.DSC.immune.DBA.iDSC.cor <- UE.DSC.immune.DBA.iDSC.cor[c(5:9,14:19,2), c(10:13,20:22)]
pheatmap(UE.DSC.immune.DBA.iDSC.cor, scale = "column", cluster_rows = F, cluster_cols = F)
UE.DSC.immune.DBA.iDSC.cor$Var1 <- rownames(UE.DSC.immune.DBA.iDSC.cor)
UE.DSC.immune.DBA.iDSC.cor.m <- melt(UE.DSC.immune.DBA.iDSC.cor, id.vars = "Var1")
pdf("20220622-Comments/05.DBA/00.figures/23.DSC-Normal_iDSC-DBA_iDSC-immune-correlation-moon-plot.pdf", width = 6, height = 6.5)
ggplot() +
  geom_moon(UE.DSC.immune.DBA.iDSC.cor.m, show.legend = T,
            mapping = aes(x = variable, y = Var1, fill = value, ratio = value^2)) + 
  scale_fill_gradientn(colors = rev(brewer.pal(n = 11, "Spectral"))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

load("20220622-Comments/05.DBA/DSC-Normal_iDSC-DBA_iDSC-immune-correlation.RData")
Normal.DBA.iDSC.cor <- UE.DSC.immune.DBA.iDSC.cor[10:13, 20:22]
Normal.DBA.iDSC.cor$Var1 <- rownames(Normal.DBA.iDSC.cor)
Normal.DBA.iDSC.cor.m <- melt(Normal.DBA.iDSC.cor, id.vars = "Var1")
pdf("20220622-Comments/05.DBA/00.figures/24.Normal_iDSC-DBA_iDSC-correlation-moon-plot.pdf", width = 3.7, height = 3.2)
ggplot() +
  geom_moon(Normal.DBA.iDSC.cor.m, show.legend = T,
            mapping = aes(x = variable, y = Var1, fill = value, ratio = value^2)) + 
  scale_fill_gradientn(colors = rev(brewer.pal(n = 11, "Spectral"))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

#Normal DSC亚群、Normal iDSC亚群、DBA iDSC亚群分别subset 500 cells做pseudotime analysis,monocle2
DSC.sub <- subset(UE.decidual.susbet, downsample = 500)
#FB.immune.SCT
#DBA.1.iDSC
library(monocle)
UE.DSC.immune.DBA.iDSC.pseudo <- merge(DSC.sub, list(FB.immune.SCT, DBA.1.iDSC))
DSC.iDSC.data <- as(as.matrix(UE.DSC.immune.DBA.iDSC.pseudo@assays$RNA@counts), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", UE.DSC.immune.DBA.iDSC.pseudo@meta.data)
fd <- data.frame(gene_short_name = row.names(DSC.iDSC.data), row.names = row.names(DSC.iDSC.data))
fd <- new("AnnotatedDataFrame", data = fd)
DSC.iDSC.cds <- newCellDataSet(DSC.iDSC.data,
                               phenoData = pd,
                               featureData = fd, 
                               lowerDetectionLimit = 0.5,
                               expressionFamily = negbinomial.size())
#过滤
DSC.iDSC.cds <- estimateSizeFactors(DSC.iDSC.cds)
DSC.iDSC.cds <- estimateDispersions(DSC.iDSC.cds)
DSC.iDSC.cds <- detectGenes(DSC.iDSC.cds, min_expr = 0.1)
expressed_genes <- rownames(subset(fData(DSC.iDSC.cds), num_cells_expressed >= 10))
length(expressed_genes)
#using variable genes
#UE.DSC.iDSC.susbet <- FindVariableFeatures(UE.DSC.iDSC.susbet, selection.method = "vst", nfeatures = 2000)
order.genes <- differentialGeneTest(DSC.iDSC.cds[expressed_genes, ],
                                    fullModelFormulaStr = "~sub_clusters",
                                    cores = 10)
order.genes.top2000 <- rownames(order.genes)[order(order.genes$qval)][1:2000]
DSC.iDSC.cds <- setOrderingFilter(DSC.iDSC.cds, ordering_genes = order.genes.top2000)
DSC.iDSC.cds <- reduceDimension(DSC.iDSC.cds, method = "DDRTree")
DSC.iDSC.cds <- orderCells(DSC.iDSC.cds)
plot_cell_trajectory(DSC.iDSC.cds, color_by = "State", show_state_number = F, cell_size = 0.5, show_branch_points = F)
table(DSC.iDSC.cds$State, DSC.iDSC.cds$sub_clusters)

#20221215 data arrange
load("20210827-figures/files/DBA.1-DSC-20221212.Rdata")
decidual.color <- c("#6A5ACD","#9ACD32","#FB8072","#8DD3C7","#08519C","#FDB462","#228B22","#FCCDE5","#D9D9D9")
pdf("20220622-Comments/05.DBA/00.figures/02.DBA1-deciuda-UMAP-color-cluster.pdf", width = 10, height = 10)
DimPlot(DBA.1.DSC, pt.size = 1, cols = decidual.color[c(3:7)])
dev.off()
DBA.1.DSC.markers <- FindAllMarkers(DBA.1.DSC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.1.DSC.markers.top10 <- DBA.1.DSC.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
pdf("20220622-Comments/05.DBA/00.figures/02.DBA1-deciuda-DEGs-top10-heatmap.pdf", width = 10, height = 10)
DoHeatmap(DBA.1.DSC, features = DBA.1.DSC.markers.top10$gene, group.colors = decidual.color[c(3:7)]) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()

load("20210827-figures/files/DBA.1-immune-20221215.RData")
immune.color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#2E8B57", "#CAB2D6", "#FF7F00", "#F4A460", "#6A3D9A")
DBA.1.immune.m$cellcluster <- factor(DBA.1.immune.m$cellcluster, levels = c("Monocyte", "Mac", "Proliferating Mac", "DC-1", "Neutrophil", "NK", "T cell", "iDSC_DBA1"))
DBA.1.immune.m@active.ident <- DBA.1.immune.m$cellcluster
pdf("20220622-Comments/05.DBA/00.figures/13.DBA1-immune-UMAP-color-cluster.pdf", width = 10, height = 9)
DimPlot(DBA.1.immune.m, pt.size = 1, cols = immune.color[c(1:4,6,7,9,10)])
dev.off()

#EC
#Normal EC compare with DBA EC
load("20220622-Comments/05.DBA/DBA-EC-clusters-20221215.RData")
DBA.1.EC$sample <- "DBA.1"
load("../20210201-integrate/UE-EC-harmony-pc30-rename.RData")
UE.EC$sample <- "Normal"
UE.Normal.EC <- merge(UE.EC, DBA.1.EC)
UE.Normal.EC@active.ident <- factor(UE.Normal.EC$sample)
UE.Normal.EC.markers <- FindAllMarkers(UE.Normal.EC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
UE.Normal.EC <- ScaleData(UE.Normal.EC, features = rownames(UE.Normal.EC))
pdf("20220622-Comments/05.DBA/00.figures/19.Normal-DBA1-EC-DEGs-heatmap.pdf", width = 10, height = 20)
DoHeatmap(UE.Normal.EC, features = UE.Normal.EC.markers$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()
UE.Normal.EC.avg <- AverageExpression(UE.Normal.EC)
UE.Normal.EC.avg <- UE.Normal.EC.avg$RNA
UE.Normal.EC.avg$change <- "no"
UE.Normal.EC.avg$gene <- rownames(UE.Normal.EC.avg)
DBA.1.up <- UE.Normal.EC.markers[UE.Normal.EC.markers$cluster == "DBA.1", ]$gene
Normal.up <- UE.Normal.EC.markers[UE.Normal.EC.markers$cluster == "Normal", ]$gene
UE.Normal.EC.avg[UE.Normal.EC.avg$gene %in% DBA.1.up, ]$change <- "DBA.1"
UE.Normal.EC.avg[UE.Normal.EC.avg$gene %in% Normal.up, ]$change <- "Normal"
UE.Normal.EC.avg.anno <- UE.Normal.EC.avg[UE.Normal.EC.avg$change != "no", ]
pdf("20220622-Comments/05.DBA/00.figures/21-Normal-DBA-EC-compare-dotplot.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(UE.Normal.EC.avg, mapping = aes(x = log2(Normal + 1), y = log2(DBA.1 + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  theme_classic()
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/21-Normal-DBA-EC-compare-dotplot-annotate.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(UE.Normal.EC.avg, mapping = aes(x = log2(Normal + 1), y = log2(DBA.1 + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  annotate(geom = "text", x = log2(UE.Normal.EC.avg.anno$Normal+1), y = log2(UE.Normal.EC.avg.anno$DBA.1+1), 
           label = UE.Normal.EC.avg.anno$gene, size = 2) + 
  theme_classic()
dev.off()

UE.Normal.EC.GO.df <- data.frame()
for (c in unique(UE.Normal.EC.markers$cluster)) {
  print(c)
  geneset <- UE.Normal.EC.markers[UE.Normal.EC.markers$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  UE.Normal.EC.GO.df <- rbind(UE.Normal.EC.GO.df, tmp)
}
UE.Normal.EC.GO.df$GO_Term <- factor(UE.Normal.EC.GO.df$GO_Term, levels = rev(unique(UE.Normal.EC.GO.df$GO_Term)))
UE.Normal.EC.GO.df.top5 <- UE.Normal.EC.GO.df %>% group_by(cluster) %>% top_n(n = 10, -log(pvalue))
pdf("20220622-Comments/05.DBA/00.figures/20.Normal-DBA1-EC-DEGs-GO-top10.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(UE.Normal.EC.GO.df.top5, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()

#DBA.1.EC reclustering
DBA.1.EC <- subset(DBA.1, subset = seurat_clusters == 3)
DBA.1.EC <- NormalizeData(DBA.1.EC, normalization.method = "LogNormalize", scale.factor = 10000)
DBA.1.EC <- FindVariableFeatures(DBA.1.EC, selection.method = "vst", nfeatures = 2000)
DBA.1.EC <- ScaleData(DBA.1.EC, features = rownames(DBA.1.EC))
DBA.1.EC <- RunPCA(DBA.1.EC, features = VariableFeatures(DBA.1.EC), npcs = 70)
ElbowPlot(DBA.1.EC, ndims = 70)
DBA.1.EC <- FindNeighbors(object = DBA.1.EC, dims = 1:30)
DBA.1.EC <- RunUMAP(DBA.1.EC, dims = 1:30, seed.use = 0912)
DBA.1.EC <- FindClusters(DBA.1.EC, resolution = 0.3)
DimPlot(DBA.1.EC, label = T)
DBA.1.EC.markers <- FindAllMarkers(DBA.1.EC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
write.table(DBA.1.EC.markers, file = "20220622-Comments/05.DBA/CBA-EC-subcluster-markers.csv", sep = ",", quote = F)
DBA.1.EC.markers.top5 <- DBA.1.EC.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
DoHeatmap(DBA.1.EC, features = DBA.1.EC.markers.top5$gene)
DBA.1.EC <- RenameIdents(DBA.1.EC, "2" = "Proliterating EC", "0" = "Angiogenic EC", "6" = "Arterial EC", "3"= "Venous EC1", "1" = "Venous EC1", "4" = "Venous EC2", "5" = "EPC")
DBA.1.EC$cellcluster <- DBA.1.EC@active.ident
save(DBA.1.EC, file = "20220622-Comments/05.DBA/DBA-EC-clusters-20221215.RData")
EC.color <- c("#D95F02", "#7B68EE", "#E7298A", "#1B9E77", "#00BFFF", "#00FA9A")
pdf("20220622-Comments/05.DBA/00.figures/09.DBA1-EC-UMAP-color-cluster.pdf", width = 10, height = 9)
DimPlot(DBA.1.EC, pt.size = 1, cols = EC.color)
dev.off()
DBA.1.EC.markers.top10 <- DBA.1.EC.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
pdf("20220622-Comments/05.DBA/00.figures/09.DBA1-EC-DEGs-heatmap.pdf", width = 10, height = 9)
DoHeatmap(DBA.1.EC, features = DBA.1.EC.markers.top10$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()

DBA.EC.GO.df <- data.frame()
for (c in unique(DBA.1.EC.markers$cluster)) {
  print(c)
  geneset <- DBA.1.EC.markers[DBA.1.EC.markers$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  DBA.EC.GO.df <- rbind(DBA.EC.GO.df, tmp)
}
DBA.EC.GO.df$GO_Term <- factor(DBA.EC.GO.df$GO_Term, levels = rev(unique(DBA.EC.GO.df$GO_Term)))
DBA.EC.GO.df$cluster <- factor(DBA.EC.GO.df$cluster, levels = c("Proliterating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2", "EPC"))
DBA.EC.GO.df.top10 <- DBA.EC.GO.df %>% group_by(cluster) %>% top_n(n = 10, -log(pvalue))
pdf("20220622-Comments/05.DBA/00.figures/09.DBA1-EC-GO-terms-barplot-top10.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(DBA.EC.GO.df.top10, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  scale_fill_manual(values = EC.color) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()
#Inflammation EC DEGs module score
DBA.1.EC <- AddModuleScore(DBA.1.EC, features = list(DBA.1.EC.markers[DBA.1.EC.markers$cluster == "EPC", ]$gene), name = "Inflammation")
pdf("20220622-Comments/05.DBA/00.figures/09.DBA1-EC-Inflammation-module-score-featureplot.pdf", width = 9.5, height = 9)
FeaturePlot(DBA.1.EC, features = "Inflammation1", pt.size = 1) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

# SingleR annoatation
library("SingleR")
cell.EC <- read.csv("Cell-EC-normalize/Data.csv", header = T, row.names = 1)
cell.EC.t <- data.frame(t(as.matrix(cell.EC)), check.names = F)
cell.EC.t$Observation <- rownames(cell.EC.t)
cell.EC.meta <- read.csv("Cell-EC-normalize/Metadata.csv")
cell.EC.m <- merge(cell.EC.meta, cell.EC.t, by = "Observation", all = T)
cell.EC.m <- aggregate(cell.EC.m[, -c(1:3)], by = list(cell.EC.m$Cluster), mean)
rownames(cell.EC.m) <- cell.EC.m$Group.1
cell.EC.m <- cell.EC.m[, -1]
cell.EC.m <- t(cell.EC.m)
cell.EC.sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts=cell.EC.m))
colData(cell.EC.sce)$Type <- colnames(cell.EC.sce)
DBA.1.EC.data <- GetAssayData(DBA.1.EC, slot = "data")
DBA.1.EC.pred <- SingleR(test = DBA.1.EC.data, ref = cell.EC.sce, labels = cell.EC.sce$Type)
DBA.1.EC.singleR.df <- data.frame(table(DBA.1.EC.pred$labels, DBA.1.EC@active.ident), check.names = F)
DBA.1.EC.singleR.df$sample <- "DBA.1"
UE.EC.singleR.df <- data.frame(table(UE.EC.pred$labels, UE.EC@active.ident), check.names = F)
UE.EC.singleR.df$sample <- "Normal"
UE.DBA.1.EC.singleR.df <- rbind(UE.EC.singleR.df, DBA.1.EC.singleR.df)
UE.DBA.1.EC.singleR.df$Var1 <- gsub("capillary.*", "capillary", UE.DBA.1.EC.singleR.df$Var1)
UE.DBA.1.EC.singleR.df$Var1 <- gsub(".*artery.*", "artery", UE.DBA.1.EC.singleR.df$Var1)
UE.DBA.1.EC.singleR.df$Var1 <- gsub(".*vein.*", "vein", UE.DBA.1.EC.singleR.df$Var1)
ggplot() + 
  geom_bar(UE.DBA.1.EC.singleR.df, mapping = aes(x = sample, y = Freq, fill = Var1), stat = "identity", position = "fill") + 
  theme_classic()

#iDSC
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221212.RData")
iDSC.color <- c("#DB7093", "#BD7FF8", "#87AD34", "#56BCC2")
pdf("20220622-Comments/05.DBA/00.figures/14.DBA-iDSC-UMAP-color-cluster.pdf", width = 10, height = 9)
DimPlot(DBA.1.iDSC, pt.size = 2, cols = iDSC.color)
dev.off()
DBA.1.iDSC.markers <- FindAllMarkers(DBA.1.iDSC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.1.iDSC.markers.top10 <- DBA.1.iDSC.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
pdf("20220622-Comments/05.DBA/00.figures/25.DBA-iDSC-subclusters-DEGs-top10-heatmap.pdf", width = 10, height = 10)
DoHeatmap(DBA.1.iDSC, features = DBA.1.iDSC.markers.top10$gene, group.colors = iDSC.color) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()

DBA.iDSC.GO.df <- data.frame()
for (c in unique(DBA.1.iDSC.markers$cluster)) {
  print(c)
  geneset <- DBA.1.iDSC.markers[DBA.1.iDSC.markers$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  DBA.iDSC.GO.df <- rbind(DBA.iDSC.GO.df, tmp)
}
DBA.iDSC.GO.df$GO_Term <- factor(DBA.iDSC.GO.df$GO_Term, levels = rev(unique(DBA.iDSC.GO.df$GO_Term)))
DBA.iDSC.GO.df.top5 <- DBA.iDSC.GO.df %>% group_by(cluster) %>% top_n(n = 5, -log(pvalue))
pdf("20220622-Comments/05.DBA/00.figures/26.DBA-iDSC-subclusters-GO-top5.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(DBA.iDSC.GO.df.top5, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  scale_fill_manual(values = iDSC.color) +
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()


#merge Normal iDSC and DBA iDSC
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221212.RData")
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
FB.immune.SCT$cellcluster <- FB.immune.SCT$subcluster
DBA.1.iDSC.normal.iDSC <- merge(FB.immune.SCT, DBA.1.iDSC)
DBA.1.iDSC.normal.iDSC@active.ident <- factor(DBA.1.iDSC.normal.iDSC$cellcluster)
DBA.1.iDSC.normal.iDSC <- NormalizeData(DBA.1.iDSC.normal.iDSC)
DBA.1.iDSC.normal.iDSC <- ScaleData(DBA.1.iDSC.normal.iDSC, features = rownames(DBA.1.iDSC.normal.iDSC))
DBA.1.iDSC.normal.iDSC.markers <- FindAllMarkers(DBA.1.iDSC.normal.iDSC, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
DBA.1.iDSC.normal.iDSC.markers.top20 <- DBA.1.iDSC.normal.iDSC.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
pdf("20220622-Comments/05.DBA/00.figures/26.Normal-DBA-iDSC-subclusters-DEGs-top20-heatmap.pdf", width = 10, height = 10)
DoHeatmap(DBA.1.iDSC.normal.iDSC, features = DBA.1.iDSC.normal.iDSC.markers.top20$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()

#cluster
load("20210827-figures/files/UE-DSC-web.RData")
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221129.RData")
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
FB.immune.SCT <- NormalizeData(FB.immune.SCT)
FB.immune.SCT <- ScaleData(FB.immune.SCT, features = rownames(FB.immune.SCT))
DefaultAssay(FB.immune.SCT) <- "RNA"
FB.immune.SCT$sub_clusters <- "Normal.iDSC"
DBA.1.iDSC.decidual <- merge(UE.decidual.susbet, list(DBA.1.iDSC, FB.immune.SCT))
DBA.1.iDSC.decidual$sub_clusters[is.na(DBA.1.iDSC.decidual$sub_clusters)] <- "DBA.1.iDSC"
DBA.1.iDSC.decidual@active.ident <- factor(DBA.1.iDSC.decidual$sub_clusters)
DBA.1.iDSC.decidual.RNA <- AverageExpression(DBA.1.iDSC.decidual, assays = "RNA", slot = "data")
pdf("20220622-Comments/05.DBA/00.figures/22.DSC-Normal_iDSC-DBA_iDSC-cluster.pdf", width = 10, height = 10)
plclust(hclust(d = dist(x = 1-cor(DBA.1.iDSC.decidual.RNA$RNA))))
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/22.DSC-DBA_iDSC-cluster.pdf", width = 10, height = 10)
plclust(hclust(d = dist(x = 1-cor(DBA.1.iDSC.decidual.RNA$RNA[, -9]))))
dev.off()

#iDSC1 and iDSC2 between Normal and DBA comparing
DBA.1.iDSC.normal.iDSC.avg <- AverageExpression(DBA.1.iDSC.normal.iDSC, assays = "RNA")
DBA.1.iDSC.normal.iDSC.avg <- DBA.1.iDSC.normal.iDSC.avg$RNA
DBA.1.iDSC.normal.iDSC.avg$gene <- rownames(DBA.1.iDSC.normal.iDSC.avg)
#iDSC0
iDSC0.DEGs <- FindMarkers(DBA.1.iDSC.normal.iDSC, ident.1 = "Normal0", ident.2 = "DBA.1.iDSC0-1", min.pct = 0.5, logfc.threshold = 0.5)
iDSC0.DEGs$cluster <- "Normal0"
iDSC0.DEGs$gene <- rownames(iDSC0.DEGs)
iDSC0.DEGs[iDSC0.DEGs$avg_logFC<0, ]$cluster <- "DBA0"
iDSC0.RNA.avg <- DBA.1.iDSC.normal.iDSC.avg[, c(8,5,1)]
iDSC0.RNA.avg$change <- "NO"
iDSC0.RNA.avg[iDSC0.RNA.avg$gene %in% rownames(iDSC0.DEGs[iDSC0.DEGs$avg_logFC>0, ]), ]$change <- "Normal0"
iDSC0.RNA.avg[iDSC0.RNA.avg$gene %in% rownames(iDSC0.DEGs[iDSC0.DEGs$avg_logFC<0, ]), ]$change <- "DBA0"
iDSC0.RNA.avg.anno <- iDSC0.RNA.avg[iDSC0.RNA.avg$change != "NO", ]

pdf("20220622-Comments/05.DBA/00.figures/27.Normal-DBA-iDSC0-compare-dotplot.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(iDSC0.RNA.avg, mapping = aes(x = log2(Normal0 + 1), y = log2(`DBA.1.iDSC0-1` + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  theme_classic()
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/27.Normal-DBA-iDSC0-compare-dotplot-annotate.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(iDSC0.RNA.avg, mapping = aes(x = log2(Normal0 + 1), y = log2(`DBA.1.iDSC0-1` + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  annotate(geom = "text", x = log2(iDSC0.RNA.avg.anno$Normal0+1), y = log2(iDSC0.RNA.avg.anno$`DBA.1.iDSC0-1`+1), 
           label = iDSC0.RNA.avg.anno$gene, size = 2) + 
  theme_classic()
dev.off()

iDSC0.RNA.GO.df <- data.frame()
for (c in unique(iDSC0.DEGs$cluster)) {
  print(c)
  geneset <- iDSC0.DEGs[iDSC0.DEGs$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  iDSC0.RNA.GO.df <- rbind(iDSC0.RNA.GO.df, tmp)
}
iDSC0.RNA.GO.df$GO_Term <- factor(iDSC0.RNA.GO.df$GO_Term, levels = rev(unique(iDSC0.RNA.GO.df$GO_Term)))
iDSC0.RNA.GO.df.top5 <- iDSC0.RNA.GO.df %>% group_by(cluster) %>% top_n(n = 10, -log(pvalue))

pdf("20220622-Comments/05.DBA/00.figures/28.Normal-DBA1-iDSC0-DEGs-GO-top10.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(iDSC0.RNA.GO.df.top5, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()

#iDSC1
iDSC1.DEGs <- FindMarkers(DBA.1.iDSC.normal.iDSC, ident.1 = "Normal1", ident.2 = "DBA.1.iDSC1", min.pct = 0.5, logfc.threshold = 0.5)
iDSC1.DEGs$cluster <- "Normal1"
iDSC1.DEGs$gene <- rownames(iDSC1.DEGs)
iDSC1.DEGs[iDSC1.DEGs$avg_logFC<0, ]$cluster <- "DBA1"
iDSC1.RNA.avg <- DBA.1.iDSC.normal.iDSC.avg[, c(8,6,3)]
iDSC1.RNA.avg$change <- "NO"
iDSC1.RNA.avg[iDSC1.RNA.avg$gene %in% rownames(iDSC1.DEGs[iDSC1.DEGs$avg_logFC>0, ]), ]$change <- "Normal1"
iDSC1.RNA.avg[iDSC1.RNA.avg$gene %in% rownames(iDSC1.DEGs[iDSC1.DEGs$avg_logFC<0, ]), ]$change <- "DBA1"
iDSC1.RNA.avg.anno <- iDSC1.RNA.avg[iDSC1.RNA.avg$change != "NO", ]

pdf("20220622-Comments/05.DBA/00.figures/27.Normal-DBA-iDSC1-compare-dotplot.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(iDSC1.RNA.avg, mapping = aes(x = log2(Normal1 + 1), y = log2(DBA.1.iDSC1 + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  theme_classic()
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/27.Normal-DBA-iDSC1-compare-dotplot-annotate.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(iDSC1.RNA.avg, mapping = aes(x = log2(Normal1 + 1), y = log2(DBA.1.iDSC1 + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  annotate(geom = "text", x = log2(iDSC1.RNA.avg.anno$Normal1+1), y = log2(iDSC1.RNA.avg.anno$DBA.1.iDSC1+1), 
           label = iDSC1.RNA.avg.anno$gene, size = 2) + 
  theme_classic()
dev.off()

iDSC1.RNA.GO.df <- data.frame()
for (c in unique(iDSC1.DEGs$cluster)) {
  print(c)
  geneset <- iDSC1.DEGs[iDSC1.DEGs$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  iDSC1.RNA.GO.df <- rbind(iDSC1.RNA.GO.df, tmp)
}
iDSC1.RNA.GO.df$GO_Term <- factor(iDSC1.RNA.GO.df$GO_Term, levels = rev(unique(iDSC1.RNA.GO.df$GO_Term)))
iDSC1.RNA.GO.df.top5 <- iDSC1.RNA.GO.df %>% group_by(cluster) %>% top_n(n = 10, -log(pvalue))
pdf("20220622-Comments/05.DBA/00.figures/28.Normal-DBA1-iDSC1-DEGs-GO-top10.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(iDSC1.RNA.GO.df.top5, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()

#iDSC2
iDSC2.DEGs <- FindMarkers(DBA.1.iDSC.normal.iDSC, ident.1 = "Normal2", ident.2 = "DBA.1.iDSC2", min.pct = 0.5, logfc.threshold = 0.5)
iDSC2.DEGs$cluster <- "Normal2"
iDSC2.DEGs$gene <- rownames(iDSC2.DEGs)
iDSC2.DEGs[iDSC2.DEGs$avg_logFC<0, ]$cluster <- "DBA2"
iDSC2.RNA.avg <- DBA.1.iDSC.normal.iDSC.avg[, c(8,7,4)]
iDSC2.RNA.avg$change <- "NO"
iDSC2.RNA.avg[iDSC2.RNA.avg$gene %in% rownames(iDSC2.DEGs[iDSC2.DEGs$avg_logFC>0, ]), ]$change <- "Normal2"
iDSC2.RNA.avg[iDSC2.RNA.avg$gene %in% rownames(iDSC2.DEGs[iDSC2.DEGs$avg_logFC<0, ]), ]$change <- "DBA2"
iDSC2.RNA.avg.anno <- iDSC2.RNA.avg[iDSC2.RNA.avg$change != "NO", ]

pdf("20220622-Comments/05.DBA/00.figures/27.Normal-DBA-iDSC2-compare-dotplot.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(iDSC2.RNA.avg, mapping = aes(x = log2(Normal2 + 1), y = log2(DBA.1.iDSC2 + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  theme_classic()
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/27.Normal-DBA-iDSC2-compare-dotplot-annotate.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(iDSC2.RNA.avg, mapping = aes(x = log2(Normal2 + 1), y = log2(DBA.1.iDSC2 + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  annotate(geom = "text", x = log2(iDSC2.RNA.avg.anno$Normal2+1), y = log2(iDSC2.RNA.avg.anno$DBA.1.iDSC2+1), 
           label = iDSC2.RNA.avg.anno$gene, size = 2) + 
  theme_classic()
dev.off()

iDSC2.RNA.GO.df <- data.frame()
for (c in unique(iDSC2.DEGs$cluster)) {
  print(c)
  geneset <- iDSC2.DEGs[iDSC2.DEGs$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  iDSC2.RNA.GO.df <- rbind(iDSC2.RNA.GO.df, tmp)
}
iDSC2.RNA.GO.df$GO_Term <- factor(iDSC2.RNA.GO.df$GO_Term, levels = rev(unique(iDSC2.RNA.GO.df$GO_Term)))
iDSC2.RNA.GO.df.top5 <- iDSC2.RNA.GO.df %>% group_by(cluster) %>% top_n(n = 10, -log(pvalue))

pdf("20220622-Comments/05.DBA/00.figures/28.Normal-DBA1-iDSC2-DEGs-GO-top10.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(iDSC2.RNA.GO.df.top5, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()

#iDSC (think as one celltype, not subclusters) difference in Normal and DBA
load("20210827-figures/files/DBA.1-iDSC-subclusters-20221212.RData")
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
FB.immune.SCT$sample <- "Normal"
DBA.1.iDSC.normal.iDSC <- merge(FB.immune.SCT, DBA.1.iDSC)
DBA.1.iDSC.normal.iDSC@active.ident <- factor(DBA.1.iDSC.normal.iDSC$sample)
DBA.1.iDSC.normal.iDSC <- ScaleData(DBA.1.iDSC.normal.iDSC, features = rownames(DBA.1.iDSC.normal.iDSC))
DBA.1.iDSC.normal.iDSC.markers <- FindAllMarkers(DBA.1.iDSC.normal.iDSC, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
write.table(DBA.1.iDSC.normal.iDSC.markers, file = '20220622-Comments/05.DBA/Normal-DBA-iDSC-DEGs.csv', sep = ",", quote = F, row.names = F)
DBA.1.iDSC.normal.iDSC.avg <- AverageExpression(DBA.1.iDSC.normal.iDSC, assays = "RNA")
DBA.1.iDSC.normal.iDSC.avg <- DBA.1.iDSC.normal.iDSC.avg$RNA
DBA.1.iDSC.normal.iDSC.avg$gene <- rownames(DBA.1.iDSC.normal.iDSC.avg)
DBA.1.iDSC.normal.iDSC.avg$change <- "NO"
DBA.1.iDSC.normal.iDSC.avg[DBA.1.iDSC.normal.iDSC.avg$gene %in% DBA.1.iDSC.normal.iDSC.markers[DBA.1.iDSC.normal.iDSC.markers$cluster=="Normal", ]$gene, ]$change <- "Normal"
DBA.1.iDSC.normal.iDSC.avg[DBA.1.iDSC.normal.iDSC.avg$gene %in% DBA.1.iDSC.normal.iDSC.markers[DBA.1.iDSC.normal.iDSC.markers$cluster=="DBA.1", ]$gene, ]$change <- "DBA.1"
DBA.1.iDSC.normal.iDSC.avg.anno <- DBA.1.iDSC.normal.iDSC.avg[DBA.1.iDSC.normal.iDSC.avg$change != "NO", ]
pdf("20220622-Comments/05.DBA/00.figures/29.Normal-DBA-iDSC-all-compare-DEGs-heatmap.pdf", width = 10, height = 10)
DoHeatmap(DBA.1.iDSC.normal.iDSC, features = DBA.1.iDSC.normal.iDSC.markers$gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()
DBA.1.iDSC.normal.iDSC.markers.top20 <- DBA.1.iDSC.normal.iDSC.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
selection.gene <- c("Rpl29", "Fez1", "Lyz2", "Cxcl2", "Ccl9", "Ptn", "C1qa", "Ifi44", "F13a1", "Fcgr2b", "S100a4", 
                    "C1qb", "Cd83", "Crabp1", "Ccl8", "Lars2", "Cda", "Bst2", "Ear2", "Il1b", 
                    "A2m", "Erv3", "Tgfbr2", "Smoc2", "Gzmd", "Rgcc", "Anxa1", "Col4a2", 
                    "Klf4", "Rock2", "Spp1", "Dcn", "Angpt2", "Klf2", "Cdh5", "Gzmf", "Tnfrsf11b", "Fabp4", "Gzmg", "Htra3")
pdf("20220622-Comments/05.DBA/00.figures/29.Normal-DBA-iDSC-all-compare-DEGs-heatmap-selection20.pdf", width = 10, height = 10)
DoHeatmap(DBA.1.iDSC.normal.iDSC, features = selection.gene) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/29.Normal-DBA-iDSC-all-compare-dotplot.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(DBA.1.iDSC.normal.iDSC.avg, mapping = aes(x = log2(Normal + 1), y = log2(DBA.1 + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  theme_classic()
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/29.Normal-DBA-iDSC-all-compare-dotplot-annotate.pdf", width = 6, height = 5.5)
ggplot() + 
  geom_point(DBA.1.iDSC.normal.iDSC.avg, mapping = aes(x = log2(Normal + 1), y = log2(DBA.1 + 1), color = change), size = 1) + 
  scale_color_manual(values = c("#E87D72", "lightgray", "#6F9BF8")) + 
  annotate(geom = "text", x = log2(DBA.1.iDSC.normal.iDSC.avg.anno$Normal+1), y = log2(DBA.1.iDSC.normal.iDSC.avg.anno$DBA.1+1), 
           label = DBA.1.iDSC.normal.iDSC.avg.anno$gene, size = 2) + 
  theme_classic()
dev.off()

iDSC.RNA.GO.df <- data.frame()
for (c in unique(DBA.1.iDSC.normal.iDSC.markers$cluster)) {
  print(c)
  geneset <- DBA.1.iDSC.normal.iDSC.markers[DBA.1.iDSC.normal.iDSC.markers$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    GeneRation = enrich.GO$GeneRatio, 
                    gene = enrich.GO$geneID)
  iDSC.RNA.GO.df <- rbind(iDSC.RNA.GO.df, tmp)
}
iDSC.RNA.GO.df$GO_Term <- factor(iDSC.RNA.GO.df$GO_Term, levels = rev(unique(iDSC.RNA.GO.df$GO_Term)))
iDSC.RNA.GO.df.top5 <- iDSC.RNA.GO.df %>% group_by(cluster) %>% top_n(n = 10, -log(pvalue))
write.table(iDSC.RNA.GO.df, file = "20220622-Comments/05.DBA/DBA-Normal-iDSC-GO.csv", sep = ",", col.names = T, row.names = F, quote = F)
pdf("20220622-Comments/05.DBA/00.figures/30.Normal-DBA1-iDSC-all-DEGs-GO-top10.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(iDSC.RNA.GO.df.top5, mapping = aes(x = -log10(pvalue), y = GO_Term, fill = cluster), stat = "identity", position = "dodge") + 
  facet_wrap(~ cluster, scales = "free", ncol = 1) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
dev.off()

#cell proportion
load("20210827-figures/files/DBA.1-DSC-20221212.Rdata")
DimPlot(DBA.1.DSC, label = T)
DBA.1.DSC.proportion <- data.frame(table(DBA.1.DSC@active.ident), check.names = F)
DSC.p <- ggplot() + 
  geom_bar(DBA.1.DSC.proportion, mapping = aes(x = 1, y = Freq, fill = Var1), stat = "identity") + 
  scale_fill_manual(values = decidual.color[3:7]) + 
  coord_polar("y", start = 0) + 
  theme_void()

load("20210827-figures/files/DBA.1-immune-20221215.RData")
DimPlot(DBA.1.immune.m, label = T)
DBA.1.immune.m.proportion <- data.frame(table(DBA.1.immune.m@active.ident), check.names = F)
DBA.1.immune.m.proportion$Var1 <- factor(DBA.1.immune.m.proportion$Var1, 
                                         levels = c("Monocyte", "Mac", "Proliferating Mac", "DC-1", "Neutrophil", "NK", "T cell", "iDSC_DBA1"))
Immune.p <- ggplot() + 
  geom_bar(DBA.1.immune.m.proportion, mapping = aes(x = 1, y = Freq, fill = Var1), stat = "identity") + 
  scale_fill_manual(values = immune.color) + 
  coord_polar("y", start = 0) + 
  theme_void()

load("20220622-Comments/05.DBA/DBA-EC-clusters-20221215.RData")
DimPlot(DBA.1.EC, label = T)
DBA.1.EC.proportion <- data.frame(table(DBA.1.EC@active.ident), check.names = F)
EC.p <- ggplot() + 
  geom_bar(DBA.1.EC.proportion, mapping = aes(x = 1, y = Freq, fill = Var1), stat = "identity") + 
  scale_fill_manual(values = EC.color) + 
  coord_polar("y", start = 0) + 
  theme_void()
pdf("20220622-Comments/05.DBA/00.figures/31.DBA-DSC-Immune-EC-subclusters-proportion-pie-plot.pdf", width = 10, height = 10)
DSC.p + Immune.p + EC.p
dev.off()

#Correlation between Normal and DBA DSC subclusters 
load("20210827-figures/files/DBA.1-DSC-20221212.Rdata")
load("20210827-figures/files/UE-DSC-web.RData")
UE.decidual.susbet$cellcluster <- UE.decidual.susbet$sub_clusters
DBA.Normal.DSC <- merge(DBA.1.DSC, UE.decidual.susbet)
DBA.Normal.DSC@active.ident <- factor(DBA.Normal.DSC$cellcluster)
DBA.Normal.DSC.avg <- AverageExpression(DBA.Normal.DSC, assays = "RNA")
DBA.Normal.DSC.avg <- DBA.Normal.DSC.avg$RNA
DBA.Normal.DSC.avg.cor <- cor(DBA.Normal.DSC.avg)
DBA.Normal.DSC.avg.cor <- DBA.Normal.DSC.avg.cor[c(3,5,8,9,11), c(4,6,7,10,12)]
pheatmap(DBA.Normal.DSC.avg.cor, scale = "column", cluster_rows = T, cluster_cols = T)

DBA.Normal.DSC.avg.D3 <- (DBA.Normal.DSC.avg$`D3_Proliferating DSC` + 1)/(DBA.Normal.DSC.avg$D3_Top2a + 1)
DBA.Normal.DSC.avg.D4 <- (DBA.Normal.DSC.avg$`D4_Biomineral-regulatory DSC` + 1)/(DBA.Normal.DSC.avg$D4_Ptn + 1)
DBA.Normal.DSC.avg.D5 <- (DBA.Normal.DSC.avg$`D5_Nourishing DSC` + 1)/(DBA.Normal.DSC.avg$D5_Gatm + 1)
DBA.Normal.DSC.avg.D6 <- (DBA.Normal.DSC.avg$`D6_Angiogenesis DSC` + 1)/(DBA.Normal.DSC.avg$D6_S100a8 + 1)
DBA.Normal.DSC.avg.D7 <- (DBA.Normal.DSC.avg$`D7_Postmature DSC` + 1)/(DBA.Normal.DSC.avg$D7_Prl8a2 + 1)
D3.cor.p <- ggplot() + 
  geom_point(DBA.Normal.DSC.avg, mapping = aes(x = log2(`D3_Proliferating DSC`+1), y = log2(D3_Top2a+1)), size = 1, alpha = 0.2, color = "darkblue") + 
  theme_classic() +
  labs(title = "D3 (R2=0.95)") +
  ylim(0, 10) + 
  xlim(0, 10)
D4.cor.p <- ggplot() + 
  geom_point(DBA.Normal.DSC.avg, mapping = aes(x = log2(`D4_Biomineral-regulatory DSC`+1), y = log2(D4_Ptn+1)), size = 1, alpha = 0.2, color = "darkblue") + 
  theme_classic() +
  labs(title = "D4 (R2=0.96)") +
  ylim(0, 10) + 
  xlim(0, 10)
D5.cor.p <- ggplot() + 
  geom_point(DBA.Normal.DSC.avg, mapping = aes(x = log2(`D5_Nourishing DSC`+1), y = log2(D5_Gatm+1)), size = 1, alpha = 0.2, color = "darkblue") + 
  theme_classic() +
  labs(title = "D5 (R2=0.96)") +
  ylim(0, 10) + 
  xlim(0, 10)
D6.cor.p <- ggplot() + 
  geom_point(DBA.Normal.DSC.avg, mapping = aes(x = log2(`D6_Angiogenesis DSC`+1), y = log2(D6_S100a8+1)), size = 1, alpha = 0.2, color = "darkblue")+ 
  theme_classic() +
  labs(title = "D6 (R2=0.95)") +
  ylim(0, 10) + 
  xlim(0, 10)
D7.cor.p <- ggplot() + 
  geom_point(DBA.Normal.DSC.avg, mapping = aes(x = log2(`D7_Postmature DSC`+1), y = log2(D7_Prl8a2+1)), size = 1, alpha = 0.2, color = "darkblue")+ 
  theme_classic() +
  labs(title = "D7 (R2=0.73)") +
  ylim(0, 10) + 
  xlim(0, 10)
D3.lm <- lm(log2(`D3_Proliferating DSC`+1) ~ log2(D3_Top2a+1), data = DBA.Normal.DSC.avg)
summary(D3.lm) # R^2=0.95
D4.lm <- lm(log2(`D4_Biomineral-regulatory DSC`+1) ~ log2(D4_Ptn+1), data = DBA.Normal.DSC.avg)
summary(D4.lm) # R^2=0.96
D5.lm <- lm(log2(`D5_Nourishing DSC`+1) ~ log2(D5_Gatm+1), data = DBA.Normal.DSC.avg)
summary(D5.lm) # R^2=0.96
D6.lm <- lm(log2(`D6_Angiogenesis DSC`+1) ~ log2(D6_S100a8+1), data = DBA.Normal.DSC.avg)
summary(D6.lm) #R^2=0.95
D7.lm <- lm(log2(`D7_Postmature DSC`+1) ~ log2(D7_Prl8a2+1), data = DBA.Normal.DSC.avg)
summary(D7.lm) #R^2=0.73
pdf("20220622-Comments/05.DBA/00.figures/33.Normal-DBA-DSC-subclusters-correlation-dotplot.pdf", width = 15, height = 10)
D3.cor.p + D4.cor.p + D5.cor.p + D6.cor.p + D7.cor.p
dev.off()

#CBA Inflammation EC gene expression
load("20220622-Comments/05.DBA/DBA-EC-clusters-20221215.RData")
DBA.1.EC <- RenameIdents(DBA.1.EC, "EPC" = "Inflammatory EC")
DBA.1.EC$cellcluster <- DBA.1.EC@active.ident
DimPlot(DBA.1.EC)
genes.1 <- c("Acta1", "Tagln", "Cnn1", "Pecam1", "Thbd", "Nos3")
FeaturePlot(DBA.1.EC, features = genes.1)
load("20220622-Comments/Normal-all-subclusters.RData")
load("20220622-Comments/05.DBA/DBA.1-merge-embryo-SCT.RData")
DBA.1.ST.inte$cellcluster <- gsub("Proliterating EC", "Proliferating EC", DBA.1.ST.inte$cellcluster)
DBA.1.ST.inte$sample <- "CBA"
DBA.all$sample <- "CBA"
Normal.DBA.all <- merge(Normal.all, DBA.1.ST.inte)
Normal.DBA.all <- NormalizeData(Normal.DBA.all)
Normal.DBA.all.EC.s <- subset(Normal.DBA.all, cellcluster %in% c("Proliferating EC", "Angiogenic EC", "Arterial EC", "Venous EC1", "Venous EC2", "EPC"))
Normal.DBA.all.immune.s <- subset(Normal.DBA.all, cellcluster %in% c("DC-1", "DC-2", "NK", "Mac", "T cell", "Neutrophil", "B", "Monocyte"))
Normal.DBA.all.iDSC.s <- subset(Normal.DBA.all, cellcluster %in% c("iDSC0", "iDSC1", "iDSC2"))

VlnPlot(Normal.DBA.all.EC.s, features = genes.1, split.by = "sample", pt.size = 0)

genes.2 <- c("Vcam1", "Icam1", "Sele", "Selp", "Fgl2")
FeaturePlot(DBA.1.EC, features = "Fgl2")
VlnPlot(Normal.DBA.all.EC.s, features = "Fgl2", split.by = "sample", pt.size = 0)

VlnPlot(Normal.DBA.all.EC.s, features = c("Hif1a", "Hif3a", "Epas1"), split.by = "sample", pt.size = 0)
VlnPlot(Normal.DBA.all.immune.s, features = c("Hif1a", "Hif3a", "Epas1"), split.by = "sample", pt.size = 0)
VlnPlot(Normal.DBA.all.iDSC.s, features = c("Hif1a", "Hif3a", "Epas1"), split.by = "sample", pt.size = 0)
VlnPlot(Normal.DBA.all.iDSC.s, features = "Fgl2", split.by = "sample", pt.size = 0)

#Cell proportion in Normal and CBA
Normal.all$time <- factor(Normal.all$time, levels = c("T55", "T65", "T75", "T85", "T95", "T105"))
#normal DSC subclusters
Normal.DSC <- subset(Normal.all, cellcluster %in% c("D1_eSF", "D2_Pre-DSC", "D3_Top2a", "D4_Ptn", "D5_Gatm", "D6_S100a8", "D7_Prl8a2"))

Normal.DSC.time.tb <- data.frame(table(Normal.DSC$time, Normal.DSC$cellcluster))
Normal.DSC.time <- ggplot() +
  geom_bar(Normal.DSC.time.tb, mapping = aes(x = Var1, y = Freq, fill = Var2), stat = "identity", position = "fill") +
  scale_fill_manual(values = decidual.color) + 
  theme_classic() + 
  ggtitle(label = "DSC subcluster proportion over time")

Normal.DSC.all.tb <- data.frame(table(Normal.DSC$cellcluster))
Normal.DSC.all <- ggplot() +
  geom_bar(Normal.DSC.all.tb, mapping = aes(x = 1, y = Freq, fill = Var1), stat = "identity", position = "fill") +
  scale_fill_manual(values = decidual.color) + 
  theme_classic() + 
  ggtitle(label = "DSC subcluster proportion in total")

Normal.DSC.E85.E95.tb <- Normal.DSC.time.tb[Normal.DSC.time.tb$Var1 %in% c("T85", "T95"), ]
Normal.DSC.E85.E95 <- ggplot() +
  geom_bar(Normal.DSC.E85.E95.tb, mapping = aes(x = 1, y = Freq, fill = Var2), stat = "identity", position = "fill") +
  scale_fill_manual(values = decidual.color) + 
  theme_classic() + 
  ggtitle(label = "DSC subcluster proportion in E8.5/E9.5")

#Normal EC
Normal.EC <- subset(Normal.all, cellcluster %in% c("Angiogenic EC","Proliferating EC","Arterial EC","Venous EC2","Venous EC1","EC-5"))
Normal.EC$cellcluster <- factor(Normal.EC$cellcluster, levels = c("Angiogenic EC","Proliferating EC","Venous EC1","Venous EC2","Arterial EC","EC-5"))
Normal.EC.time.tb <- data.frame(table(Normal.EC$time, Normal.EC$cellcluster))
Normal.EC.time <- ggplot() +
  geom_bar(Normal.EC.time.tb, mapping = aes(x = Var1, y = Freq, fill = Var2), stat = "identity", position = "fill") +
  scale_fill_manual(values = EC.color[c(1,3,2,4,5,6)]) + 
  theme_classic() + 
  ggtitle(label = "EC subcluster proportion over time")
Normal.EC.all.tb <- data.frame(table(Normal.EC$cellcluster))
Normal.EC.all <- ggplot() +
  geom_bar(Normal.EC.all.tb, mapping = aes(x = 1, y = Freq, fill = Var1), stat = "identity", position = "fill") +
  scale_fill_manual(values = EC.color[c(1,3,2,4,5,6)]) + 
  theme_classic() +
  ggtitle(label = "EC subcluster proportion in all")
Normal.EC.E85.E95.tb <- Normal.EC.time.tb[Normal.EC.time.tb$Var1 %in% c("T85", "T95"), ]
Normal.EC.E85.E95 <- ggplot() +
  geom_bar(Normal.EC.E85.E95.tb, mapping = aes(x = 1, y = Freq, fill = Var2), stat = "identity", position = "fill") +
  scale_fill_manual(values = EC.color[c(1,3,2,4,5,6)]) + 
  theme_classic() + 
  ggtitle(label = "EC subcluster proportion in E8.5/E9.5")
#Normal immune
Normal.immune <- subset(Normal.all, cellcluster %in% immune.levels[1:9])
Normal.immune$cellcluster <- factor(Normal.immune$cellcluster, levels = immune.levels)
Normal.immune.time.tb <- data.frame(table(Normal.immune$time, Normal.immune$cellcluster))
Normal.immune.time <- ggplot() +
  geom_bar(Normal.immune.time.tb, mapping = aes(x = Var1, y = Freq, fill = Var2), stat = "identity", position = "fill") +
  scale_fill_manual(values = immune.color) + 
  theme_classic() + 
  ggtitle(label = "immune subcluster proportion over time")
Normal.immune.all.tb <- data.frame(table(Normal.immune$cellcluster))
Normal.immune.all <- ggplot() +
  geom_bar(Normal.immune.all.tb, mapping = aes(x = 1, y = Freq, fill = Var1), stat = "identity", position = "fill") +
  scale_fill_manual(values = immune.color) + 
  theme_classic() + 
  ggtitle(label = "immune subcluster proportion in all")
Normal.immune.E85.E95.tb <- Normal.immune.time.tb[Normal.immune.time.tb$Var1 %in% c("T85", "T95"), ]
Normal.immune.E85.E95 <- ggplot() +
  geom_bar(Normal.immune.E85.E95.tb, mapping = aes(x = 1, y = Freq, fill = Var2), stat = "identity", position = "fill") +
  scale_fill_manual(values = immune.color) + 
  theme_classic() + 
  ggtitle(label = "immune subcluster proportion at E8.5/E9.5")

#CBA DSC subclusters
CBA.DSC <- subset(DBA.1.ST.inte, cellcluster %in% c("D1_eSF", "D2_Pre-DSC", "D3_Top2a", "D4_Ptn", "D5_Gatm", "D6_S100a8", "D7_Prl8a2"))
CBA.DSC <- subset(CBA.DSC, time =="DBA.1")
CBA.DSC.all.tb <- data.frame(table(CBA.DSC$cellcluster))
CBA.DSC.all <- ggplot() +
  geom_bar(CBA.DSC.all.tb, mapping = aes(x = 1, y = Freq, fill = Var1), stat = "identity", position = "fill") +
  scale_fill_manual(values = decidual.color[3:7]) + 
  theme_classic() + 
  ggtitle(label = "CBA_DSC subcluster proportion")
#CBA EC
CBA.EC <- subset(DBA.1.ST.inte, cellcluster %in% c("Angiogenic EC","Proliferating EC","Arterial EC","Venous EC2","Venous EC1","EPC"))
CBA.EC$cellcluster <- factor(CBA.EC$cellcluster, levels = c("Angiogenic EC","Proliferating EC","Venous EC1","Venous EC2","Arterial EC","EPC"))
CBA.EC <- subset(CBA.EC, time == "DBA.1")
CBA.EC.all.tb <- data.frame(table(CBA.EC$cellcluster))
CBA.EC.all <- ggplot() +
  geom_bar(CBA.EC.all.tb, mapping = aes(x = 1, y = Freq, fill = Var1), stat = "identity", position = "fill") +
  scale_fill_manual(values = c(EC.color[c(1,3,2,4,5)], "#00FA9A")) + 
  theme_classic() +
  ggtitle(label = "CBA_EC subcluster proportion")
#CBA immune
CBA.immune <- subset(DBA.1.ST.inte, cellcluster %in% c("Monocyte", "Mac", "Proliferating Mac", 'DC-1', "Neutrophil", "NK", "T cell"))
CBA.immune$cellcluster <- factor(CBA.immune$cellcluster, levels = c("Monocyte", "Mac", "Proliferating Mac", 'DC-1', "Neutrophil", "NK", "T cell"))
CBA.immune.all.tb <- data.frame(table(CBA.immune$cellcluster))
CBA.immune.all <- ggplot() +
  geom_bar(CBA.immune.all.tb, mapping = aes(x = 1, y = Freq, fill = Var1), stat = "identity", position = "fill") +
  scale_fill_manual(values = immune.color[c(1:4,6,7,9)]) + 
  theme_classic() + 
  ggtitle(label = "CBA_immune subcluster proportion")
pdf("20220622-Comments/05.DBA/00.figures/60.Normal_CBA_DSC_EC_immune_subcluster_proportion-barplot.pdf", width = 20, height = 20)
Normal.DSC.time + Normal.DSC.all + Normal.DSC.E85.E95 + 
  Normal.EC.time + Normal.EC.all + Normal.EC.E85.E95 + 
  Normal.immune.time + Normal.immune.all + Normal.immune.E85.E95 + 
  CBA.DSC.all + CBA.EC.all + CBA.immune.all
dev.off()





