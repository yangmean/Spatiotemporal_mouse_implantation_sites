library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(patchwork)
color <- brewer.pal(11, "RdBu")
color <- color[c(2:5, 7:10)]
dec.vitro <- read.table("../../20220818-decidualization-in-vitro/result/Merge_count.csv", sep = ",", header = T,row.names = 1, check.names = F)
dec.vitro <- dec.vitro[, c(3,4,5,6,9,8)]
write.table(dec.vitro, file = "/Data2/Users/yangmin/Jennie/20230219_data_upload/04.files/DSC_in_vitro_decidulaization.csv", sep = ",", quote = F, row.names = T, col.names = T)
dec.vitro.nor <- data.frame(t(apply(dec.vitro, 1, function(x){x/colSums(dec.vitro)* 1000000})), check.names = F)
write.table(dec.vitro.nor, file = "/Data2/Users/yangmin/Jennie/20230219_data_upload/04.files/DSC_in_vitro_decidulaization-nor.csv", sep = ",", quote = F, row.names = T, col.names = T)
dec.vitro.nor <- dec.vitro.nor[rowSums(dec.vitro.nor)>6, ]
#dec.vitro.nor <- dec.vitro.nor[, c(2,1,3,7,4,5,6,9,8)]
dec.vitro.nor.fc <- (dec.vitro.nor[, 2:6] + 5) / (dec.vitro.nor[, 1] + 5)
dec.vitro.nor.fc2 <- dec.vitro.nor[apply(dec.vitro.nor.fc, 1, max)>2, ]
dec.vitro.nor.fc2.p <- pheatmap(dec.vitro.nor.fc2, scale = "row", cluster_rows = F, cluster_cols = F, show_rownames = F, 
         color = colorRampPalette(colors = rev(color))(100), main = "Upregulated (FC>2, N = 1066)")
dec.vitro.nor.fc.5 <- dec.vitro.nor[apply(dec.vitro.nor.fc, 1, min)<0.5, ]
dec.vitro.nor.fc.5.p <- pheatmap(dec.vitro.nor.fc.5, scale = "row", cluster_rows = F, cluster_cols = F, show_rownames = F, 
         color = colorRampPalette(colors = rev(color))(100), main = "Downregulated (FC<0.5, N = 865)")
save(dec.vitro.nor.fc2, file = "20220622-Comments/03.vitro-decidualization/Dec-vitro-nor-fc2.RData")
save(dec.vitro.nor.fc.5, file = "20220622-Comments/03.vitro-decidualization/Dec-vitro-nor-fc05.RData")

pdf("20220622-Comments/03.vitro-decidualization/04.Up-Down-regulated-heatmap.pdf", width = 10, height = 10)
as.ggplot(dec.vitro.nor.fc2.p) + as.ggplot(dec.vitro.nor.fc.5.p)
dev.off()
#load scRNA-seq
load("../20210201-integrate/UE-decidual-harmony-pc30-subset-rename.RData")
UE.decidual.susbet$cellname <- gsub("Postn-D", "Lum-D", UE.decidual.susbet$cellname)
UE.decidual.susbet$cellname <- gsub("Ifit1-D", "S100a8-D", UE.decidual.susbet$cellname)
UE.decidual.susbet$cellname <- gsub("Hist1h2ap-D", "Top2a-D", UE.decidual.susbet$cellname)
UE.decidual.susbet@active.ident <- factor(UE.decidual.susbet$cellname)
UE.decidual.susbet <- AddModuleScore(UE.decidual.susbet, features = list(rownames(dec.vitro.nor.fc2)), name = "vitro.up.reg")
UE.decidual.susbet <- AddModuleScore(UE.decidual.susbet, features = list(rownames(dec.vitro.nor.fc.5)), name = "vitro.down.reg")
UE.vitro.up.p <- FeaturePlot(UE.decidual.susbet, features = c("vitro.up.reg1"), min.cutoff = 0.085, pt.size = 0.5) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
UE.vitro.down.p <- FeaturePlot(UE.decidual.susbet, features = c("vitro.down.reg1"), pt.size = 0.5, min.cutoff = 0.04, max.cutoff = 0.15) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
pdf("20220622-Comments/03.vitro-decidualization/05.scRNA_seq-vitro-up-down-regulated-featureplot.pdf", width = 10, height = 5)
UE.vitro.up.p + UE.vitro.down.p + plot_layout(ncol = 2)
dev.off()
#load spatial RNA-seq
load("UE-decidual-MET-ST-monolce3-20211105.RData")
decidual.MET.mono3.psdeudotime <- read.table("UE-decidual-MET-ST-monolce3-pseudotime-df-20211105.csv", header = T, check.names = F, sep = ",")
load("20210821-decidual-pseudotime/Decidual-ST-downsample-1000-seurat.RData")
Decidual.bin50.subsample <- AddModuleScore(Decidual.bin50.subsample, assay = "SCT", features = list(rownames(dec.vitro.nor.fc2)), name = "vitro.up.reg")
Decidual.bin50.subsample <- AddModuleScore(Decidual.bin50.subsample, assay = "SCT", features = list(rownames(dec.vitro.nor.fc.5)), name = "vitro.down.reg")
Decidual.ST.vitro.score.df <- data.frame(cellname = colnames(Decidual.bin50.subsample), check.names = F,
                                         vitro.up = Decidual.bin50.subsample$vitro.up.reg1, 
                                         vitro.down = Decidual.bin50.subsample$vitro.down.reg1)
decidual.MET.mono3.psdeudotime.vitro.score <- merge(decidual.MET.mono3.psdeudotime, Decidual.ST.vitro.score.df, by = "cellname", all = F)
decidual.MET.mono3.psdeudotime.vitro.score[decidual.MET.mono3.psdeudotime.vitro.score$vitro.up<0.09, ]$vitro.up <- 0.09
decidual.MET.mono3.psdeudotime.vitro.score[decidual.MET.mono3.psdeudotime.vitro.score$vitro.down<0.02, ]$vitro.down <- 0.01
ST.umap.up.p <- ggplot() + 
  geom_point(decidual.MET.mono3.psdeudotime.vitro.score, mapping = aes(x = UMAP_1, y = UMAP_2, color = vitro.up), size = 0.5) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + 
  theme_classic()
ST.umap.down.p <- ggplot() + 
  geom_point(decidual.MET.mono3.psdeudotime.vitro.score, mapping = aes(x = UMAP_1, y = UMAP_2, color = vitro.down), size = 0.5) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + 
  theme_classic()
pdf("20220622-Comments/03.vitro-decidualization/06.ST-vitro-up-down-regulated-featureplot.pdf", width = 10, height = 5)
ST.umap.up.p + ST.umap.down.p + plot_layout(ncol = 2)
dev.off()
load("20210818-spatial/B6-bin50-spatial.RData")
B6.bin50 <- AddModuleScore(B6.bin50, assay = "SCT", features = list(rownames(dec.vitro.nor.fc2)), name = "vitro.up.reg")
B6.bin50 <- AddModuleScore(B6.bin50, assay = "SCT", features = list(rownames(dec.vitro.nor.fc.5)), name = "vitro.down.reg")
B6.bin50.vitro.up.p <- SpatialFeaturePlot(B6.bin50, features = "vitro.up.reg1", ncol = 1, min.cutoff = 0.06, stroke = NA, pt.size.factor = 2.25) +
  scale_fill_gradientn(colours = brewer.pal(9, "PuBu")) +
  NoLegend()
B6.bin50.vitro.down.p <- SpatialFeaturePlot(B6.bin50, features = "vitro.down.reg1", ncol = 1, min.cutoff = 0.01, stroke = NA, pt.size.factor = 2.25) +
  scale_fill_gradientn(colours = brewer.pal(9, "PuBu")) + 
  NoLegend()
pdf("20220622-Comments/03.vitro-decidualization/07.E85-B6-ST-vitro-up-down-regulated-featureplot.pdf")
B6.bin50.vitro.up.p + B6.bin50.vitro.down.p + plot_layout(ncol = 2)
dev.off()
#E65.2
load("20210818-spatial/E65-2-bin50-spatial-phago-map-remove.RData")
E65.2.bin50 <- AddModuleScore(E65.2.bin50, assay = "SCT", features = list(rownames(dec.vitro.nor.fc2)), name = "vitro.up.reg")
E65.2.bin50 <- AddModuleScore(E65.2.bin50, assay = "SCT", features = list(rownames(dec.vitro.nor.fc.5)), name = "vitro.down.reg")
E65.2.bin50.vitro.up.p <- SpatialFeaturePlot(E65.2.bin50, features = "vitro.up.reg1", ncol = 1, min.cutoff = 0.07, stroke = NA, pt.size.factor = 2.6) +
  scale_fill_gradientn(colours = brewer.pal(9, "PuBu")) +
  NoLegend()
E65.2.bin50.vitro.down.p <- SpatialFeaturePlot(E65.2.bin50, features = "vitro.down.reg1", ncol = 1, min.cutoff = 0.01, stroke = NA, pt.size.factor = 2.6) +
  scale_fill_gradientn(colours = brewer.pal(9, "PuBu")) + 
  NoLegend()
pdf("20220622-Comments/03.vitro-decidualization/07.E65-2-ST-vitro-up-down-regulated-featureplot.pdf")
E65.2.bin50.vitro.up.p + E65.2.bin50.vitro.down.p + plot_layout(ncol = 2)
dev.off()
#E75.1
load("20210818-spatial/E75-1-bin50-spatial-phago-map-remove.RData")
E75.1.bin50 <- AddModuleScore(E75.1.bin50, assay = "SCT", features = list(rownames(dec.vitro.nor.fc2)), name = "vitro.up.reg")
E75.1.bin50 <- AddModuleScore(E75.1.bin50, assay = "SCT", features = list(rownames(dec.vitro.nor.fc.5)), name = "vitro.down.reg")
E75.1.bin50.vitro.up.p <- SpatialFeaturePlot(E75.1.bin50, features = "vitro.up.reg1", ncol = 1, min.cutoff = 0.09, stroke = NA, pt.size.factor = 2) +
  scale_fill_gradientn(colours = brewer.pal(9, "PuBu")) +
  NoLegend()
E75.1.bin50.vitro.down.p <- SpatialFeaturePlot(E75.1.bin50, features = "vitro.down.reg1", ncol = 1, min.cutoff = 0.01, stroke = NA, pt.size.factor = 2) +
  scale_fill_gradientn(colours = brewer.pal(9, "PuBu")) + 
  NoLegend()
pdf("20220622-Comments/03.vitro-decidualization/07.E75-2-ST-vitro-up-down-regulated-featureplot.pdf")
E75.1.bin50.vitro.up.p + E75.1.bin50.vitro.down.p + plot_layout(ncol = 2)
dev.off()
load("20210818-spatial/E95-bin50-spatial.RData")
E95.bin50 <- AddModuleScore(E95.bin50, assay = "SCT", features = list(rownames(dec.vitro.nor.fc2)), name = "vitro.up.reg")
E95.bin50 <- AddModuleScore(E95.bin50, assay = "SCT", features = list(rownames(dec.vitro.nor.fc.5)), name = "vitro.down.reg")
E95.bin50.vitro.up.p <- SpatialFeaturePlot(E95.bin50, features = "vitro.up.reg1", ncol = 1, min.cutoff = 0.07, stroke = NA, pt.size.factor = 1.5) +
  scale_fill_gradientn(colours = brewer.pal(9, "PuBu")) +
  NoLegend()
E95.bin50.vitro.down.p <- SpatialFeaturePlot(E95.bin50, features = "vitro.down.reg1", ncol = 1, min.cutoff = 0.01, stroke = NA, pt.size.factor = 1.5) +
  scale_fill_gradientn(colours = brewer.pal(9, "PuBu")) + 
  NoLegend()
pdf("20220622-Comments/03.vitro-decidualization/07.E95-ST-vitro-up-down-regulated-featureplot.pdf")
E95.bin50.vitro.up.p + E95.bin50.vitro.down.p + plot_layout(ncol = 2)
dev.off()


decidua.DEGs.top20 <- read.table("20210827-figures/files/UE-harmony-decidual-markers-top20.csv", sep = ",", header = T, check.names = F)
decidua.DEGs <- read.table("20210827-figures/files/UE-harmony-decidual-markers.csv", sep = ",", header = T, check.names = F)
#eSF
dec.vitro.nor.eSF.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Lum-D", ]$gene), ]
dec.vitro.nor.eSF.top20 <- na.omit(dec.vitro.nor.eSF.top20)
pheatmap(dec.vitro.nor.eSF.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.eSF <- dec.vitro.nor[as.character(decidua.DEGs[decidua.DEGs$cluster=="Lum-D", ]$gene), ]
dec.vitro.nor.eSF <- na.omit(dec.vitro.nor.eSF)
dec.vitro.nor.eSF.h <- pheatmap(dec.vitro.nor.eSF, scale = "row", cluster_rows = T, cluster_cols = F, 
         show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.eSF.h.c <- cutree(dec.vitro.nor.eSF.h$tree_row, k = 2)
dec.vitro.nor.eSF.h.d <- names(dec.vitro.nor.eSF.h.c[dec.vitro.nor.eSF.h.c==1])
dec.vitro.nor.eSF.d <- dec.vitro.nor.eSF[dec.vitro.nor.eSF.h.d, ]
dec.vitro.nor.eSF.d.p <- pheatmap(dec.vitro.nor.eSF.d, scale = "row", cluster_rows = F, show_rownames = F,
                                  cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                  main = "D1_eSF DEG expression (87/142)")
#D2
dec.vitro.nor.D2.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Sfrp4-D", ]$gene), ]
dec.vitro.nor.D2.top20 <- na.omit(dec.vitro.nor.D2.top20)
pheatmap(dec.vitro.nor.D2.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D2 <- dec.vitro.nor[as.character(decidua.DEGs[decidua.DEGs$cluster=="Sfrp4-D", ]$gene), ]
dec.vitro.nor.D2 <- na.omit(dec.vitro.nor.D2)
dec.vitro.nor.D2.h <- pheatmap(dec.vitro.nor.D2, scale = "row", cluster_rows = T, cluster_cols = F, 
                                show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.D2.h.c <- cutree(dec.vitro.nor.D2.h$tree_row, k = 2)
dec.vitro.nor.D2.h.d <- names(dec.vitro.nor.D2.h.c[dec.vitro.nor.D2.h.c==1])
dec.vitro.nor.D2.d <- dec.vitro.nor.D2[dec.vitro.nor.D2.h.d, ]
dec.vitro.nor.D2.d.p <- pheatmap(dec.vitro.nor.D2.d, scale = "row", cluster_rows = F, show_rownames = F,
                                  cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                  main = "D2_Pre-DSC DEG expression (56/91)")
#D3
dec.vitro.nor.D3.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Top2a-D", ]$gene), ]
dec.vitro.nor.D3.top20 <- na.omit(dec.vitro.nor.D3.top20)
pheatmap(dec.vitro.nor.D3.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D3 <- dec.vitro.nor[as.character(decidua.DEGs[decidua.DEGs$cluster=="Top2a-D", ]$gene), ]
dec.vitro.nor.D3 <- na.omit(dec.vitro.nor.D3)
dec.vitro.nor.D3.h <- pheatmap(dec.vitro.nor.D3, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.D3.h.c <- cutree(dec.vitro.nor.D3.h$tree_row, k = 2)
dec.vitro.nor.D3.h.d <- names(dec.vitro.nor.D3.h.c[dec.vitro.nor.D3.h.c==1])
dec.vitro.nor.D3.d <- dec.vitro.nor.D3[dec.vitro.nor.D3.h.d, ]
dec.vitro.nor.D3.d.p <- pheatmap(dec.vitro.nor.D3.d, scale = "row", cluster_rows = F, show_rownames = F,
                                  cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                  main = "D3_proliferation-DSC DEG expression (58/85)")
#D4
dec.vitro.nor.D4.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Ptn-D", ]$gene), ]
dec.vitro.nor.D4.top20 <- na.omit(dec.vitro.nor.D4.top20)
pheatmap(dec.vitro.nor.D4.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D4 <- dec.vitro.nor[as.character(decidua.DEGs[decidua.DEGs$cluster=="Ptn-D", ]$gene), ]
dec.vitro.nor.D4 <- na.omit(dec.vitro.nor.D4)
dec.vitro.nor.D4.h <- pheatmap(dec.vitro.nor.D4, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.D4.h.c <- cutree(dec.vitro.nor.D4.h$tree_row, k = 2)
dec.vitro.nor.D4.h.d <- names(dec.vitro.nor.D4.h.c[dec.vitro.nor.D4.h.c==2])
dec.vitro.nor.D4.d <- dec.vitro.nor.D4[dec.vitro.nor.D4.h.d, ]
dec.vitro.nor.D4.d.p <- pheatmap(dec.vitro.nor.D4.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "D4_DSC DEG expression (11/21)")
#D5
dec.vitro.nor.D5.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Gatm-D", ]$gene), ]
dec.vitro.nor.D5.top20 <- na.omit(dec.vitro.nor.D5.top20)
pheatmap(dec.vitro.nor.D5.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D5 <- dec.vitro.nor[as.character(decidua.DEGs[decidua.DEGs$cluster=="Gatm-D", ]$gene), ]
dec.vitro.nor.D5 <- na.omit(dec.vitro.nor.D5)
dec.vitro.nor.D5.h <- pheatmap(dec.vitro.nor.D5, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.D5.h.c <- cutree(dec.vitro.nor.D5.h$tree_row, k = 2)
dec.vitro.nor.D5.h.d <- names(dec.vitro.nor.D5.h.c[dec.vitro.nor.D5.h.c==2])
dec.vitro.nor.D5.d <- dec.vitro.nor.D5[dec.vitro.nor.D5.h.d, ]
dec.vitro.nor.D5.d.p <- pheatmap(dec.vitro.nor.D5.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "D5_Nourishing-DSC DEG expression (11/19)")

#D6
dec.vitro.nor.D6.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="S100a8-D", ]$gene), ]
dec.vitro.nor.D6.top20 <- na.omit(dec.vitro.nor.D6.top20)
pheatmap(dec.vitro.nor.D6.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D6 <- dec.vitro.nor[as.character(decidua.DEGs[decidua.DEGs$cluster=="S100a8-D", ]$gene), ]
dec.vitro.nor.D6 <- na.omit(dec.vitro.nor.D6)
dec.vitro.nor.D6.h <- pheatmap(dec.vitro.nor.D6, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.D6.h.c <- cutree(dec.vitro.nor.D6.h$tree_row, k = 2)
dec.vitro.nor.D6.h.d <- names(dec.vitro.nor.D6.h.c[dec.vitro.nor.D6.h.c==1])
dec.vitro.nor.D6.d <- dec.vitro.nor.D6[dec.vitro.nor.D6.h.d, ]
dec.vitro.nor.D6.d.p <- pheatmap(dec.vitro.nor.D6.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "D6_Angiogenic-DSC DEG expression (163/238)")
#D7
dec.vitro.nor.D7.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Prl8a2-D", ]$gene), ]
dec.vitro.nor.D7.top20 <- na.omit(dec.vitro.nor.D7.top20)
pheatmap(dec.vitro.nor.D7.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D7 <- dec.vitro.nor[as.character(decidua.DEGs[decidua.DEGs$cluster=="Prl8a2-D", ]$gene), ]
dec.vitro.nor.D7 <- na.omit(dec.vitro.nor.D7)
dec.vitro.nor.D7.h <- pheatmap(dec.vitro.nor.D7, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.D7.h.c <- cutree(dec.vitro.nor.D7.h$tree_row, k = 2)
dec.vitro.nor.D7.h.d <- names(dec.vitro.nor.D7.h.c[dec.vitro.nor.D7.h.c==1])
dec.vitro.nor.D7.d <- dec.vitro.nor.D7[dec.vitro.nor.D7.h.d, ]
dec.vitro.nor.D7.d.p <- pheatmap(dec.vitro.nor.D7.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "D7_Postmature-DSC DEG expression (42/77)")
pdf("20220622-Comments/03.vitro-decidualization/01.DSC-subcluster-DEG-expression-vitro-RNA-seq-heatmap.pdf", width = 10, height = 20)
as.ggplot(dec.vitro.nor.eSF.d.p) + as.ggplot(dec.vitro.nor.D2.d.p) + as.ggplot(dec.vitro.nor.D3.d.p) + as.ggplot(dec.vitro.nor.D4.d.p) + 
  as.ggplot(dec.vitro.nor.D5.d.p) + as.ggplot(dec.vitro.nor.D6.d.p)+ as.ggplot(dec.vitro.nor.D7.d.p) + plot_layout(ncol = 2)
dev.off()
#scatterplot
#D1:D1-FBS
dec.vitro.nor.D1 <- dec.vitro.nor[, c(1, 2)]
dec.vitro.nor.D1$color <- 0
dec.vitro.nor.D1[(dec.vitro.nor.D1$`B2-D1` + 5)/(dec.vitro.nor.D1$`B2-2-percent-FBS-D1` + 5)>2, ]$color <- 1
dec.vitro.nor.D1[(dec.vitro.nor.D1$`B2-D1` + 5)/(dec.vitro.nor.D1$`B2-2-percent-FBS-D1` + 5)<0.5, ]$color <- -1
dec.vitro.nor.D1$color <- factor(dec.vitro.nor.D1$color, levels = c(-1, 0, 1))
dec.vitro.nor.D1.p <- ggplot() + 
  geom_point(dec.vitro.nor.D1, mapping = aes(x = log(dec.vitro.nor.D1$`B2-2-percent-FBS-D1` + 1), 
                                          y = log(dec.vitro.nor.D1$`B2-D1` + 1), color = color), size = 1) + 
  scale_color_manual(values = c("darkblue", "lightgrey", "darkred")) +
  theme(panel.background = element_rect(color = "black", fill = "white"), 
        panel.grid = element_blank()) + 
  labs(title = "B2_D1 vs B2_FBS_D1")
#D2:D1-FBS
dec.vitro.nor.D2 <- dec.vitro.nor[, c(1, 3)]
dec.vitro.nor.D2$color <- 0
dec.vitro.nor.D2[(dec.vitro.nor.D2$`B2-D2` + 5)/(dec.vitro.nor.D2$`B2-2-percent-FBS-D1` + 5)>2, ]$color <- 1
dec.vitro.nor.D2[(dec.vitro.nor.D2$`B2-D2` + 5)/(dec.vitro.nor.D2$`B2-2-percent-FBS-D1` + 5)<0.5, ]$color <- -1
dec.vitro.nor.D2$color <- factor(dec.vitro.nor.D2$color, levels = c(-1, 0, 1))
dec.vitro.nor.D2.p <- ggplot() + 
  geom_point(dec.vitro.nor.D2, mapping = aes(x = log(dec.vitro.nor.D2$`B2-2-percent-FBS-D1` + 1), 
                                             y = log(dec.vitro.nor.D2$`B2-D2` + 1), color = color), size = 1) + 
  scale_color_manual(values = c("darkblue", "lightgrey", "darkred")) +
  theme(panel.background = element_rect(color = "black", fill = "white"), 
        panel.grid = element_blank()) + 
  labs(title = "B2_D2 vs B2_FBS_D1")
#D3:D1-FBS
dec.vitro.nor.D3 <- dec.vitro.nor[, c(1, 4)]
dec.vitro.nor.D3$color <- 0
dec.vitro.nor.D3[(dec.vitro.nor.D3$`B2-D3` + 5)/(dec.vitro.nor.D3$`B2-2-percent-FBS-D1` + 5)>2, ]$color <- 1
dec.vitro.nor.D3[(dec.vitro.nor.D3$`B2-D3` + 5)/(dec.vitro.nor.D3$`B2-2-percent-FBS-D1` + 5)<0.5, ]$color <- -1
dec.vitro.nor.D3$color <- factor(dec.vitro.nor.D3$color, levels = c(-1, 0, 1))
dec.vitro.nor.D3.p <- ggplot() + 
  geom_point(dec.vitro.nor.D3, mapping = aes(x = log(dec.vitro.nor.D3$`B2-2-percent-FBS-D1` + 1), 
                                             y = log(dec.vitro.nor.D3$`B2-D3` + 1), color = color), size = 1) + 
  scale_color_manual(values = c("darkblue", "lightgrey", "darkred")) +
  theme(panel.background = element_rect(color = "black", fill = "white"), 
        panel.grid = element_blank()) + 
  labs(title = "B2_D3 vs B2_FBS_D1")
#D4:D1-FBS
dec.vitro.nor.D4 <- dec.vitro.nor[, c(1, 5)]
dec.vitro.nor.D4$color <- 0
dec.vitro.nor.D4[(dec.vitro.nor.D4$`B2-D4` + 5)/(dec.vitro.nor.D4$`B2-2-percent-FBS-D1` + 5)>2, ]$color <- 1
dec.vitro.nor.D4[(dec.vitro.nor.D4$`B2-D4` + 5)/(dec.vitro.nor.D4$`B2-2-percent-FBS-D1` + 5)<0.5, ]$color <- -1
dec.vitro.nor.D4$color <- factor(dec.vitro.nor.D4$color, levels = c(-1, 0, 1))
dec.vitro.nor.D4.p <- ggplot() + 
  geom_point(dec.vitro.nor.D4, mapping = aes(x = log(dec.vitro.nor.D4$`B2-2-percent-FBS-D1` + 1), 
                                             y = log(dec.vitro.nor.D4$`B2-D4` + 1), color = color), size = 1) + 
  scale_color_manual(values = c("darkblue", "lightgrey", "darkred")) +
  theme(panel.background = element_rect(color = "black", fill = "white"), 
        panel.grid = element_blank()) + 
  labs(title = "B2_D4 vs B2_FBS_D1")
#D5:D1-FBS
dec.vitro.nor.D5 <- dec.vitro.nor[, c(1, 6)]
dec.vitro.nor.D5$color <- 0
dec.vitro.nor.D5[(dec.vitro.nor.D5$`B2-D5` + 5)/(dec.vitro.nor.D5$`B2-2-percent-FBS-D1` + 5)>2, ]$color <- 1
dec.vitro.nor.D5[(dec.vitro.nor.D5$`B2-D5` + 5)/(dec.vitro.nor.D5$`B2-2-percent-FBS-D1` + 5)<0.5, ]$color <- -1
dec.vitro.nor.D5$color <- factor(dec.vitro.nor.D5$color, levels = c(-1, 0, 1))
highlight.down.genes <- c("Mmp14", "Lgals1", "Igfbp7", "Col1a1", "Col14a1", "Col5a1", "Sfrp1", "Egr1", "Sparcl1", "Col6a3", 
                     "Col5a2", "Junb", "Col12a1", "Col6a1", "Ogn", "Eln", "Col1a2", "Col3a1")
highlight.up.genes <- c("Lgals3", "Fabp4", "Cstb", "Pdgfra", 'Htra1', "Htra3", "S100a6", "Prl8a2", "Mt1", 
                        "Procr", "Lamc1", "Lamb2", "Anxa8", "Hmgn5", "Hand2")
dec.vitro.nor.D5.annotate <- dec.vitro.nor.D5[c(highlight.down.genes, highlight.up.genes),]
dec.vitro.nor.D5.p.lable <- ggplot() + 
  geom_point(dec.vitro.nor.D5, mapping = aes(x = log(dec.vitro.nor.D5$`B2-2-percent-FBS-D1` + 1), 
                                             y = log(dec.vitro.nor.D5$`B2-D5` + 1), color = color), size = 1) + 
  scale_color_manual(values = c("darkblue", "lightgrey", "darkred")) +
  theme(panel.background = element_rect(color = "black", fill = "white"), 
        panel.grid = element_blank()) + 
  labs(title = "B2_D5 vs B2_FBS_D1 (label)") + 
  geom_point(dec.vitro.nor.D5.annotate, mapping = aes(x = log(dec.vitro.nor.D5.annotate$`B2-2-percent-FBS-D1` + 1), 
                                                      y = log(dec.vitro.nor.D5.annotate$`B2-D5` + 1)), color = "orange", size = 2) + 
  annotate(geom = "text", label = rownames(dec.vitro.nor.D5.annotate), color = "black",
           x = log(dec.vitro.nor.D5.annotate$`B2-2-percent-FBS-D1` + 1),
           y = log(dec.vitro.nor.D5.annotate$`B2-D5` + 1))
pdf("20220622-Comments/03.vitro-decidualization/02.Day5-vs-Day1-scatter-plot-label.pdf", width = 10, height = 10)
dec.vitro.nor.D5.p.lable
dev.off()

dec.vitro.nor.D5.p <- ggplot() + 
  geom_point(dec.vitro.nor.D5, mapping = aes(x = log(dec.vitro.nor.D5$`B2-2-percent-FBS-D1` + 1), 
                                             y = log(dec.vitro.nor.D5$`B2-D5` + 1), color = color), size = 1) + 
  scale_color_manual(values = c("darkblue", "lightgrey", "darkred")) +
  theme(panel.background = element_rect(color = "black", fill = "white"), 
        panel.grid = element_blank()) + 
  labs(title = "B2_D5 vs B2_FBS_D1")
pdf("20220622-Comments/03.vitro-decidualization/02.Day5-vs-Day1-scatter-plot.pdf", width = 6, height = 5.5)
dec.vitro.nor.D5.p
dev.off()

#calculation DSC subcluster DEG enrichment score in in-vitro-decidualization samples
dec.vitro.nor.min_max <- data.frame(t(apply(dec.vitro.nor, 1, function(x){(x-min(x))/(max(x)-min(x))})), check.names = F)
dec.vitro.nor.min_max$gene <- rownames(dec.vitro.nor.min_max)
decidua.DEGs.RNA.expr.min_max <- merge(dec.vitro.nor.min_max, decidua.DEGs, by = "gene", all = F)
decidua.DEGs.RNA.expr.min_max <- decidua.DEGs.RNA.expr.min_max[, c(1:7, 13)]
decidua.DEGs.RNA.expr.min_max.m <- melt(decidua.DEGs.RNA.expr.min_max, id.vars = c("gene", "cluster"))
pdf("20220622-Comments/03.vitro-decidualization/03.DSC-DEGs-enrich-score-vitro-decidualization-sample-boxplot.pdf", width = 20, height = 10)
ggplot() +
  geom_boxplot(decidua.DEGs.RNA.expr.min_max.m, mapping = aes(x = variable, y = value, fill = variable)) + 
  facet_wrap(~ cluster) + 
  theme(panel.background = element_rect(color = "black", fill = "white"), 
        panel.grid = element_blank())
dev.off()


#Postn-D expression
dec.vitro.nor.postn <- dec.vitro.nor[c("Igfbp6", "Postn", "Sparcl1", "Dusp1", "Penk", "Ctla2a", "Calml4", "Atf3", "Krtdap", "Lox", "Dcn", 
                                       "Nfkbia", "Klf4", "H2-D1", "Car2", "Ifi27l2a", "Rprm", "Fgl2", "Hspa1a", "Hmgcs1", "Col5a1", "Igfbp5"), ]
pheatmap(dec.vitro.nor.postn, scale = "row", cluster_rows = T, cluster_cols = F)


#20230101 使用单细胞中每个DSC亚群的差异基因看在in vitro decidualization中的表达变化。min.ptc=0.25
UE.MET <- subset(UE.decidual, subset = cellname %in% c("Lum-D", "Sfrp4-D", "Gatm-D", "Prl8a2-D", "S100a8-D", "Top2a-D"))
UE.MET.markers <- FindAllMarkers(UE.MET, min.pct = 0.25, logfc.threshold = 0.5, only.pos = F)
UE.MET.markers.top10 <- UE.MET.markers %>% group_by(cluster) %>% top_n(n = 10, avg_logFC)
DoHeatmap(UE.MET, features = UE.MET.markers.top10$gene)
#eSF
dec.vitro.nor.eSF.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Lum-D", ]$gene), ]
dec.vitro.nor.eSF.top20 <- na.omit(dec.vitro.nor.eSF.top20)
pheatmap(dec.vitro.nor.eSF.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.eSF <- dec.vitro.nor[as.character(UE.MET.markers[UE.MET.markers$cluster=="Lum-D", ]$gene), ]
dec.vitro.nor.eSF <- na.omit(dec.vitro.nor.eSF)
dec.vitro.nor.eSF.h <- pheatmap(dec.vitro.nor.eSF, scale = "row", cluster_rows = T, cluster_cols = F, 
                                show_rownames = F, cutree_rows = 7, clustering_method = "ward.D")
dec.vitro.nor.eSF.h.c <- cutree(dec.vitro.nor.eSF.h$tree_row, k = 7)
dec.vitro.nor.eSF.h.d <- names(dec.vitro.nor.eSF.h.c[dec.vitro.nor.eSF.h.c %in% c(1,2,4,7)])
dec.vitro.nor.eSF.d <- dec.vitro.nor.eSF[dec.vitro.nor.eSF.h.d, ]
dec.vitro.nor.eSF.d.p <- pheatmap(dec.vitro.nor.eSF.d, scale = "row", cluster_rows = F, show_rownames = F,
                                  cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                  main = "D1_eSF DEG expression (281/449)")
#D2
dec.vitro.nor.D2.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Sfrp4-D", ]$gene), ]
dec.vitro.nor.D2.top20 <- na.omit(dec.vitro.nor.D2.top20)
pheatmap(dec.vitro.nor.D2.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D2 <- dec.vitro.nor[as.character(UE.MET.markers[UE.MET.markers$cluster=="Sfrp4-D", ]$gene), ]
dec.vitro.nor.D2 <- na.omit(dec.vitro.nor.D2)
dec.vitro.nor.D2.h <- pheatmap(dec.vitro.nor.D2, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 5, clustering_method = "ward.D")
dec.vitro.nor.D2.h.c <- cutree(dec.vitro.nor.D2.h$tree_row, k = 5)
dec.vitro.nor.D2.h.d <- names(dec.vitro.nor.D2.h.c[dec.vitro.nor.D2.h.c %in% c(2,3,4)])
dec.vitro.nor.D2.d <- dec.vitro.nor.D2[dec.vitro.nor.D2.h.d, ]
dec.vitro.nor.D2.d.p <- pheatmap(dec.vitro.nor.D2.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "D2_Pre-DSC DEG expression (197/344)")
#D3
dec.vitro.nor.D3.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Top2a-D", ]$gene), ]
dec.vitro.nor.D3.top20 <- na.omit(dec.vitro.nor.D3.top20)
pheatmap(dec.vitro.nor.D3.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D3 <- dec.vitro.nor[as.character(UE.MET.markers[UE.MET.markers$cluster=="Top2a-D", ]$gene), ]
dec.vitro.nor.D3 <- na.omit(dec.vitro.nor.D3)
dec.vitro.nor.D3.h <- pheatmap(dec.vitro.nor.D3, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 6, clustering_method = "ward.D")
dec.vitro.nor.D3.h.c <- cutree(dec.vitro.nor.D3.h$tree_row, k = 6)
dec.vitro.nor.D3.h.d <- names(dec.vitro.nor.D3.h.c[dec.vitro.nor.D3.h.c %in% c(1,2,3)])
dec.vitro.nor.D3.d <- dec.vitro.nor.D3[dec.vitro.nor.D3.h.d, ]
dec.vitro.nor.D3.d.p <- pheatmap(dec.vitro.nor.D3.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "D3_proliferation-DSC DEG expression (110/197)")
#D5
dec.vitro.nor.D5.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Gatm-D", ]$gene), ]
dec.vitro.nor.D5.top20 <- na.omit(dec.vitro.nor.D5.top20)
pheatmap(dec.vitro.nor.D5.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D5 <- dec.vitro.nor[as.character(UE.MET.markers[UE.MET.markers$cluster=="Gatm-D", ]$gene), ]
dec.vitro.nor.D5 <- na.omit(dec.vitro.nor.D5)
dec.vitro.nor.D5.h <- pheatmap(dec.vitro.nor.D5, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 4, clustering_method = "ward.D")
dec.vitro.nor.D5.h.c <- cutree(dec.vitro.nor.D5.h$tree_row, k = 4)
dec.vitro.nor.D5.h.d <- names(dec.vitro.nor.D5.h.c[dec.vitro.nor.D5.h.c %in% c(3,4)])
dec.vitro.nor.D5.d <- dec.vitro.nor.D5[dec.vitro.nor.D5.h.d, ]
dec.vitro.nor.D5.d.p <- pheatmap(dec.vitro.nor.D5.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "D5_Nourishing-DSC DEG expression (72/122)")

#D6
dec.vitro.nor.D6.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="S100a8-D", ]$gene), ]
dec.vitro.nor.D6.top20 <- na.omit(dec.vitro.nor.D6.top20)
pheatmap(dec.vitro.nor.D6.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D6 <- dec.vitro.nor[as.character(UE.MET.markers[UE.MET.markers$cluster=="S100a8-D", ]$gene), ]
dec.vitro.nor.D6 <- na.omit(dec.vitro.nor.D6)
dec.vitro.nor.D6.h <- pheatmap(dec.vitro.nor.D6, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 7, clustering_method = "ward.D")
dec.vitro.nor.D6.h.c <- cutree(dec.vitro.nor.D6.h$tree_row, k = 7)
dec.vitro.nor.D6.h.d <- names(dec.vitro.nor.D6.h.c[dec.vitro.nor.D6.h.c %in% c(1,3,4,7)])
dec.vitro.nor.D6.d <- dec.vitro.nor.D6[dec.vitro.nor.D6.h.d, ]
dec.vitro.nor.D6.d.p <- pheatmap(dec.vitro.nor.D6.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "D6_Angiogenic-DSC DEG expression (265/471)")
#D7
dec.vitro.nor.D7.top20 <- dec.vitro.nor[as.character(decidua.DEGs.top20[decidua.DEGs.top20$cluster=="Prl8a2-D", ]$gene), ]
dec.vitro.nor.D7.top20 <- na.omit(dec.vitro.nor.D7.top20)
pheatmap(dec.vitro.nor.D7.top20, scale = "row", cluster_rows = T, cluster_cols = F)
dec.vitro.nor.D7 <- dec.vitro.nor[as.character(UE.MET.markers[UE.MET.markers$cluster=="Prl8a2-D", ]$gene), ]
dec.vitro.nor.D7 <- na.omit(dec.vitro.nor.D7)
dec.vitro.nor.D7.h <- pheatmap(dec.vitro.nor.D7, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.D7.h.c <- cutree(dec.vitro.nor.D7.h$tree_row, k = 2)
dec.vitro.nor.D7.h.d <- names(dec.vitro.nor.D7.h.c[dec.vitro.nor.D7.h.c==2])
dec.vitro.nor.D7.d <- dec.vitro.nor.D7[dec.vitro.nor.D7.h.d, ]
dec.vitro.nor.D7.d.p <- pheatmap(dec.vitro.nor.D7.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "D7_Postmature-DSC DEG expression (234/381)")

pdf("20220622-Comments/03.vitro-decidualization/01.DSC-subcluster-DEG-expression-vitro-RNA-seq-heatmap-20230101.pdf", width = 10, height = 20)
as.ggplot(dec.vitro.nor.eSF.d.p) + as.ggplot(dec.vitro.nor.D2.d.p) + as.ggplot(dec.vitro.nor.D3.d.p) + 
  as.ggplot(dec.vitro.nor.D5.d.p) + as.ggplot(dec.vitro.nor.D6.d.p)+ as.ggplot(dec.vitro.nor.D7.d.p) + plot_layout(ncol = 3)
dev.off()

load("../20210201-integrate/UE-decidual-harmony-pc30-rename.RData")
UE.decidual <- subset(UE.decidual, subset = cellname != "Prg-")
UE.decidual <- RenameIdents(UE.decidual, "Postn-D" = "Lum-D")
UE.decidual <- RenameIdents(UE.decidual, "Ifit1-D" = "S100a8-D")
UE.decidual <- RenameIdents(UE.decidual, "Hist1h2ap-D" = "Top2a-D")

UE.MET <- subset(UE.decidual, subset = cellname %in% c("Lum-D", "Sfrp4-D", "Gatm-D", "Prl8a2-D", "S100a8-D", "Top2a-D"))
UE.MET.markers <- FindAllMarkers(UE.MET, min.pct = 0.25, logfc.threshold = 0.5, only.pos = F)
UE.MET.change <- UE.MET.avg[as.character(unique(UE.MET.markers$gene)), 1:6]
UE.MET.change <- UE.MET.change[, c(3,4,1,5,2,6)]
UE.MET.avg <- AverageExpression(UE.MET, slot = "data", assays = "RNA")
UE.MET.avg <- UE.MET.avg$RNA
UE.MET.avg <- UE.MET.avg[rowSums(UE.MET.avg)>0, ]
UE.MET.avg <- UE.MET.avg * 100

#使用figure4J的8个pattern
UE.MET.change.heatmap.cluster <- read.table("20210827-figures/UE-MET-DEGs-8-clusters.csv", sep = ",", header = T)#this file code in 20210827-figures.Rmd file
UE.MET.change.heatmap.cluster$cluster <- gsub(8, 1, UE.MET.change.heatmap.cluster$cluster)
UE.MET.change.heatmap.cluster$cluster <- gsub(6, 2, UE.MET.change.heatmap.cluster$cluster)
UE.MET.change.heatmap.cluster$cluster <- gsub(5, 4, UE.MET.change.heatmap.cluster$cluster)
min_max <- function(x){(x-min(x))/(max(x)-min(x))}

tmp.df <- UE.MET.change[as.character(UE.MET.change.heatmap.cluster[UE.MET.change.heatmap.cluster$cluster %in% c(1,8), ]$gene), ]
pheatmap(tmp.df[, c(3,4,1,5,2,6)], scale = "row", cluster_rows = F, cluster_cols = F)


UE.MET.change.boxplot.list <- list()
for (i in c(1,2,3,4,7)) {
  tmp.df <- UE.MET.change[as.character(UE.MET.change.heatmap.cluster[UE.MET.change.heatmap.cluster$cluster %in% c(i), ]$gene), ]
  tmp.df <- tmp.df[rowSums(tmp.df)>0, ]
  tmp.df.scale <- data.frame(t(apply(tmp.df, 1, scale)), check.names = F)
  colnames(tmp.df.scale) <- c("Lum-D", "Sfrp4-D", "Top2a-D", "Gatm-D", "S100a8-D", "Prl8a2-D")
  tmp.df.scale <- melt(tmp.df.scale)
  UE.MET.change.boxplot.list[[i]] <- ggplot() + 
    geom_boxplot(tmp.df.scale, mapping = aes(x = variable, y = value, fill = variable)) +
    scale_fill_manual(values = decidual.color[c(1,2,3,5,6,7)]) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", color = "black"))
}
M.genes <- c("Fn1", "Vim", "Col3a1", "Col1a1", "Col1a2", "Zeb1", "Igfbp5", "Sfrp4", "Igf1", "Tcf4", "Cdh11", "Fscn1", "Wnt5a", "F2r")
E.genes <- c("Jun", "Bmp2", "Lama5", "Angpt4", "Icam1", "Adam19", "Trem2", "Cdh5", "Gja1")
intersect(UE.MET.change.heatmap.cluster[UE.MET.change.heatmap.cluster$cluster == 2, ]$gene, M.genes)

pdf("20210827-figures/figure2/Figure2-13-MET-average-comparing-boxplot-cluster-20230103.pdf", width = 15, height = 15)
cowplot::plot_grid(UE.MET.change.boxplot.list[[1]], UE.MET.change.boxplot.list[[2]], UE.MET.change.boxplot.list[[3]], 
                   UE.MET.change.boxplot.list[[4]], UE.MET.change.boxplot.list[[7]], ncol=3)
dev.off()
#GO
UE.MET.change.heatmap.cluster.GO <- data.frame()
for (c in c(1,2,3,4,7)) {
  print(c)
  geneset <- UE.MET.change.heatmap.cluster[UE.MET.change.heatmap.cluster$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    gene = enrich.GO$geneID)
  UE.MET.change.heatmap.cluster.GO <- rbind(UE.MET.change.heatmap.cluster.GO, tmp)
}
#write.table(UE.MET.change.heatmap.cluster.GO, file = "")

dec.vitro.nor.P1 <- dec.vitro.nor[as.character(UE.MET.change.heatmap.cluster[UE.MET.change.heatmap.cluster$cluster %in% c(1,8), ]$gene), ]
dec.vitro.nor.P1 <- na.omit(dec.vitro.nor.P1)
dec.vitro.nor.P1.h <- pheatmap(dec.vitro.nor.P1, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 6, clustering_method = "ward.D")
dec.vitro.nor.P1.h.c <- cutree(dec.vitro.nor.P1.h$tree_row, k = 6)
dec.vitro.nor.P1.h.d <- names(dec.vitro.nor.P1.h.c[dec.vitro.nor.P1.h.c %in% c(1,2,3)])
dec.vitro.nor.P1.d <- dec.vitro.nor.P1[dec.vitro.nor.P1.h.d, ]
dec.vitro.nor.P1.d.p <- pheatmap(dec.vitro.nor.P1.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "P1_Postmature-DSC DEG expression (92/183)")

dec.vitro.nor.P2 <- dec.vitro.nor[as.character(UE.MET.change.heatmap.cluster[UE.MET.change.heatmap.cluster$cluster %in% c(7), ]$gene), ]
dec.vitro.nor.P2 <- na.omit(dec.vitro.nor.P2)
dec.vitro.nor.P2.h <- pheatmap(dec.vitro.nor.P2, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 4, clustering_method = "ward.D")
dec.vitro.nor.P2.h.c <- cutree(dec.vitro.nor.P2.h$tree_row, k = 4)
dec.vitro.nor.P2.h.d <- names(dec.vitro.nor.P2.h.c[dec.vitro.nor.P2.h.c %in% c(2,3)])
dec.vitro.nor.P2.d <- dec.vitro.nor.P2[dec.vitro.nor.P2.h.d, ]
dec.vitro.nor.P2.d.p <- pheatmap(dec.vitro.nor.P2.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "P2_Angiogenic-DSC DEG expression (133/203)")

dec.vitro.nor.P3 <- dec.vitro.nor[as.character(UE.MET.change.heatmap.cluster[UE.MET.change.heatmap.cluster$cluster %in% c(3), ]$gene), ]
dec.vitro.nor.P3 <- na.omit(dec.vitro.nor.P3)
dec.vitro.nor.P3.h <- pheatmap(dec.vitro.nor.P3, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.P3.h.c <- cutree(dec.vitro.nor.P3.h$tree_row, k = 2)
dec.vitro.nor.P3.h.d <- names(dec.vitro.nor.P3.h.c[dec.vitro.nor.P3.h.c %in% c(1)])
dec.vitro.nor.P3.d <- dec.vitro.nor.P3[dec.vitro.nor.P3.h.d, ]
dec.vitro.nor.P3.d.p <- pheatmap(dec.vitro.nor.P3.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "P3_Proliferating-DSC DEG expression (66/123)")

dec.vitro.nor.P4 <- dec.vitro.nor[as.character(UE.MET.change.heatmap.cluster[UE.MET.change.heatmap.cluster$cluster %in% c(2,6), ]$gene), ]
dec.vitro.nor.P4 <- na.omit(dec.vitro.nor.P4)
dec.vitro.nor.P4.h <- pheatmap(dec.vitro.nor.P4, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 2, clustering_method = "ward.D")
dec.vitro.nor.P4.h.c <- cutree(dec.vitro.nor.P4.h$tree_row, k = 2)
dec.vitro.nor.P4.h.d <- names(dec.vitro.nor.P4.h.c[dec.vitro.nor.P4.h.c %in% c(2)])
dec.vitro.nor.P4.d <- dec.vitro.nor.P4[dec.vitro.nor.P4.h.d, ]
dec.vitro.nor.P4.d.p <- pheatmap(dec.vitro.nor.P4.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "P1-eSF_P2-Pre-DSC_ DEG expression (157/246)")

dec.vitro.nor.P5 <- dec.vitro.nor[as.character(UE.MET.change.heatmap.cluster[UE.MET.change.heatmap.cluster$cluster %in% c(4,5), ]$gene), ]
dec.vitro.nor.P5 <- na.omit(dec.vitro.nor.P5)
dec.vitro.nor.P5.h <- pheatmap(dec.vitro.nor.P5, scale = "row", cluster_rows = T, cluster_cols = F, 
                               show_rownames = F, cutree_rows = 5, clustering_method = "ward.D")
dec.vitro.nor.P5.h.c <- cutree(dec.vitro.nor.P5.h$tree_row, k = 5)
dec.vitro.nor.P5.h.d <- names(dec.vitro.nor.P5.h.c[dec.vitro.nor.P5.h.c %in% c(1,2)])
dec.vitro.nor.P5.d <- dec.vitro.nor.P5[dec.vitro.nor.P5.h.d, ]
dec.vitro.nor.P5.d.p <- pheatmap(dec.vitro.nor.P5.d, scale = "row", cluster_rows = F, show_rownames = F,
                                 cluster_cols = F, color = colorRampPalette(colors = rev(color))(100), border_color = NA, 
                                 main = "P5_Proliferating-DSC DEG expression (90/224)")
pdf("20220622-Comments/03.vitro-decidualization/01.DSC-subcluster-DEG-expression-vitro-RNA-seq-heatmap-20230101.pdf", width = 10, height = 20)
as.ggplot(dec.vitro.nor.P1.d.p) + as.ggplot(dec.vitro.nor.P2.d.p) + as.ggplot(dec.vitro.nor.P3.d.p) + 
  as.ggplot(dec.vitro.nor.P4.d.p) + as.ggplot(dec.vitro.nor.P5.d.p) + plot_layout(ncol = 3)
dev.off()


