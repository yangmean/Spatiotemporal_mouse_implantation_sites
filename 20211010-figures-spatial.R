library(ggplot2)
library(Seurat)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(patchwork)
library(ggpubr)
#color, levels
source("../../../01.scripts/StackedVlnPlot.R")
major.color <- c("#5F9EA0", "#D2691E", "#20B2AA", "#6B8E23", "#006400", "#6495ED", "#8B0000", "#4B0082", "#CD5C5C", "#DAA520", "#9370DB", "#8B4513", "#FFD700")
decidual.color <- c("#6A5ACD","#9ACD32","#FB8072","#8DD3C7","#08519C","#FDB462","#228B22","#FCCDE5","#D9D9D9")
immune.color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#2E8B57", "#CAB2D6", "#FF7F00", "#FDBF6F", "#6A3D9A")
EC.color <- c("#7B68EE", "#1B9E77", "#D95F02", "#00BFFF", "#E7298A", "#CCEBC5")
tropho.color <- c("#3CB371", "#DA70D6", "#4682B4", "#483D8B", "#008B8B", "#FF8C00")
spatial_clusters <- c(brewer.pal(9, "Greys"))
spatial_clusters <- colorRampPalette(spatial_clusters[4:9])(17)
spatial_clusters <- sample(spatial_clusters, 17)
major.levels <- c("Decidual", "Immune", "EC", "SMC", "Epithelial", "TP", "Erythroid", "Mesenchymal stem cell", "Mesenchymal-epithelial transition", "Blood P", "Viseral endoderm", "Allantois mesodermal", "Fetal EC")
decidual.levels <- c("Lum-D", "Sfrp4-D", "Top2a-D", "Ptn-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "cluster9-1", "Prg-")
immune.levels <- c("Mono", "Mac", "Pro-Mac", "DC-1", "DC-2", "Neutrophil", "NK", "B", "T", "FB-immune")
EC.levels <- c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "EC-5")
tropho.levels <- c("tropho-0", "tropho-1", "tropho-2", "tropho-3", "tropho-4", "tropho-5")

annotation.color <- list(cell_type = c("Lum-D"="#6A5ACD", "Sfrp4-D"="#9ACD32", "Top2a-D"="#FB8072","Ptn-D"="#8DD3C7","Gatm-D"="#08519C",
                                       "S100a8-D"="#FDB462","Prl8a2-D"="#228B22", "cluster9-1"="#FCCDE5", "Mono"="#A6CEE3",
                                       "Mac"="#1F78B4","Pro-Mac"="#B2DF8A","DC-1"="#33A02C", "DC-2" = "#FB9A99", "NK"="#CAB2D6",
                                       "B"="#FF7F00","T"="#FDBF6F", "Neutrophil" = "#2E8B57", "FB-immune"="#6A3D9A", "EC-0"="#7B68EE", 
                                       "EC-1"="#1B9E77","EC-2"="#D95F02","EC-3"="#00BFFF", "EC-4" = "#E7298A", "EC-5" = "#CCEBC5",
                                       "Epithelial"="#006400","Erythroid"="#8B0000", "SMC"="#6B8E23",
                                       "Mesenchymal stem cell"="#4B0082","Mesenchymal-epithelial transition"="#CD5C5C",
                                       "Blood P"="#DAA520", "Allantois mesodermal"="#8B4513","Viseral endoderm"="#9370DB", 
                                       "Fetal EC" = "#FFD700", "tropho-0"="#3CB371","tropho-1"="#DA70D6", 
                                       "tropho-2"="#4682B4","tropho-3"="#483D8B","tropho-4"="#008B8B", "tropho-5" = "#FF8C00"), 
                         time = c("E65.2.bin50" = "#D6604D", "E65.1.bin50" = "#D6604D", "E75.2.bin50" = "#F4A582", "E75.1.bin50" = "#F4A582", 
                                  "E85.B6.bin50" = "#92C5DE", "E85.C1.bin50" = "#92C5DE", "E95.bin50" = "#4393C3"))
##############################################################################################################################
#E65.1
load("20210818-spatial/E65-1-bin50-spatial-phago-map-remove.RData")
E65.1.bin50.celltype.df <- read.table("20210818-spatial/E65-1-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E65.1.bin50.celltype.df <- E65.1.bin50.celltype.df[E65.1.bin50.celltype.df$cell_type != "Erythroid", ]
E65.1.cell_type.levels <- c(decidual.levels[1:7], immune.levels[c(1:7,9,10)], EC.levels, major.levels[4:5])
E65.1.color <- c(decidual.color[1:7], immune.color[c(1:7,9,10)], EC.color, major.color[4:5])
E65.1.bin50.celltype.df$cell_type <- factor(E65.1.bin50.celltype.df$cell_type, levels = E65.1.cell_type.levels)
pdf("20210827-figures/figure2-spatial/Figure2-01-E65-1-bin50-all-types.pdf", width = 5.2, height = 3.5)
ggplot() + 
  geom_point(E65.1.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = E65.1.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
#subtype
plots <- list()
for (c in seq(1, 24, 1)) {
  print(E65.1.cell_type.levels[c])
  c.df <- E65.1.bin50.celltype.df[E65.1.bin50.celltype.df$cell_type==E65.1.cell_type.levels[c], ]
  plot <- ggplot() + 
    geom_point(E65.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.1) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = E65.1.color[c], size = 0.3) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(E65.1.cell_type.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E65-1-bin50-subtype.pdf", width = 12, height = 11)
wrap_plots(plots, ncol = 6)
dev.off()

#E65.2
load("20210818-spatial/E65-2-bin50-spatial-phago-map-remove.RData")
E65.2.bin50.celltype.df <- read.table("20210818-spatial/E65-2-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E65.2.bin50.celltype.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type != "Erythroid", ]
E65.2.cell_type.levels <- c(decidual.levels[1:7], immune.levels[c(1:10)], EC.levels, major.levels[4:5])
E65.2.color <- c(decidual.color[1:7], immune.color[c(1:10)], EC.color, major.color[4:5])
E65.2.bin50.celltype.df$cell_type <- factor(E65.2.bin50.celltype.df$cell_type, levels = E65.2.cell_type.levels)
pdf("20210827-figures/figure2-spatial/Figure2-01-E65-2-bin50-all-types.pdf", width = 5.5, height = 3.5)
ggplot() + 
  geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = E65.2.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
#subtype
plots <- list()
for (c in seq(1, 25, 1)) {
  print(E65.2.cell_type.levels[c])
  c.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type==E65.2.cell_type.levels[c], ]
  plot <- ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.1) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = E65.2.color[c], size = 0.3) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(E65.2.cell_type.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E65-2-bin50-subtype.pdf", width = 12, height = 11)
wrap_plots(plots, ncol = 6)
dev.off()

#E75.1
load("20210818-spatial/E75-1-bin50-spatial-phago-map-remove.RData")
E75.1.bin50.celltype.df <- read.table("20210818-spatial/E75-1-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E75.1.bin50.celltype.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type != "Erythroid", ]
E75.1.cell_type.levels <- c(decidual.levels[c(1:3,5:7)], immune.levels, EC.levels[c(1,3:5)], major.levels[4:5])
E75.1.color <- c(decidual.color[c(1:3,5:7)], immune.color, EC.color[c(1,3:5)], major.color[4:5])
E75.1.bin50.celltype.df$cell_type <- factor(E75.1.bin50.celltype.df$cell_type, levels = E75.1.cell_type.levels)
pdf("20210827-figures/figure2-spatial/Figure2-01-E75-1-bin50-all-types.pdf", width = 6.3, height = 4.8)
ggplot() + 
  geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = E75.1.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#subtype
plots <- list()
for (c in seq(1, 22, 1)) {
  print(E75.1.cell_type.levels[c])
  c.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type==E75.1.cell_type.levels[c], ]
  plot <- ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.1) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = E75.1.color[c], size = 0.3) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(E75.1.cell_type.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E75-1-bin50-subtype.pdf", width = 13, height = 12)
wrap_plots(plots, ncol = 6)
dev.off()

#E75.2
load("20210818-spatial/E75-2-bin50-spatial-phago-map-remove.RData")
E75.2.bin50.celltype.df <- read.table("20210818-spatial/E75-2-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E75.2.bin50.celltype.df <- E75.2.bin50.celltype.df[E75.2.bin50.celltype.df$cell_type != "Erythroid", ]
E75.2.cell_type.levels <- c(decidual.levels[1:7], immune.levels[c(1:4,6:10)], EC.levels[1:5], major.levels[4:5], tropho.levels[5])
E75.2.color <- c(decidual.color[1:7], immune.color[c(1:4,6:10)], EC.color[1:5], major.color[4:5], tropho.color[5])
E75.2.bin50.celltype.df$cell_type <- factor(E75.2.bin50.celltype.df$cell_type, levels = E75.2.cell_type.levels)
pdf("20210827-figures/figure2-spatial/Figure2-01-E75-2-bin50-all-types.pdf", width = 7.3, height = 5.3)
ggplot() + 
  geom_point(E75.2.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = E75.2.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#subtype
plots <- list()
for (c in seq(1, 24, 1)) {
  print(E75.2.cell_type.levels[c])
  c.df <- E75.2.bin50.celltype.df[E75.2.bin50.celltype.df$cell_type==E75.2.cell_type.levels[c], ]
  plot <- ggplot() + 
    geom_point(E75.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.1) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = E75.2.color[c], size = 0.3) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(E75.2.cell_type.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E75-2-bin50-subtype.pdf", width = 13.5, height = 12)
wrap_plots(plots, ncol = 6)
dev.off()

#E85-B6
load("20210818-spatial/E85-B6-bin50-spatial-phago-map-remove.RData")
B6.bin50.celltype.df <- read.table("20210818-spatial/E85-B6-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
B6.bin50.celltype.df <- B6.bin50.celltype.df[B6.bin50.celltype.df$cell_type != "Erythroid", ]
B6.cell_type.levels <- c(decidual.levels[1:7], immune.levels[c(1:4,7:10)], EC.levels[1:4], major.levels[c(4:5,11)], tropho.levels[c(3:5)], "embryo")
B6.bin50.celltype.df$cell_type <- factor(B6.bin50.celltype.df$cell_type, levels = B6.cell_type.levels)
B6.color <- c(decidual.color[1:7], immune.color[c(1:4,7:10)], EC.color[1:4], major.color[c(4:5,11)], tropho.color[c(3:5)], "lightgray")
pdf("20210827-figures/figure2-spatial/Figure2-01-E85-B6-bin50-all-types.pdf", width = 6.3, height = 4.5)
ggplot() + 
  geom_point(B6.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = B6.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#subtype
plots <- list()
for (c in seq(1, 26, 1)) {
  print(B6.cell_type.levels[c])
  c.df <- B6.bin50.celltype.df[B6.bin50.celltype.df$cell_type==B6.cell_type.levels[c], ]
  plot <- ggplot() + 
    geom_point(B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.1) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = B6.color[c], size = 0.3) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(B6.cell_type.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E85-B6-bin50-subtype.pdf", width = 12, height = 12)
wrap_plots(plots, ncol = 6)
dev.off()

#E85-C1
load("20210818-spatial/E85-C1-bin50-spatial-phago-map-remove.RData")
C1.bin50.celltype.df <- read.table("20210818-spatial/E85-C1-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
C1.bin50.celltype.df <- C1.bin50.celltype.df[C1.bin50.celltype.df$cell_type != "Erythroid", ]

C1.cell_type.levels <- c(decidual.levels[c(1:3,5:7)], immune.levels[c(1:4,7:10)], EC.levels[2:4], major.levels[c(4:5)],tropho.levels[c(3,5)])
C1.bin50.celltype.df$cell_type <- factor(C1.bin50.celltype.df$cell_type, levels = C1.cell_type.levels)
C1.color <- c(decidual.color[c(1:3,5:7)], immune.color[c(1:4,7:10)], EC.color[2:4], major.color[c(4:5)], tropho.color[c(3,5)])
pdf("20210827-figures/figure2-spatial/Figure2-01-E85-C1-bin50-all-types.pdf", width = 5.8, height = 4.2)
ggplot() + 
  geom_point(C1.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = C1.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
#subtype
plots <- list()
for (c in seq(1, 21, 1)) {
  print(c)
  c.df <- C1.bin50.celltype.df[C1.bin50.celltype.df$cell_type==C1.cell_type.levels[c], ]
  plot <- ggplot() + 
    geom_point(C1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.1) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = C1.color[c], size = 0.3) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(C1.cell_type.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E85-C1-bin50-subtype.pdf", width = 12, height = 11)
wrap_plots(plots, ncol = 6)
dev.off()


#E95
load("20210818-spatial/E95-bin50-spatial-phago-map-remove.RData")
E95.bin50.celltype.df <- read.table("20210818-spatial/E95-bin50-spatial-remove-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E95.bin50.celltype.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type != "Erythroid", ]
E95.cell_type.levels <- c(decidual.levels[c(1:7)], immune.levels[c(1:7,9,10)], EC.levels[1:5], major.levels[c(4:5,11)], tropho.levels[c(1:5)], "embryo")
E95.color <- c(decidual.color[c(1:7)], immune.color[c(1:7,9,10)], EC.color[1:5], major.color[c(4:5,11)], tropho.color[c(1:5)], "lightgray")
E95.bin50.celltype.df$cell_type <- factor(E95.bin50.celltype.df$cell_type, levels = E95.cell_type.levels)
pdf("20210827-figures/figure2-spatial/Figure2-01-E95-bin50-all-types.pdf", width = 9, height = 6.5)
ggplot() + 
  geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = E95.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()
#subtype
plots <- list()
for (c in seq(1, 30, 1)) {
  print(E95.cell_type.levels[c])
  c.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type==E95.cell_type.levels[c], ]
  plot <- ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.1) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = E95.color[c], size = 0.3) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(E95.cell_type.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E95-bin50-subtype.pdf", width = 20, height = 20)
wrap_plots(plots, ncol = 6)
dev.off()
#E105
load("20210818-spatial/E105-bin50-spatial-phago-map.RData")
E105.bin50.celltype.df <- read.table("20210818-spatial/E105-bin50-spatial-phago-celltype-df.csv", sep = ",", check.names = F, row.names = 1, header = T)
E105.cell_type.levels <- c(decidual.levels[c(1:8)], immune.levels[c(1:8,10)], EC.levels[c(1:2, 4:5)], major.levels[4:5], major.levels[c(7:13)], tropho.levels[c(1:5)])
E105.color <- c(decidual.color[c(1:8)], immune.color[c(1:8,10)], EC.color[c(1:2, 4:5)], major.color[c(4:5,7:13)], tropho.color[c(1:5)])
E105.bin50.celltype.df$cell_type <- factor(E105.bin50.celltype.df$cell_type, levels = E105.cell_type.levels)
ggplot() + 
  geom_point(E105.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.01) +
  scale_color_manual(values = E105.color) + 
  theme_classic() + 
  theme(axis.text = element_blank()) + 
  facet_wrap(~cell_type, ncol = 7)
ggplot() + 
  geom_point(E105.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = E105.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())

#color spatial_clusters, UMAP and spatial slides
spatial.object <- list(E65.1.bin50, E65.2.bin50, E75.1.bin50, E75.2.bin50, B6.bin50, C1.bin50, E95.bin50)
names(spatial.object) <- c("E65.1.bin50", "E65.2.bin50", "E75.1.bin50", "E75.2.bin50", "E85.B6.bin50", "E85.C1.bin50", "E95.bin50")
plots <- list()
for (o in seq(1, 7, 1)) {
  print(o)
  plot <- DimPlot(spatial.object[[o]], group.by = "seurat_clusters", pt.size = 1) + 
    ggtitle(names(spatial.object)[[o]])
  plots[[o]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-02-UMAP-color-spatial-cluster.pdf", width = 21, height = 40)
wrap_plots(plots, ncol = 2)
dev.off()

pdf("20210827-figures/figure2-spatial/Figure2-02-E65-1-color-spatial-cluster-gradient.pdf", width = 3, height = 2)
SpatialDimPlot(E65.1.bin50, stroke = 0, pt.size.factor = 3) + NoLegend()
dev.off()
pdf("20210827-figures/figure2-spatial/Figure2-02-E65-2-color-spatial-cluster-gradient.pdf", width = 6, height = 5)
SpatialDimPlot(E65.2.bin50, stroke = 0, pt.size.factor = 2.5) + NoLegend()
dev.off()
pdf("20210827-figures/figure2-spatial/Figure2-02-E75-1-color-spatial-cluster-gradient.pdf", width = 6, height = 5)
SpatialDimPlot(E75.1.bin50, stroke = 0, pt.size.factor = 2) + NoLegend()
dev.off()
pdf("20210827-figures/figure2-spatial/Figure2-02-E75-2-color-spatial-cluster-gradient.pdf", width = 6, height = 5)
SpatialDimPlot(E75.2.bin50, stroke = 0, pt.size.factor = 1.8) + NoLegend()
dev.off()
load("20210818-spatial/B6-bin50-spatial.RData")
pdf("20210827-figures/figure2-spatial/Figure2-02-E85-B6-color-spatial-cluster-gradient.pdf", width = 6, height = 5)
SpatialDimPlot(B6.bin50, stroke = 0, pt.size.factor = 2.1) + NoLegend()
dev.off()
pdf("20210827-figures/figure2-spatial/Figure2-02-E85-c1-color-spatial-cluster-gradient.pdf", width = 6, height = 5)
SpatialDimPlot(C1.bin50, stroke = 0, pt.size.factor = 2.2) + NoLegend()
dev.off()
pdf("20210827-figures/figure2-spatial/Figure2-02-E95-color-spatial-cluster.pdf", width = 8, height = 11)
SpatialDimPlot(E95.bin50, stroke = 0, pt.size.factor = 1.4) + NoLegend()
dev.off()

#做bin20的图，各个小图
#E65.2
load("20210818-spatial/E65-2-bin20-spatial-phago-map.RData")
E65.2.bin20.celltype.df <- read.table("20210818-spatial/E65-2-bin20-spatial-phago-celltype-df.csv", sep = ",", header = T, check.names = F)
#E65.2.bin20.cell_type.levels <- c(decidual.levels[1:8], immune.levels[c(1:10)], EC.levels[1:5], major.levels[4:5], major.levels[c(7:13)], tropho.levels)
#E65.2.color <- c(decidual.color[1:8], immune.color[c(1:10)], EC.color[1:5], major.color[c(4:5,7:13)], tropho.color)
#E65.2.bin20.celltype.df$cell_type <- factor(E65.2.bin20.celltype.df$cell_type, levels = E65.2.bin20.cell_type.levels)
E65.2.bin20.celltype.df.immune <- E65.2.bin20.celltype.df[E65.2.bin20.celltype.df$cell_type %in% immune.levels[c(1:10)], ]
E65.2.bin20.celltype.df.immune$cell_type <- factor(E65.2.bin20.celltype.df.immune$cell_type, levels = immune.levels[c(1:10)])
#subtype
plots <- list()
for (c in seq(1, 10, 1)) {
  c.df <- E65.2.bin20.celltype.df.immune[E65.2.bin20.celltype.df.immune$cell_type==immune.levels[c], ]
  plot <- ggplot() + 
    geom_point(E65.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = immune.color[c], size = 0.4) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E65-2-bin20-subtype-immune-pt4.pdf", width = 15, height = 15)
wrap_plots(plots, ncol = 4)
dev.off()

#E65.1
load("20210818-spatial/E65-1-bin20-spatial-phago-map.RData")
E65.1.bin20.celltype.df <- read.table("20210818-spatial/E65-1-bin20-spatial-phago-celltype-df.csv", sep = ",", header = T, check.names = F)
E65.1.bin20.cell_type.levels <- c(decidual.levels[1:8], immune.levels[c(1:10)], EC.levels[1:5], major.levels[4:5], major.levels[c(7:13)], tropho.levels)
E65.1.color <- c(decidual.color[1:8], immune.color[c(1:10)], EC.color[1:5], major.color[c(4:5,7:13)], tropho.color)
E65.1.bin20.celltype.df$cell_type <- factor(E65.1.bin20.celltype.df$cell_type, levels = E65.1.bin20.cell_type.levels)
E65.1.bin20.celltype.df.immune <- E65.1.bin20.celltype.df[E65.1.bin20.celltype.df$cell_type %in% immune.levels[c(1:10)], ]
E65.1.bin20.celltype.df.immune$cell_type <- factor(E65.1.bin20.celltype.df.immune$cell_type, levels = immune.levels[c(1:10)])
#subtype
plots <- list()
for (c in seq(1, 10, 1)) {
  c.df <- E65.1.bin20.celltype.df.immune[E65.1.bin20.celltype.df.immune$cell_type==immune.levels[c], ]
  plot <- ggplot() + 
    geom_point(E65.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = immune.color[c], size = 0.4) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E65-1-bin20-subtype-immune-pt4.pdf", width = 15, height = 15)
wrap_plots(plots, ncol = 4)
dev.off()

#E75.2
load("20210818-spatial/E75-2-bin20-spatial-phago-map.RData")
E75.2.bin20.celltype.df <- read.table("20210818-spatial/E75-2-bin20-spatial-phago-celltype-df.csv", sep = ",", header = T, check.names = F)
E75.2.bin20.cell_type.levels <- c(decidual.levels[1:8], immune.levels[c(1:10)], EC.levels[1:5], major.levels[4:5], major.levels[c(7:13)], tropho.levels)
E75.2.color <- c(decidual.color[1:8], immune.color[c(1:10)], EC.color[1:5], major.color[c(4:5,7:13)], tropho.color)
E75.2.bin20.celltype.df$cell_type <- factor(E75.2.bin20.celltype.df$cell_type, levels = E75.2.bin20.cell_type.levels)
E75.2.bin20.celltype.df.immune <- E75.2.bin20.celltype.df[E75.2.bin20.celltype.df$cell_type %in% immune.levels[c(1:10)], ]
E75.2.bin20.celltype.df.immune$cell_type <- factor(E75.2.bin20.celltype.df.immune$cell_type, levels = immune.levels[c(1:10)])

#subtype
plots <- list()
for (c in seq(1, 10, 1)) {
  c.df <- E75.2.bin20.celltype.df.immune[E75.2.bin20.celltype.df.immune$cell_type==immune.levels[c], ]
  plot <- ggplot() + 
    geom_point(E75.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = immune.color[c], size = 0.4) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E75-2-bin20-subtype-immune-pt4.pdf", width = 15, height = 15)
wrap_plots(plots, ncol = 4)
dev.off()

#E75.1.bin20
load("20210818-spatial/E75-1-bin20-spatial-phago-map.RData")
E75.1.bin20.celltype.df <- read.table("20210818-spatial/E75-1-bin20-spatial-phago-celltype-df.csv", sep = ",", header = T, check.names = F)
E75.1.bin20.cell_type.levels <- c(decidual.levels[1:8], immune.levels[c(1:10)], EC.levels[1:5], major.levels[4:5], major.levels[c(7:13)], tropho.levels[1:5])
E75.1.color <- c(decidual.color[1:8], immune.color[c(1:10)], EC.color[1:5], major.color[c(4:5,7:13)], tropho.color[1:5])
E75.1.bin20.celltype.df$cell_type <- factor(E75.1.bin20.celltype.df$cell_type, levels = E75.1.bin20.cell_type.levels)
E75.1.bin20.celltype.df.immune <- E75.1.bin20.celltype.df[E75.1.bin20.celltype.df$cell_type %in% immune.levels[c(1:10)], ]
E75.1.bin20.celltype.df.immune$cell_type <- factor(E75.1.bin20.celltype.df.immune$cell_type, levels = immune.levels[c(1:10)])

#subtype
plots <- list()
for (c in seq(1, 10, 1)) {
  c.df <- E75.1.bin20.celltype.df.immune[E75.1.bin20.celltype.df.immune$cell_type==immune.levels[c], ]
  plot <- ggplot() + 
    geom_point(E75.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = immune.color[c], size = 0.4) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E75-1-bin20-subtype-immune-pt4.pdf", width = 15, height = 15)
wrap_plots(plots, ncol = 4)
dev.off()

#E85.B6
load("20210818-spatial/E85-B6-bin20-spatial-phago-map.RData")
E85.B6.bin20.celltype.df <- read.table("20210818-spatial/E85-B6-bin20-spatial-phago-celltype-df.csv", sep = ",", header = T, check.names = F)
E85.B6.bin20.cell_type.levels <- c(decidual.levels[1:8], immune.levels[c(1:10)], EC.levels[1:5], major.levels[4:5], major.levels[c(7:11, 13)], tropho.levels[1:5])
E85.B6.color <- c(decidual.color[1:8], immune.color[c(1:10)], EC.color[1:5], major.color[c(4:5,7:11, 13)], tropho.color[1:5])
E85.B6.bin20.celltype.df$cell_type <- factor(E85.B6.bin20.celltype.df$cell_type, levels = E85.B6.bin20.cell_type.levels)
E85.B6.bin20.celltype.df.immune <- E85.B6.bin20.celltype.df[E85.B6.bin20.celltype.df$cell_type %in% immune.levels[c(1:10)], ]
E85.B6.bin20.celltype.df.immune$cell_type <- factor(E85.B6.bin20.celltype.df.immune$cell_type, levels = immune.levels[c(1:10)])

#subtype
plots <- list()
for (c in seq(1, 10, 1)) {
  c.df <- E85.B6.bin20.celltype.df.immune[E85.B6.bin20.celltype.df.immune$cell_type==immune.levels[c], ]
  plot <- ggplot() + 
    geom_point(E85.B6.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = immune.color[c], size = 0.4) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E85-B6-bin20-subtype-immune-pt4.pdf", width = 15, height = 15)
wrap_plots(plots, ncol = 4)
dev.off()


#E85.C1
load("20210818-spatial/E85-C1-bin20-spatial-phago-map.RData")
E85.C1.bin20.celltype.df <- read.table("20210818-spatial/E85-C1-bin20-spatial-phago-celltype-df.csv", sep = ",", header = T, check.names = F)
E85.C1.bin20.cell_type.levels <- c(decidual.levels[1:8], immune.levels[c(1:4, 6:10)], EC.levels[1:5], major.levels[4:5], major.levels[c(7:11, 13)], tropho.levels[1:5])
E85.C1.color <- c(decidual.color[1:8], immune.color[c(1:4, 6:10)], EC.color[1:5], major.color[c(4:5,7:11, 13)], tropho.color[1:5])
E85.C1.bin20.celltype.df$cell_type <- factor(E85.C1.bin20.celltype.df$cell_type, levels = E85.C1.bin20.cell_type.levels)
E85.C1.bin20.celltype.df.immune <- E85.C1.bin20.celltype.df[E85.C1.bin20.celltype.df$cell_type %in% immune.levels[c(1:10)], ]
E85.C1.bin20.celltype.df.immune$cell_type <- factor(E85.C1.bin20.celltype.df.immune$cell_type, levels = immune.levels[c(1:10)])

#subtype
plots <- list()
for (c in seq(1, 10, 1)) {
  c.df <- E85.C1.bin20.celltype.df.immune[E85.C1.bin20.celltype.df.immune$cell_type==immune.levels[c], ]
  plot <- ggplot() + 
    geom_point(E85.C1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = immune.color[c], size = 0.4) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E85-C1-bin20-subtype-immune-pt4.pdf", width = 15, height = 15)
wrap_plots(plots, ncol = 4)
dev.off()


#E95
load("20210818-spatial/E85-B6-bin20-spatial-phago-map.RData")
E95.bin20.celltype.df <- read.table("20210818-spatial/E95-bin20-spatial-phago-celltype-df.csv", sep = ",", header = T, check.names = F)
E95.bin20.cell_type.levels <- c(decidual.levels[1:8], immune.levels[c(1:10)], EC.levels[1:5], major.levels[4:5], major.levels[c(7:13)], tropho.levels)
E95.color <- c(decidual.color[1:8], immune.color[c(1:10)], EC.color[1:5], major.color[c(4:5,7:13)], tropho.color)
E95.bin20.celltype.df$cell_type <- factor(E95.bin20.celltype.df$cell_type, levels = E95.bin20.cell_type.levels)

E95.bin20.celltype.df.immune <- E95.bin20.celltype.df[E95.bin20.celltype.df$cell_type %in% immune.levels[c(1:10)], ]
E95.bin20.celltype.df.immune$cell_type <- factor(E95.bin20.celltype.df.immune$cell_type, levels = immune.levels[c(1:10)])

#subtype
plots <- list()
for (c in seq(1, 10, 1)) {
  c.df <- E95.bin20.celltype.df.immune[E95.bin20.celltype.df.immune$cell_type==immune.levels[c], ]
  plot <- ggplot() + 
    geom_point(E95.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(c.df, mapping = aes(x = x, y = y), color = immune.color[c], size = 0.4) + 
    theme_classic() + 
    theme(axis.text = element_blank()) +
    ggtitle(immune.levels[c])
  plots[[c]] <- plot
}
pdf("20210827-figures/figure2-spatial/Figure2-01-E95-bin20-subtype-immune-pt4.pdf", width = 20, height = 15)
wrap_plots(plots, ncol = 4)
dev.off()

#bin50各个片子的基因的表达
decidual.genes <- c("Lum", "Dcn", "Sfrp4", "Igfbp3", "Top2a", "Ube2c", "Ptn", "Bmp2", 
                    "Gatm", "Anxa8", "S100a8", "Angpt2", "Prl8a2", "Erv3", "Rxfp1", "Wnt4", "Hand2")
SpatialFeaturePlot(E65.2.bin50, features = decidual.genes)
SpatialFeaturePlot(E75.2.bin50, features = decidual.genes)
SpatialFeaturePlot(B6.bin50, features = decidual.genes)
SpatialFeaturePlot(E95.bin50, features = decidual.genes)

SpatialFeaturePlot(E65.2.bin20, features = decidual.genes)
SpatialFeaturePlot(E75.2.bin20, features = decidual.genes)
SpatialFeaturePlot(B6.bin20, features = decidual.genes)
SpatialFeaturePlot(E95.bin20, features = decidual.genes)

EC.genes <- c("Plvap", "Igfbp7", "Akr1c18", "Cldn5", "Ube2c", "Hmgb2", "Gpx3", "Tdo2", "Slc6a2", "Stmn2", "Lyve1", "Fgl2")
SpatialFeaturePlot(E65.2.bin50, features = EC.genes)
SpatialFeaturePlot(E75.2.bin50, features = EC.genes)
SpatialFeaturePlot(B6.bin50, features = EC.genes)
SpatialFeaturePlot(E95.bin50, features = EC.genes)

#QC 所有片子
pdf("20210827-figures/figure2-spatial/Figure2-03-QC-E65-1.pdf", width = 20, height = 8)
SpatialFeaturePlot(E65.1.bin50, features = c("nFeature_Spatial", "nCount_Spatial"), pt.size.factor = 2.5)
dev.off()
pdf("20210827-figures/figure2-spatial/Figure2-03-QC-E65-2.pdf", width = 20, height = 8)
SpatialFeaturePlot(E65.2.bin50, features = c("nFeature_Spatial", "nCount_Spatial"), pt.size.factor = 2)
dev.off()
pdf("20210827-figures/figure2-spatial/Figure2-03-QC-E75-1.pdf", width = 20, height = 8)
SpatialFeaturePlot(E75.1.bin50, features = c("nFeature_Spatial", "nCount_Spatial"))
dev.off()
pdf("20210827-figures/figure2-spatial/Figure2-03-QC-E75-2.pdf", width = 20, height = 8)
SpatialFeaturePlot(E75.2.bin50, features = c("nFeature_Spatial", "nCount_Spatial"))
dev.off()
pdf("20210827-figures/figure2-spatial/Figure2-03-QC-E95.pdf", width = 20, height = 8)
SpatialFeaturePlot(E95.bin50, features = c("nFeature_Spatial", "nCount_Spatial"))
dev.off()


#所有片子的灰色值
E65.1.bin50.celltype.df$sample <- "E65.1"
E65.2.bin50.celltype.df$sample <- "E65.2"
E75.1.bin50.celltype.df$sample <- "E75.1"
E75.2.bin50.celltype.df$sample <- "E75.2"
B6.bin50.celltype.df$sample <- "E85.B6"
C1.bin50.celltype.df$sample <- "E85.C1"
E95.bin50.celltype.df$sample <- "E95"
all.slides.def <- Reduce(rbind, list(E65.1.bin50.celltype.df, E65.2.bin50.celltype.df, E75.1.bin50.celltype.df, E75.2.bin50.celltype.df,
                                     B6.bin50.celltype.df, C1.bin50.celltype.df, E95.bin50.celltype.df))
write.table(all.slides.def, file = "20210818-spatial/All-slides-bin50-spatial-remove-celltype-df.csv", sep = ",", col.names = T, row.names = F, quote = F)
pdf("20210827-figures/figure2-spatial/Figure2-04-all-slides-gray-background.pdf", width = 8, height = 10)
ggplot() + 
  geom_point(all.slides.def, mapping = aes(x = x, y = y), size = 0.3, color = "#E7E7E7") + 
  facet_wrap(~ sample, ncol = 3, scales = "free") +
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#所有片子的enrichment score bin50的
E65.1.enrichment.score <- function(){
  #E65.1
  #enrichment score
  DefaultAssay(E65.1.bin50) <- "predictions"
  E65.1.bin50.celltype.score <- GetAssayData(E65.1.bin50, slot = "data")
  E65.1.bin50.celltype.max <- apply(E65.1.bin50.celltype.score, 2, function(x) rownames(E65.1.bin50.celltype.score)[which.max(x)])
  E65.1.bin50.celltype.sort <- E65.1.bin50.celltype.score[, names(sort(E65.1.bin50.celltype.max))]
  col_anno <- data.frame(cell_type = E65.1.bin50.celltype.max)
  rownames(col_anno) <- names(E65.1.bin50.celltype.max)
  col_anno$cell_type <- factor(col_anno$cell_type, levels = E65.1.cell_type.levels)
  col_anno$cell_name <- rownames(col_anno)
  col_anno.sort <- data.frame()
  for (n in E65.1.cell_type.levels) {
    col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
  }
  col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type)
  rownames(col_anno.df) <- rownames(col_anno.sort)
  color <- brewer.pal(n =11, name = "RdBu")
  color <- color[c(1:4, 7:11)]
  color = colorRampPalette(rev(color))(100)
  E65.1.bin50.celltype.sort.s <- E65.1.bin50.celltype.sort[E65.1.cell_type.levels, rownames(col_anno.df)]
  pdf("20210827-figures/figure2-spatial/Figure2-E65.1-seurat-score-heatmap.pdf", width = 20, height = 10)
  pheatmap(E65.1.bin50.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
           scale = "column", show_colnames = F, annotation_col = col_anno.df, 
           color = color)
  dev.off()
}
E65.2.enrichment.score <- function(){
  DefaultAssay(E65.2.bin50) <- "predictions"
  E65.2.bin50.celltype.score <- GetAssayData(E65.2.bin50, slot = "data")
  E65.2.bin50.celltype.max <- apply(E65.2.bin50.celltype.score, 2, function(x) rownames(E65.2.bin50.celltype.score)[which.max(x)])
  E65.2.bin50.celltype.sort <- E65.2.bin50.celltype.score[, names(sort(E65.2.bin50.celltype.max))]
  col_anno <- data.frame(cell_type = E65.2.bin50.celltype.max)
  rownames(col_anno) <- names(E65.2.bin50.celltype.max)
  col_anno$cell_type <- factor(col_anno$cell_type, levels = E65.2.cell_type.levels)
  col_anno$cell_name <- rownames(col_anno)
  col_anno.sort <- data.frame()
  for (n in E65.2.cell_type.levels) {
    col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
  }
  col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type)
  rownames(col_anno.df) <- rownames(col_anno.sort)
  color <- brewer.pal(n =11, name = "RdBu")
  color <- color[c(1:4, 7:11)]
  color = colorRampPalette(rev(color))(100)
  E65.2.bin50.celltype.sort.s <- E65.2.bin50.celltype.sort[E65.2.cell_type.levels, rownames(col_anno.df)]
  pdf("20210827-figures/figure2-spatial/Figure2-E65.2-seurat-score-heatmap.pdf", width = 20, height = 10)
  pheatmap(E65.2.bin50.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
           scale = "column", show_colnames = F, annotation_col = col_anno.df, 
           color = color)
  dev.off()
}
E75.1.enrichment.score <- function(){
  #E75.1
  #enrichment score
  DefaultAssay(E75.1.bin50) <- "predictions"
  E75.1.bin50.celltype.score <- GetAssayData(E75.1.bin50, slot = "data")
  E75.1.bin50.celltype.max <- apply(E75.1.bin50.celltype.score, 2, function(x) rownames(E75.1.bin50.celltype.score)[which.max(x)])
  E75.1.bin50.celltype.sort <- E75.1.bin50.celltype.score[, names(sort(E75.1.bin50.celltype.max))]
  col_anno <- data.frame(cell_type = E75.1.bin50.celltype.max)
  rownames(col_anno) <- names(E75.1.bin50.celltype.max)
  col_anno$cell_type <- factor(col_anno$cell_type, levels = E75.1.cell_type.levels)
  col_anno$cell_name <- rownames(col_anno)
  col_anno.sort <- data.frame()
  for (n in E75.1.cell_type.levels) {
    col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
  }
  col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type)
  rownames(col_anno.df) <- rownames(col_anno.sort)
  color <- brewer.pal(n =11, name = "RdBu")
  color <- color[c(1:4, 7:11)]
  color = colorRampPalette(rev(color))(100)
  E75.1.bin50.celltype.sort.s <- E75.1.bin50.celltype.sort[E75.1.cell_type.levels, rownames(col_anno.df)]
  pdf("20210827-figures/figure2-spatial/Figure2-E75.1-seurat-score-heatmap.pdf", width = 20, height = 10)
  pheatmap(E75.1.bin50.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
           scale = "column", show_colnames = F, annotation_col = col_anno.df, 
           color = color)
  dev.off()
}
E75.2.enrichment.score <- function(){
  DefaultAssay(E75.2.bin50) <- "predictions"
  E75.2.bin50.celltype.score <- GetAssayData(E75.2.bin50, slot = "data")
  E75.2.bin50.celltype.max <- apply(E75.2.bin50.celltype.score, 2, function(x) rownames(E75.2.bin50.celltype.score)[which.max(x)])
  E75.2.bin50.celltype.sort <- E75.2.bin50.celltype.score[, names(sort(E75.2.bin50.celltype.max))]
  col_anno <- data.frame(cell_type = E75.2.bin50.celltype.max)
  rownames(col_anno) <- names(E75.2.bin50.celltype.max)
  col_anno$cell_type <- factor(col_anno$cell_type, levels = E75.2.cell_type.levels)
  col_anno$cell_name <- rownames(col_anno)
  col_anno.sort <- data.frame()
  for (n in E75.2.cell_type.levels) {
    col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
  }
  col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type)
  rownames(col_anno.df) <- rownames(col_anno.sort)
  color <- brewer.pal(n =11, name = "RdBu")
  color <- color[c(1:4, 7:11)]
  color = colorRampPalette(rev(color))(100)
  E75.2.bin50.celltype.sort.s <- E75.2.bin50.celltype.sort[E75.2.cell_type.levels, rownames(col_anno.df)]
  pdf("20210827-figures/figure2-spatial/Figure2-E75.2-seurat-score-heatmap.pdf", width = 20, height = 10)
  pheatmap(E75.2.bin50.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
           scale = "column", show_colnames = F, annotation_col = col_anno.df, 
           color = color)
  dev.off()
}
E95.enrichment.score <- function(){
  #E95
  #enrichment score
  DefaultAssay(E95.bin50) <- "predictions"
  E95.bin50.celltype.score <- GetAssayData(E95.bin50, slot = "data")
  E95.bin50.celltype.max <- apply(E95.bin50.celltype.score, 2, function(x) rownames(E95.bin50.celltype.score)[which.max(x)])
  E95.bin50.celltype.sort <- E95.bin50.celltype.score[, names(sort(E95.bin50.celltype.max))]
  col_anno <- data.frame(cell_type = E95.bin50.celltype.max)
  rownames(col_anno) <- names(E95.bin50.celltype.max)
  col_anno$cell_type <- factor(col_anno$cell_type, levels = E95.cell_type.levels)
  col_anno$cell_name <- rownames(col_anno)
  col_anno.sort <- data.frame()
  for (n in E95.cell_type.levels) {
    col_anno.sort <- rbind(col_anno.sort, col_anno[col_anno$cell_type==n, ])
  }
  col_anno.df <- data.frame(cell_type = col_anno.sort$cell_type)
  rownames(col_anno.df) <- rownames(col_anno.sort)
  color <- brewer.pal(n =11, name = "RdBu")
  color <- color[c(1:4, 7:11)]
  color = colorRampPalette(rev(color))(100)
  E95.bin50.celltype.sort.s <- E95.bin50.celltype.sort[E95.cell_type.levels, rownames(col_anno.df)]
  pdf("20210827-figures/figure2-spatial/Figure2-E95-seurat-score-heatmap.pdf", width = 20, height = 10)
  pheatmap(E95.bin50.celltype.sort.s, cluster_rows = F, cluster_cols = F, annotation_colors = annotation.color,
           scale = "column", show_colnames = F, annotation_col = col_anno.df, 
           color = color)
  dev.off()
}
all.slide.score.2 <- merge(E65.2.bin50, list(E65.1.bin50, E75.2.bin50, E75.1.bin50, B6.bin50, C1.bin50, E95.bin50), 
                           add.cell.ids = c("E65.2.bin50", "E65.1.bin50", "E75.2.bin50", "E75.1.bin50", "E85.B6.bin50", "E85.C1.bin50", "E95.bin50"))
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


#每个时间点的空间上每个zone的比例统计
#decidual, EC, immune, FB-immune, trophoblast
#统计E65.2, E75.1, E85.B6, E95
all.slides.def.subset <- all.slides.def[all.slides.def$sample %in% c("E65.2",  "E75.1",  "E85.B6", "E95"), ]
#MET
MET.cell <- c("Lum-D", "Sfrp4-D")
all.slides.def.subset.MET <- all.slides.def.subset[all.slides.def.subset$cell_type %in% MET.cell, ]
all.slides.def.subset.MET.df <- table(all.slides.def.subset.MET$cell_type, all.slides.def.subset.MET$sample)
all.slides.def.subset.MET.df <- all.slides.def.subset.MET.df[MET.cell, ]
all.slides.def.subset.MET.df <- apply(all.slides.def.subset.MET.df, 1, function(x){x/table(all.slides.def.subset$sample) * 100})
all.slides.def.subset.MET.df <- melt(all.slides.def.subset.MET.df)
all.slides.def.subset.MET.df$zone <- "MET"

#cycling
cycling.cell <- c("Top2a-D")
all.slides.def.subset.cycling <- all.slides.def.subset[all.slides.def.subset$cell_type %in% cycling.cell, ]
all.slides.def.subset.cycling.df <- table(all.slides.def.subset.cycling$cell_type, all.slides.def.subset.cycling$sample)
all.slides.def.subset.cycling.df <- all.slides.def.subset.cycling.df[cycling.cell, ]
all.slides.def.subset.cycling.df <- all.slides.def.subset.cycling.df/table(all.slides.def.subset$sample) * 100
all.slides.def.subset.cycling.df <- data.frame(Var1 = names(all.slides.def.subset.cycling.df), Var2 = "Top2a-D", 
                                          value = as.numeric(all.slides.def.subset.cycling.df), zone = "cycling DSC")
#nourish
nourish.cell <- c("Gatm-D", "Prl8a2-D")
all.slides.def.subset.nourish <- all.slides.def.subset[all.slides.def.subset$cell_type %in% nourish.cell, ]
all.slides.def.subset.nourish.df <- table(all.slides.def.subset.nourish$cell_type, all.slides.def.subset.nourish$sample)
all.slides.def.subset.nourish.df <- all.slides.def.subset.nourish.df[nourish.cell, ]
all.slides.def.subset.nourish.df <- apply(all.slides.def.subset.nourish.df, 1, function(x){x/table(all.slides.def.subset$sample) * 100})
all.slides.def.subset.nourish.df <- melt(all.slides.def.subset.nourish.df)
all.slides.def.subset.nourish.df$zone <- "nourish"


immune.cell <- c("Mac", "Pro-Mac", "Neutrophil", "DC-1", "DC-2", "Mono", "NK")
all.slides.def.subset.immune <- all.slides.def.subset[all.slides.def.subset$cell_type %in% immune.cell, ]
all.slides.def.subset.immune.df <- table(all.slides.def.subset.immune$cell_type, all.slides.def.subset.immune$sample)
all.slides.def.subset.immune.df <- all.slides.def.subset.immune.df[immune.cell, ]
all.slides.def.subset.immune.df <- apply(all.slides.def.subset.immune.df, 1, function(x){x/table(all.slides.def.subset$sample) * 100})
all.slides.def.subset.immune.df <- melt(all.slides.def.subset.immune.df)
all.slides.def.subset.immune.df$zone <- "immune"

EC.cell <- c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D")
all.slides.def.subset.EC <- all.slides.def.subset[all.slides.def.subset$cell_type %in% EC.cell, ]
all.slides.def.subset.EC.df <- table(all.slides.def.subset.EC$cell_type, all.slides.def.subset.EC$sample)
all.slides.def.subset.EC.df <- all.slides.def.subset.EC.df[EC.cell, ]
all.slides.def.subset.EC.df <- apply(all.slides.def.subset.EC.df, 1, function(x){x/table(all.slides.def.subset$sample) * 100})
all.slides.def.subset.EC.df <- melt(all.slides.def.subset.EC.df)
all.slides.def.subset.EC.df$zone <- "EC"

all.slides.def.subset.FB.immune <- all.slides.def.subset[all.slides.def.subset$cell_type %in% c("FB-immune"), ]
all.slides.def.subset.FB.immune.df <- table(all.slides.def.subset.FB.immune$cell_type, all.slides.def.subset.FB.immune$sample)
all.slides.def.subset.FB.immune.df <- all.slides.def.subset.FB.immune.df["FB-immune", ]/table(all.slides.def.subset$sample) * 100
all.slides.def.subset.FB.immune.df <- data.frame(Var1 = names(all.slides.def.subset.FB.immune.df), Var2 = "FB-immune", 
                                                 value = as.numeric(all.slides.def.subset.FB.immune.df), zone = "FB-immune")

placenta.cell <- c("tropho-1", "tropho-2", "tropho-3", "tropho-0")
all.slides.def.subset.placenta <- all.slides.def.subset[all.slides.def.subset$cell_type %in% placenta.cell, ]
all.slides.def.subset.placenta.df <- table(all.slides.def.subset.placenta$cell_type, all.slides.def.subset.placenta$sample)
all.slides.def.subset.placenta.df <- all.slides.def.subset.placenta.df[placenta.cell, ]
all.slides.def.subset.placenta.df<- cbind(E65.2 = rep(0, 4), E75.1 = rep(0, 4), all.slides.def.subset.placenta.df)
all.slides.def.subset.placenta.df <- apply(all.slides.def.subset.placenta.df, 1, function(x){x/table(all.slides.def.subset$sample) * 100})
all.slides.def.subset.placenta.df <- melt(all.slides.def.subset.placenta.df)
all.slides.def.subset.placenta.df$zone <- "placenta"

all.slides.def.subset.giant <- all.slides.def.subset[all.slides.def.subset$cell_type %in% c("tropho-4"), ]
all.slides.def.subset.giant.df <- table(all.slides.def.subset.giant$cell_type, all.slides.def.subset.giant$sample)
all.slides.def.subset.giant.df <- c(0, 0, 243, 1330)
names(all.slides.def.subset.giant.df) <- c("E65.2", "E75.1", "E85.B6", "E95")
all.slides.def.subset.giant.df <- all.slides.def.subset.giant.df/table(all.slides.def.subset$sample) * 100
all.slides.def.subset.giant.df <- data.frame(Var1 = names(all.slides.def.subset.giant.df), Var2 = "tropho-4", 
                                          value = as.numeric(all.slides.def.subset.giant.df), zone = "Trophobalst Giant cell")


all.slides.def.subset.VE <- all.slides.def.subset[all.slides.def.subset$cell_type %in% c("Viseral endoderm"), ]
all.slides.def.subset.VE.df <- table(all.slides.def.subset.VE$cell_type, all.slides.def.subset.VE$sample)
all.slides.def.subset.VE.df <- c(0, 0, 153, 1857)
names(all.slides.def.subset.VE.df) <- c("E65.2", "E75.1", "E85.B6", "E95")
all.slides.def.subset.VE.df <- all.slides.def.subset.VE.df/table(all.slides.def.subset$sample) * 100
all.slides.def.subset.VE.df <- data.frame(Var1 = names(all.slides.def.subset.VE.df), Var2 = "Viseral endoderm", 
                                                 value = as.numeric(all.slides.def.subset.VE.df), zone = "Viseral endoderm")

all.slides.def.subset.zone.df <- Reduce(rbind, list(all.slides.def.subset.MET.df, all.slides.def.subset.cycling.df,
                                                    all.slides.def.subset.nourish.df, all.slides.def.subset.EC.df, 
                                                    all.slides.def.subset.FB.immune.df, all.slides.def.subset.immune.df, 
                                                    all.slides.def.subset.placenta.df, 
                                                    all.slides.def.subset.giant.df, all.slides.def.subset.VE.df))
all.slides.def.subset.zone.df$Var2 <- factor(all.slides.def.subset.zone.df$Var2, levels = c(decidual.levels, immune.levels, EC.levels, tropho.levels, "Viseral endoderm"))

pdf("20210827-figures/figure2-spatial/Figure2-05-spatial-zone-ratio.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(all.slides.def.subset.zone.df, mapping = aes(x = Var1, y = value, fill = Var2), stat = "identity", position = "stack") + 
  facet_wrap(~ zone, scales = "free") + 
  scale_fill_manual(values = c(decidual.color[c(1:3,5:7)], immune.color[c(1:7,10)], EC.color[1:5], tropho.color[1:5], "#9370DB")) + 
  theme_classic()
dev.off()

#统计E65.1, E75.1, E85.C1
all.slides.def.subset <- all.slides.def[all.slides.def$sample %in% c("E65.1",  "E75.2",  "E85.C1"), ]
#MET
MET.cell <- c("Lum-D", "Sfrp4-D")
all.slides.def.subset.MET <- all.slides.def.subset[all.slides.def.subset$cell_type %in% MET.cell, ]
all.slides.def.subset.MET.df <- table(all.slides.def.subset.MET$cell_type, all.slides.def.subset.MET$sample)
all.slides.def.subset.MET.df <- all.slides.def.subset.MET.df[MET.cell, ]
all.slides.def.subset.MET.df <- apply(all.slides.def.subset.MET.df, 1, function(x){x/table(all.slides.def.subset$sample) * 100})
all.slides.def.subset.MET.df <- melt(all.slides.def.subset.MET.df)
all.slides.def.subset.MET.df$zone <- "MET"

#cycling
cycling.cell <- c("Top2a-D")
all.slides.def.subset.cycling <- all.slides.def.subset[all.slides.def.subset$cell_type %in% cycling.cell, ]
all.slides.def.subset.cycling.df <- table(all.slides.def.subset.cycling$cell_type, all.slides.def.subset.cycling$sample)
all.slides.def.subset.cycling.df <- all.slides.def.subset.cycling.df[cycling.cell, ]
all.slides.def.subset.cycling.df <- all.slides.def.subset.cycling.df/table(all.slides.def.subset$sample) * 100
all.slides.def.subset.cycling.df <- data.frame(Var1 = names(all.slides.def.subset.cycling.df), Var2 = "Top2a-D", 
                                               value = as.numeric(all.slides.def.subset.cycling.df), zone = "cycling DSC")
#nourish
nourish.cell <- c("Gatm-D", "Prl8a2-D")
all.slides.def.subset.nourish <- all.slides.def.subset[all.slides.def.subset$cell_type %in% nourish.cell, ]
all.slides.def.subset.nourish.df <- table(all.slides.def.subset.nourish$cell_type, all.slides.def.subset.nourish$sample)
all.slides.def.subset.nourish.df <- all.slides.def.subset.nourish.df[nourish.cell, ]
all.slides.def.subset.nourish.df <- apply(all.slides.def.subset.nourish.df, 1, function(x){x/table(all.slides.def.subset$sample) * 100})
all.slides.def.subset.nourish.df <- melt(all.slides.def.subset.nourish.df)
all.slides.def.subset.nourish.df$zone <- "nourish"

all.slides.def.subset.immune <- all.slides.def.subset[all.slides.def.subset$cell_type %in% immune.cell, ]
all.slides.def.subset.immune.df <- table(all.slides.def.subset.immune$cell_type, all.slides.def.subset.immune$sample)
all.slides.def.subset.immune.df <- all.slides.def.subset.immune.df[immune.cell, ]
all.slides.def.subset.immune.df <- apply(all.slides.def.subset.immune.df, 1, function(x){x/table(all.slides.def.subset$sample) * 100})
all.slides.def.subset.immune.df <- melt(all.slides.def.subset.immune.df)
all.slides.def.subset.immune.df$zone <- "immune"

all.slides.def.subset.EC <- all.slides.def.subset[all.slides.def.subset$cell_type %in% EC.cell, ]
all.slides.def.subset.EC.df <- table(all.slides.def.subset.EC$cell_type, all.slides.def.subset.EC$sample)
all.slides.def.subset.EC.df <- all.slides.def.subset.EC.df[EC.cell, ]
all.slides.def.subset.EC.df <- apply(all.slides.def.subset.EC.df, 1, function(x){x/table(all.slides.def.subset$sample) * 100})
all.slides.def.subset.EC.df <- melt(all.slides.def.subset.EC.df)
all.slides.def.subset.EC.df$zone <- "EC"

all.slides.def.subset.FB.immune <- all.slides.def.subset[all.slides.def.subset$cell_type %in% c("FB-immune"), ]
all.slides.def.subset.FB.immune.df <- table(all.slides.def.subset.FB.immune$cell_type, all.slides.def.subset.FB.immune$sample)
all.slides.def.subset.FB.immune.df <- all.slides.def.subset.FB.immune.df["FB-immune", ]/table(all.slides.def.subset$sample) * 100
all.slides.def.subset.FB.immune.df <- data.frame(Var1 = names(all.slides.def.subset.FB.immune.df), Var2 = "FB-immune", 
                                                 value = as.numeric(all.slides.def.subset.FB.immune.df), zone = "FB-immune")

all.slides.def.subset.tropho <- all.slides.def.subset[all.slides.def.subset$cell_type %in% tropho.cell, ]
all.slides.def.subset.tropho.df <- table(all.slides.def.subset.tropho$cell_type, all.slides.def.subset.tropho$sample)
all.slides.def.subset.tropho.df <- all.slides.def.subset.tropho.df[tropho.cell, ]
all.slides.def.subset.tropho.df <- c(0, 320, 103)
names(all.slides.def.subset.tropho.df) <- c("E65.1", "E75.2", "E85.C1")
all.slides.def.subset.tropho.df <- all.slides.def.subset.tropho.df/table(all.slides.def.subset$sample) * 100
all.slides.def.subset.tropho.df <- data.frame(Var1 = names(all.slides.def.subset.tropho.df), Var2 = "tropho-4", 
                                              value = as.numeric(all.slides.def.subset.tropho.df), zone = "Trophoblast giant cell")

all.slides.def.subset.zone.df <- Reduce(rbind, list(all.slides.def.subset.MET.df, all.slides.def.subset.cycling.df, 
                                                    all.slides.def.subset.nourish.df, all.slides.def.subset.EC.df, 
                                                    all.slides.def.subset.FB.immune.df, all.slides.def.subset.immune.df, 
                                                    all.slides.def.subset.EC.df, all.slides.def.subset.tropho.df))
all.slides.def.subset.zone.df$Var2 <- factor(all.slides.def.subset.zone.df$Var2, levels = c(decidual.levels, immune.levels, EC.levels, tropho.levels))

pdf("20210827-figures/figure2-spatial/Figure2-05-spatial-zone-ratio-another-slide.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(all.slides.def.subset.zone.df, mapping = aes(x = Var1, y = value, fill = Var2), stat = "identity", position = "stack") + 
  facet_wrap(~ zone, scales = "free") + 
  scale_fill_manual(values = c(decidual.color[c(1:3,5:7)], immune.color[c(1:7,10)], EC.color[1:5], tropho.color[5])) + 
  theme_classic()
dev.off()

#bin20 immune cell
E65.1.bin20.celltype.df$sample <- "E65.1"
E65.2.bin20.celltype.df$sample <- "E65.2"
E75.1.bin20.celltype.df$sample <- "E75.1"
E75.2.bin20.celltype.df$sample <- "E75.2"
E85.B6.bin20.celltype.df$sample <- "E85.B6"
E85.C1.bin20.celltype.df$sample <- "E85.C1"
E85.C1.bin20.celltype.df <- E85.C1.bin20.celltype.df[, -c(7,8)]
E95.bin20.celltype.df$sample <- "E95"
all.slides.def.bin20 <- Reduce(rbind, list(E65.1.bin20.celltype.df, E65.2.bin20.celltype.df, E75.1.bin20.celltype.df, E75.2.bin20.celltype.df,
                                     E85.B6.bin20.celltype.df, E85.C1.bin20.celltype.df, E95.bin20.celltype.df))
#统计"E65.2",  "E75.1",  "E85.B6", "E95"
all.slides.def.bin20.subset <- all.slides.def.bin20[all.slides.def.bin20$sample %in% c("E65.2",  "E75.1",  "E85.B6", "E95"), ]
all.slides.def.bin20.subset.immune <- all.slides.def.bin20.subset[all.slides.def.bin20.subset$cell_type %in% c("Mac", "Pro-Mac", "Neutrophil", "DC-1", "DC-2", "Mono", "NK"), ]
all.slides.def.bin20.subset.immune.df <- table(all.slides.def.bin20.subset.immune$cell_type, all.slides.def.bin20.subset.immune$sample)
all.slides.def.bin20.subset.immune.df <- all.slides.def.bin20.subset.immune.df[c(9:15), ]
all.slides.def.bin20.subset.immune.df <- apply(all.slides.def.bin20.subset.immune.df, 1, function(x){x/table(all.slides.def.bin20.subset$sample) * 100})
all.slides.def.bin20.subset.immune.df <- melt(all.slides.def.bin20.subset.immune.df)
all.slides.def.bin20.subset.immune.df$zone <- "immune"
pdf("20210827-figures/figure2-spatial/Figure2-05-spatial-zone-ratio-immune-bin20.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(all.slides.def.bin20.subset.immune.df, mapping = aes(x = Var1, y = value, fill = Var2), stat = "identity", position = "stack") + 
  scale_fill_manual(values = c(immune.color[c(1:7,10)])) + 
  theme_classic()
dev.off()

#统计"E65.1",  "E75.2",  "E85.C1"
all.slides.def.bin20.subset <- all.slides.def.bin20[all.slides.def.bin20$sample %in% c("E65.1",  "E75.2",  "E85.C1"), ]
all.slides.def.bin20.subset.immune <- all.slides.def.bin20.subset[all.slides.def.bin20.subset$cell_type %in% c("Mac", "Pro-Mac", "Neutrophil", "DC-1", "DC-2", "Mono", "NK"), ]
all.slides.def.bin20.subset.immune.df <- table(all.slides.def.bin20.subset.immune$cell_type, all.slides.def.bin20.subset.immune$sample)
all.slides.def.bin20.subset.immune.df <- all.slides.def.bin20.subset.immune.df[c(9:15), ]
all.slides.def.bin20.subset.immune.df <- apply(all.slides.def.bin20.subset.immune.df, 1, function(x){x/table(all.slides.def.bin20.subset$sample) * 100})
all.slides.def.bin20.subset.immune.df <- melt(all.slides.def.bin20.subset.immune.df)
all.slides.def.bin20.subset.immune.df$zone <- "immune"
pdf("20210827-figures/figure2-spatial/Figure2-05-spatial-zone-ratio-another-slide-immung-bin20.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(all.slides.def.bin20.subset.immune.df, mapping = aes(x = Var1, y = value, fill = Var2), stat = "identity", position = "stack") + 
  scale_fill_manual(values = c(immune.color[c(1:7,10)])) + 
  theme_classic()
dev.off()

MET.gene <- c("Lum", "Dcn", "Col3a1", "Col1a1", "Sfrp4", "Igfbp3", "Igfbp5", "Bmp2", "Ptn", "Ube2c", "Hoxa11", 
              "Hoxa10", "Wnt4", "Erv3", "Gatm", "Gja1", "Lamc1", "Lamb2")
color <- brewer.pal(9, "PuBu")
pdf("20210827-figures/figure2/Figure2-15-MET-E75-gene-featureplot.pdf", width = 7.5, height = 12)
my_spatial_plot(E75.2.bin50, features = MET.gene, color = color)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-MET-E75-1-gene-featureplot.pdf", width = 6.5, height = 11)
my_spatial_plot(E75.1.bin50, features = MET.gene, color = color)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-MET-E85-gene-featureplot.pdf", width = 6, height = 10)
my_spatial_plot(B6.bin50, features = MET.gene, color = color)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-MET-E85-C1-gene-featureplot.pdf", width = 6, height = 10)
my_spatial_plot(C1.bin50, features = MET.gene, color = color)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-MET-E65-gene-featureplot.pdf", width = 6, height = 10)
my_spatial_plot(E65.2.bin50, features = MET.gene, color = color)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-MET-E65-1-gene-featureplot.pdf", width = 5.5, height = 9)
my_spatial_plot(E65.1.bin50, features = MET.gene, color = color)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-MET-E95-gene-featureplot.pdf", width = 9, height = 14)
my_spatial_plot(E95.bin50, features = MET.gene, color = color)
dev.off()

#decidual每个cluster的基因的表达，各个时间点
decidual.gene <- c("Lum", "Dcn", "Postn", "Sfrp4", "Igfbp3", "Hsd11b2", "Top2a", "Ube2c", "Mki67", "Hist1h1b", "Ptn", "Bmp2", "Hmgn5",
                    "Gatm", "Anxa8", "Cxcl14", "S100a8", "Angpt2", "Psca", "Prl8a2", "Erv3", "Htra1", "Rxfp1", "Wnt4", "Hand2", "Twist1", "Zeb1")
pdf("20210827-figures/figure2/Figure2-15-decidual-E75-gene-featureplot.pdf", width = 12.5, height = 10)
my_spatial_plot(E75.2.bin50, features = decidual.gene, color = color, ncol = 5)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-decidual-E75-1-gene-featureplot.pdf", width = 10, height = 9)
my_spatial_plot(E75.1.bin50, features = decidual.gene, color = color, ncol = 5)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-decidual-E85-gene-featureplot.pdf", width = 11, height = 8)
my_spatial_plot(B6.bin50, features = decidual.gene, color = color, ncol = 5)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-decidual-E85-C1-gene-featureplot.pdf", width = 10, height = 8)
my_spatial_plot(C1.bin50, features = decidual.gene, color = color, ncol = 5)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-decidual-E65-gene-featureplot.pdf", width = 10, height = 8)
my_spatial_plot(E65.2.bin50, features = decidual.gene, color = color, ncol = 5)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-decidual-E65-1-gene-featureplot.pdf", width = 9.5, height = 7)
my_spatial_plot(E65.1.bin50, features = decidual.gene, color = color, ncol = 5)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-decidual-E95-gene-featureplot.pdf", width = 15, height = 12)
my_spatial_plot(E95.bin50, features = decidual.gene, color = color, ncol = 5)
dev.off()
pdf("20210827-figures/figure2/Figure2-15-decidual-UMAP-featureplot.pdf", width = 25, height = 30)
FeaturePlot(UE.decidual, features = decidual.gene, pt.size = 0.5, ncol = 5) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()
#画一个局部图
#E7.5-2
E75.2.bin50.celltype.df <- read.table("20210818-spatial/E75-2-bin50-spatial-remove-celltype-df.csv", sep = ",", header = T)
#E75.2.bin50.celltype.MET.df <- E75.2.bin50.celltype.df[E75.2.bin50.celltype.df$cell_type %in% 
#                                                         c("Lum-D", "Sfrp4-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "Top2a-D"), ]

pdf("20210827-figures/figure2-spatial/Figure2-01-E75-2-bin50-MET-rect-label.pdf", width = 6.2, height = 5)
ggplot() + 
  #geom_point(E75.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) +
  geom_point(E75.2.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) + 
  scale_color_manual(values = E75.2.color) + 
  theme_classic() + 
  geom_hline(aes(yintercept = 16000)) +
  geom_hline(aes(yintercept = 17800)) + 
  geom_vline(aes(xintercept = 17000)) +
  geom_vline(aes(xintercept = 15500)) + 
  theme(axis.text = element_blank())
dev.off()

E75.2.bin50.celltype.df.annotation <- E75.2.bin50.celltype.df[E75.2.bin50.celltype.df$x >= 15500 & E75.2.bin50.celltype.df$x <= 17000 & 
                                                                E75.2.bin50.celltype.df$y >= 16000 & E75.2.bin50.celltype.df$y <= 17800, ]
E75.2.ann.color <- E75.2.color[match(intersect(levels(E75.2.bin50.celltype.df$cell_type), unique(E75.2.bin50.celltype.df.annotation$cell_type)), 
                                     levels(E75.2.bin50.celltype.df.annotation$cell_type))]
pdf("20210827-figures/figure2-spatial/Figure2-01-E75-2-bin50-MET-rect-annotation.pdf", width = 3.4, height = 3.6)
ggplot() + 
  geom_point(E75.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 2) +
  #geom_point(E75.2.bin50.celltype.df.annotation.MET, mapping = aes(x = x, y = y, color = cell_type), size = 2) +
  scale_color_manual(values = E75.2.ann.color) + 
  theme_classic() + 
  theme(legend.position = 'none')
dev.off()

#E7.5-2
#E75.1.bin50.celltype.df <- read.table("20210818-spatial/E75-1-bin50-spatial-remove-celltype-df.csv", sep = ",", header = T)
#E75.1.bin50.celltype.MET.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% 
#                                                         c("Lum-D", "Sfrp4-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "Top2a-D"), ]

pdf("20210827-figures/figure2-spatial/Figure2-01-E75-1-bin50-MET-rect.pdf", width = 6.2, height = 5)
ggplot() + 
  #geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) +
  geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) + 
  scale_color_manual(values = E75.1.color) + 
  theme_classic() + 
  geom_hline(aes(yintercept = 6900)) +
  geom_hline(aes(yintercept = 8700)) + 
  geom_vline(aes(xintercept = 6500)) +
  geom_vline(aes(xintercept = 8000)) + 
  theme(axis.text = element_blank())
dev.off()

E75.1.bin50.celltype.df.annotation <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$x >= 6500 & E75.1.bin50.celltype.df$x <= 8000 & 
                                                                E75.1.bin50.celltype.df$y >= 6900 & E75.1.bin50.celltype.df$y <= 8700, ]
E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin50.celltype.df$cell_type), unique(E75.1.bin50.celltype.df.annotation$cell_type)), 
                                     levels(E75.1.bin50.celltype.df.annotation$cell_type))]
pdf("20210827-figures/figure2-spatial/Figure2-01-E75-1-bin50-MET-rect-annotation.pdf", width = 3.4, height = 3.6)
ggplot() + 
  geom_point(E75.1.bin50.celltype.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 2) +
  #geom_point(E75.1.bin50.celltype.df.annotation.MET, mapping = aes(x = x, y = y, color = cell_type), size = 2) +
  scale_color_manual(values = E75.1.ann.color) + 
  theme_classic() + 
  theme(legend.position = 'none')
dev.off()

#E6.5-2
#E65.2.bin50.celltype.MET.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% 
#                                                         c("Lum-D", "Sfrp4-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "Top2a-D"), ]
pdf("20210827-figures/figure2-spatial/Figure2-01-E65-2-bin50-MET-rect.pdf", width = 3.7, height = 3.3)
ggplot() + 
  geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = E65.2.color) + 
  theme_classic() + 
  geom_hline(aes(yintercept = 11200)) +
  geom_hline(aes(yintercept = 13000)) + 
  geom_vline(aes(xintercept = 6300)) +
  geom_vline(aes(xintercept = 7800)) + 
  theme(axis.text = element_blank()) +
  theme(legend.position = "none")
dev.off()
E65.2.bin50.celltype.df.annotation <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$y <= 13000 & E65.2.bin50.celltype.df$y >= 11200 & 
                                                                E65.2.bin50.celltype.df$x <= 7800 & E65.2.bin50.celltype.df$x >= 6300, ]
E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin50.celltype.df$cell_type), unique(E65.2.bin50.celltype.df.annotation$cell_type)), 
                                     levels(E65.2.bin50.celltype.df.annotation$cell_type))]
pdf("20210827-figures/figure2-spatial/Figure2-01-E65-2-bin50-MET-rect-annotation.pdf", width = 3.4, height = 3.6)
ggplot() + 
  geom_point(E65.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 2) +
  scale_color_manual(values = E65.2.ann.color) + 
  theme_classic() +
  theme(legend.position = "none")
dev.off()

#E8.5-B6
#E85.B6.bin50.celltype.df <- read.table("20210818-spatial/E85-B6-bin50-spatial-remove-celltype-df.csv", sep = ",", header = T)
#E85.B6.bin50.celltype.MET.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% 
#                                                         c("Lum-D", "Sfrp4-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "Top2a-D"), ]
pdf("20210827-figures/figure2-spatial/Figure2-01-E85-B6-bin50-MET-rect-label.pdf", width = 5, height = 3.9)
ggplot() + 
  geom_point(B6.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = B6.color) + 
  theme_classic() + 
  geom_hline(aes(yintercept = -7500)) +
  geom_hline(aes(yintercept = -9300)) + 
  geom_vline(aes(xintercept = 4200)) +
  geom_vline(aes(xintercept = 5700)) + 
  theme(axis.text = element_blank())
dev.off()
E85.B6.bin50.celltype.df.annotation <- B6.bin50.celltype.df[B6.bin50.celltype.df$y >= -9300 & B6.bin50.celltype.df$y <= -7500 & 
                                                              B6.bin50.celltype.df$x <= 5700 & B6.bin50.celltype.df$x >= 4200, ]
B6.color.ann.color <- B6.color[match(intersect(levels(B6.bin50.celltype.df$cell_type), unique(E85.B6.bin50.celltype.df.annotation$cell_type)), 
                                     levels(E85.B6.bin50.celltype.df.annotation$cell_type))]

pdf("20210827-figures/figure2-spatial/Figure2-01-E85-B6-bin50-MET-rect-annotation.pdf", width = 3.4, height = 3.4)
ggplot() + 
  geom_point(E85.B6.bin50.celltype.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 2) +
  scale_color_manual(values = B6.color.ann.color) + 
  theme_classic() +
  theme(legend.position = "none")
dev.off()

#E95
#E95.bin50.celltype.MET.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% 
#                                                           c("Lum-D", "Sfrp4-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "Top2a-D"), ]
pdf("20210827-figures/figure2-spatial/Figure2-01-E95-bin50-MET-rect-label.pdf", width = 9, height = 7)
ggplot() + 
  geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = E95.color) + 
  theme_classic() + 
  geom_hline(aes(yintercept = 17700)) +
  geom_hline(aes(yintercept = 19200)) + 
  geom_vline(aes(xintercept = 23100)) +
  geom_vline(aes(xintercept = 24900)) + 
  theme(axis.text = element_blank())
dev.off()
E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df[E95.bin50.celltype.df$y >= 17700 & E95.bin50.celltype.df$y <= 19200 & 
                                                                  E95.bin50.celltype.df$x <= 24900 & E95.bin50.celltype.df$x >= 23100, ]
E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df.annotation[E95.bin50.celltype.df.annotation$cell_type != "Top2a-D", ]
E95.color.ann.color <- E95.color[match(intersect(levels(E95.bin50.celltype.df$cell_type), unique(E95.bin50.celltype.df.annotation$cell_type)), 
                                      levels(E95.bin50.celltype.df.annotation$cell_type))]
pdf("20210827-figures/figure2-spatial/Figure2-01-E95-bin50-MET-rect-annotation.pdf", width = 3.6, height = 3.2)
ggplot() + 
  geom_point(E95.bin50.celltype.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 2) +
  scale_color_manual(values = E95.color.ann.color) + 
  theme_classic() +
  theme(legend.position = "none")
dev.off()

#average的cluster选一个基因来做featureplot
avg.gene <- c("Col1a2", "Lamb1", "Rpl6", "Ube2c", "Tdo2", "Angpt2", "Hand2", "Plvap")
color <- brewer.pal(9, "PuBu")
pdf("20210827-figures/figure2/Figure2-15-AVG-E75-gene-featureplot.pdf", width = 7.5, height = 6)
my_spatial_plot(E75.2.bin50, features = avg.gene, color = color)
dev.off()

all.slides.def <- read.table("20210818-spatial/All-slides-bin50-spatial-remove-celltype-df.csv", sep = ",", header = T, check.names = F)

#in situ genes
library("patchwork")
source("my_spatial_plot.R")
color <- brewer.pal(9, "PuBu")
load("20210818-spatial/E85-B6-bin50-spatial-phago-map-remove.RData")
in.situ.gene <- c("Dcn", "Des", "Hand2", "Ptprc", "Angpt2", "S100a8", "Jcad", "Runx1", "Smoc2", 
                  "C1qa", "C1qb", "Ccl4", "Cxcl12", "Cd74", "Cxcl2", "Prf1", "Gzma", "Gzmb", "Il2rb", "Tnfrsf9", "Nkg7")
pdf("20220622-Comments/10.Adding/04.In-situ-gene-E85-spatial.pdf", width = 7, height = 17)
my_spatial_plot(B6.bin50, features = in.situ.gene, color = color)
dev.off()


avg.gene <- c("Col1a2", "Lamb1", "Rpl6", "Ube2c", "Tdo2", "Angpt2", "Hand2", "Plvap")

pdf("20210827-figures/figure2/Figure2-15-AVG-E75-gene-featureplot.pdf", width = 7.5, height = 6)
my_spatial_plot(E75.2.bin50, features = avg.gene, color = color)
dev.off()








