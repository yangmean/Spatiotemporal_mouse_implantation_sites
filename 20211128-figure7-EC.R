load("../20210201-integrate/UE-EC-harmony-pc30-rename.RData")
UE.EC <- RenameIdents(UE.EC, "EC-2" = "Cycling EC", "EC-0" = "Angiogeneic EC", "EC-4" = "Arterial EC", "EC-1" = "Venous EC1", "EC-3" = "Venous EC2", "EC-5" = "Lympatic EC")
UE.EC.markers <- FindAllMarkers(UE.EC, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
write.table(UE.EC.markers, file = "20210827-figures/files/UE-EC-markers.csv", sep = ",", row.names = F)
UE.EC.markers.top5 <- UE.EC.markers %>% group_by(cluster) %>% top_n(n = 5, avg_logFC)
pdf("20210827-figures/figure5/figure5-EC-top5-markers-dotplot.pdf", width = 14, height = 4)
DotPlot(UE.EC, features = UE.EC.markers.top5$gene, cols = "Spectral") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()
source("../../../01.scripts/StackedVlnPlot.R")
UE.EC$active.ident <- UE.EC@active.ident
pdf("20210827-figures/figure5/figure5-EC-top5-markers-VlnPlot.pdf", width = 8, height = 30)
StackedVlnPlot(UE.EC, features = UE.EC.markers.top5$gene, cols = EC.color[c(3,1,5,2,4,6)])
dev.off()
enrich.GO.function <- function(geneset){
  enrich.GO <- enrichGO(gene = geneset,
                        OrgDb = "org.Mm.eg.db",
                        keyType = "SYMBOL",
                        ont = "BP",
                        pAdjustMethod = "BH")
#  enrich.GO <- simplify(enrich.GO, cutoff = 0.4, by = "p.adjust", select_fun = min)
#  enrich.GO <- enrich.GO@result[enrich.GO@result$pvalue<0.01 & enrich.GO@result$qvalue<0.05, ]
  return(enrich.GO)
}
EC.GO.df <- data.frame()
for (c in levels(UE.EC.markers$cluster)) {
  print(c)
  geneset <- UE.EC.markers[UE.EC.markers$cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  #  enrich.GO <- enrich.GO[1:10, ]
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    gene = enrich.GO$geneID)
  EC.GO.df <- rbind(EC.GO.df, tmp)
}


write.table(EC.GO.df, file = "20210827-figures/files/EC-cells-cluster-GO.csv", sep = ",", row.names = F)
EC.GO.df.top5 <- EC.GO.df %>% group_by(cluster) %>% top_n(n = 5, wt = -log10(pvalue))
EC.GO.df.top5$cluster <- factor(EC.GO.df.top5$cluster, levels = c("Cycling EC","Angiogeneic EC","Arterial EC","Venous EC1","Venous EC2","Lympatic EC"))
EC.GO.df.top5 <- EC.GO.df.top5[rev(order(EC.GO.df.top5$cluster)), ]
EC.GO.df.top5$GO_Term <- factor(EC.GO.df.top5$GO_Term, levels = unique(EC.GO.df.top5$GO_Term))
pdf("20210827-figures/figure7/Figure7-EC-cells-cluster-top5-GO.pdf", width = 8, height = 10)
ggplot() + 
  geom_point(EC.GO.df.top5, mapping = aes(x = cluster, y = GO_Term, color = -log10(qvalue), size = -log10(qvalue))) + 
  scale_color_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 15) + 
  scale_size(range = c(5, 10)) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank()) + 
  labs(title = "EC cells Top5 GO Terms")
dev.off()

UE.EC.count <- data.frame(table(UE.EC$time, UE.EC@active.ident), check.names = F)
pdf("20210827-figures/figure7/Figure7-EC-count-percentage-barplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(UE.EC.count, mapping = aes(x = Var1, y = Freq, fill = Var2), stat = "identity", position = "fill") + 
  scale_fill_manual(values = EC.color[c(3,1,5,2,4,6)]) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", color = "black"))
dev.off()
UE.harmony.pc30$tdTomato <- UE.harmony.pc30@assays$RNA@counts["tdTomato", ]>0
UE.harmony.pc30.tdTomato.neg <- subset(UE.harmony.pc30, subset = tdTomato == FALSE)
total.count <- table(UE.harmony.pc30.tdTomato.neg$time)
total.count <- total.count[c(2:6,1)]
UE.EC.count.table <- table(UE.EC$time, UE.EC@active.ident)
UE.EC.count.table <- data.frame(apply(UE.EC.count.table, 2, function(x){x/total.count * 100}), check.names = F)
UE.EC.count.table$time <- factor(rownames(UE.EC.count.table), levels = c("T55", "T65", "T75", "T85", "T95", "T105"))
UE.EC.count.table <- melt(UE.EC.count.table, id.vars = "time")

pdf("20210827-figures/figure7/figure7-1-EC-proportion-dodge.pdf", width = 10, height = 9)
ggplot() + 
  geom_bar(UE.EC.count.table, mapping = aes(x = variable, y = value, fill = time), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = time.color) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank()) + 
  labs(title = "EC cells proportion")
dev.off()

#计算EC周围最多的几种细胞类型
E65.2.bin50.closed <- CalcSlideClosedCells(E65.2.bin50.celltype.df.sub, bin = 50)
E75.2.bin50.closed <- CalcSlideClosedCells(E75.2.bin50.celltype.df.sub, bin = 50)
E85.B6.bin50.closed <- CalcSlideClosedCells(E85.B6.bin50.celltype.df.sub, bin = 50)
E95.bin50.closed <- CalcSlideClosedCells(E95.bin50.celltype.df.sub, bin = 50)
E65.1.bin50.closed <- CalcSlideClosedCells(E65.1.bin50.celltype.df.sub, bin = 50)
E75.1.bin50.closed <- CalcSlideClosedCells(E75.1.bin50.celltype.df.sub, bin = 50)
E85.C1.bin50.closed <- CalcSlideClosedCells(E85.C1.bin50.celltype.df.sub, bin = 50)

E65.2.bin50.closed.r <- apply(E65.2.bin50.closed, 2, function(x){x/rowSums(E65.2.bin50.closed)})
E65.2.bin50.closed.r <- melt(E65.2.bin50.closed.r)
colnames(E65.2.bin50.closed.r) <- c("source", "bin50.closed", "value")
E65.2.bin50.closed.r$sample <- "E65.2.bin50"

E65.1.bin50.closed.r <- apply(E65.1.bin50.closed, 2, function(x){x/rowSums(E65.1.bin50.closed)})
E65.1.bin50.closed.r <- melt(E65.1.bin50.closed.r)
colnames(E65.1.bin50.closed.r) <- c("source", "bin50.closed", "value")
E65.1.bin50.closed.r$sample <- "E65.1.bin50"

E75.2.bin50.closed.r <- apply(E75.2.bin50.closed, 2, function(x){x/rowSums(E75.2.bin50.closed)})
E75.2.bin50.closed.r <- melt(E75.2.bin50.closed.r)
colnames(E75.2.bin50.closed.r) <- c("source", "bin50.closed", "value")
E75.2.bin50.closed.r$sample <- "E75.2.bin50"

E75.1.bin50.closed.r <- apply(E75.1.bin50.closed, 2, function(x){x/rowSums(E75.1.bin50.closed)})
E75.1.bin50.closed.r <- melt(E75.1.bin50.closed.r)
colnames(E75.1.bin50.closed.r) <- c("source", "bin50.closed", "value")
E75.1.bin50.closed.r$sample <- "E75.1.bin50"

E85.B6.bin50.closed.r <- apply(E85.B6.bin50.closed, 2, function(x){x/rowSums(E85.B6.bin50.closed)})
E85.B6.bin50.closed.r <- melt(E85.B6.bin50.closed.r)
colnames(E85.B6.bin50.closed.r) <- c("source", "bin50.closed", "value")
E85.B6.bin50.closed.r$sample <- "E85.B6.bin50"

E85.C1.bin50.closed.r <- apply(E85.C1.bin50.closed, 2, function(x){x/rowSums(E85.C1.bin50.closed)})
E85.C1.bin50.closed.r <- melt(E85.C1.bin50.closed.r)
colnames(E85.C1.bin50.closed.r) <- c("source", "bin50.closed", "value")
E85.C1.bin50.closed.r$sample <- "E85.C1.bin50"

E95.bin50.closed.r <- apply(E95.bin50.closed, 2, function(x){x/rowSums(E95.bin50.closed)})
E95.bin50.closed.r <- melt(E95.bin50.closed.r)
colnames(E95.bin50.closed.r) <- c("source", "bin50.closed", "value")
E95.bin50.closed.r$sample <- "E95.bin50"

all.bin50.closed <- Reduce(rbind, list(E65.1.bin50.closed.r, E65.2.bin50.closed.r, E75.1.bin50.closed.r, E75.2.bin50.closed.r, E85.C1.bin50.closed.r, E85.B6.bin50.closed.r, E95.bin50.closed.r))
all.bin50.closed <- all.bin50.closed[all.bin50.closed$bin50.closed != "Erythroid", ]
write.table(all.bin50.closed, file = "20210827-figures/files/all-slides-bin50-closed-cells-calc.csv", sep = ",", row.names = F, quote = F)
all.bin50.closed.EC <- all.bin50.closed[all.bin50.closed$source %in% c(EC.levels, "FB_immune_0", "FB_immune_1", "FB_immune_2"), ]
all.bin50.closed.EC$source <- factor(all.bin50.closed.EC$source, levels = c(EC.levels[c(3,1,5,2,4,6)], "FB_immune_0", "FB_immune_1", "FB_immune_2"))
all.bin50.closed.EC$bin50.closed <- factor(all.bin50.closed.EC$bin50.closed, 
                                                levels = c(decidual.levels[1:7], immune.levels[1:9], "FB_immune_0", "FB_immune_1", "FB_immune_2", 
                                                           EC.levels[c(3,1,5,2,4,6)], major.levels[c(4,5,11)], tropho.levels[1:5], "embryo"))
pdf("20210827-figures/figure7/Figure7-01-EC-FB-immune-distance.pdf", width = 25, height = 15)
ggplot() + 
  geom_bar(all.bin50.closed.EC, mapping = aes(x = bin50.closed, y = value, fill = bin50.closed), stat = "identity", width = 1) + 
  facet_grid(sample ~ source, space = "free") +
  scale_fill_manual(values = c(decidual.color[1:7], immune.color[c(1:9)], "#F0746F", "#39B54E", "#6292CC", EC.color[c(3,1,5,2,4,6)], major.color[c(4,5,11)], tropho.color[1:5], "lightgray")) +
  theme_bw() + 
  theme(axis.text.x = element_blank())
dev.off()

#EC与DSC, immune cell, iDSC之间的距离
EC.select.cell.DSC <- c(decidual.levels[c(1:3,5:7)])
all.bin50.closed.EC.DSC <- all.bin50.closed[all.bin50.closed$source %in% c(EC.levels) &
                                            all.bin50.closed$bin50.closed %in% EC.select.cell.DSC, ]
all.bin50.closed.EC.DSC$source <- factor(all.bin50.closed.EC.DSC$source, levels = c(EC.levels[c(3,1,5,2,4,6)]))
all.bin50.closed.EC.DSC$bin50.closed <- factor(all.bin50.closed.EC.DSC$bin50.closed, levels = EC.select.cell)
pdf("20210827-figures/figure7/Figure7-01-EC-DSC-distance.pdf", width = 15, height = 15)
ggplot() + 
  geom_bar(all.bin50.closed.EC.DSC, mapping = aes(x = bin50.closed, y = value, fill = bin50.closed), stat = "identity", width = 1) + 
  facet_grid(sample ~ source, space = "free") +
  scale_fill_manual(values = c(decidual.color[c(1:3,5:7)], immune.color[c(1:7)], "#F0746F", "#39B54E", "#6292CC")) +
  theme_bw() + 
  theme(axis.text.x = element_blank())
dev.off()
EC.select.cell.i <- c(immune.levels[1:7], "FB_immune_0", "FB_immune_1", "FB_immune_2")
all.bin50.closed.EC.i <- all.bin50.closed[all.bin50.closed$source %in% c(EC.levels) &
                                              all.bin50.closed$bin50.closed %in% EC.select.cell.i, ]
all.bin50.closed.EC.i$source <- factor(all.bin50.closed.EC.i$source, levels = c(EC.levels[c(3,1,5,2,4,6)]))
all.bin50.closed.EC.i$bin50.closed <- factor(all.bin50.closed.EC.i$bin50.closed, levels = EC.select.cell.i)
pdf("20210827-figures/figure7/Figure7-01-EC-immune-distance.pdf", width = 15, height = 15)
ggplot() + 
  geom_bar(all.bin50.closed.EC.i, mapping = aes(x = bin50.closed, y = value, fill = bin50.closed), stat = "identity", width = 1) + 
  facet_grid(sample ~ source, space = "free") +
  scale_fill_manual(values = c(immune.color[c(1:7)], "#F0746F", "#39B54E", "#6292CC")) +
  theme_bw() + 
  theme(axis.text.x = element_blank())
dev.off()
########################################看空间中的EC的分布(EC0-5, FB-immune0, S100a8-D, NK)##############################
cells <- c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "EC-5", "S100a8-D", "NK", "FB_immune_0")
#E65.2
E65.2.bin50.celltype.df.sub.EC <- E65.2.bin50.celltype.df.sub[E65.2.bin50.celltype.df.sub$cell_type %in% cells, ]
E65.2.bin50.celltype.df.sub.EC$cell_type <- factor(E65.2.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E65-2-bin50-EC-decidual-immune.pdf", width = 3.3, height = 4.2)
ggplot() + 
  geom_point(E65.2.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E65.2.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color, decidual.color[6], immune.color[7], "#F0746F")) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

#E75.1
E75.1.bin50.celltype.df.sub.EC <- E75.1.bin50.celltype.df.sub[E75.1.bin50.celltype.df.sub$cell_type %in% cells, ]
E75.1.bin50.celltype.df.sub.EC$cell_type <- factor(E75.1.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E75-1-bin50-EC-decidual-immune.pdf", width = 4.5, height = 5.2)
ggplot() + 
  geom_point(E75.1.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E75.1.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color[c(1,3:5)], decidual.color[6], immune.color[7], "#F0746F")) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

#E75.2
E75.2.bin50.celltype.df.sub.EC <- E75.2.bin50.celltype.df.sub[E75.2.bin50.celltype.df.sub$cell_type %in% cells, ]
E75.2.bin50.celltype.df.sub.EC$cell_type <- factor(E75.2.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E75-2-bin50-EC-decidual-immune.pdf", width = 5, height = 5.8)
ggplot() + 
  geom_point(E75.2.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E75.2.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color[1:5], decidual.color[6], immune.color[7], "#F0746F")) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

#E85.2
E85.B6.bin50.celltype.df.sub.EC <- E85.B6.bin50.celltype.df.sub[E85.B6.bin50.celltype.df.sub$cell_type %in% cells, ]
E85.B6.bin50.celltype.df.sub.EC$cell_type <- factor(E85.B6.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E85-B6-bin50-EC-decidual-immune.pdf", width = 4, height = 4.5)
ggplot() + 
  geom_point(E85.B6.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E85.B6.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color[1:4], decidual.color[6], immune.color[7], "#F0746F")) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

#E95
E95.bin50.celltype.df.sub.EC <- E95.bin50.celltype.df.sub[E95.bin50.celltype.df.sub$cell_type %in% cells, ]
E95.bin50.celltype.df.sub.EC$cell_type <- factor(E95.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E95-bin50-EC-decidual-immune.pdf", width = 9, height = 7)
ggplot() + 
  geom_point(E95.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E95.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color[1:5], decidual.color[6], immune.color[7], "#F0746F")) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

cells <- c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "EC-5")
#E65.2
E65.2.bin50.celltype.df.sub.EC <- E65.2.bin50.celltype.df.sub[E65.2.bin50.celltype.df.sub$cell_type %in% cells, ]
E65.2.bin50.celltype.df.sub.EC$cell_type <- factor(E65.2.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E65-2-bin50-EC.pdf", width = 3.3, height = 4.2)
ggplot() + 
  geom_point(E65.2.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E65.2.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color)) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

#E75.1
E75.1.bin50.celltype.df.sub.EC <- E75.1.bin50.celltype.df.sub[E75.1.bin50.celltype.df.sub$cell_type %in% cells, ]
E75.1.bin50.celltype.df.sub.EC$cell_type <- factor(E75.1.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E75-1-bin50-EC.pdf", width = 4.4, height = 5.1)
ggplot() + 
  geom_point(E75.1.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E75.1.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color[c(1,3,4,5)])) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

#E75.2
E75.2.bin50.celltype.df.sub.EC <- E75.2.bin50.celltype.df.sub[E75.2.bin50.celltype.df.sub$cell_type %in% cells, ]
E75.2.bin50.celltype.df.sub.EC$cell_type <- factor(E75.2.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E75-2-bin50-EC.pdf", width = 5, height = 5.8)
ggplot() + 
  geom_point(E75.2.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E75.2.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color[1:5])) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

#E85.2
E85.B6.bin50.celltype.df.sub.EC <- E85.B6.bin50.celltype.df.sub[E85.B6.bin50.celltype.df.sub$cell_type %in% cells, ]
E85.B6.bin50.celltype.df.sub.EC$cell_type <- factor(E85.B6.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E85-B6-bin50-EC.pdf", width = 4, height = 4.5)
ggplot() + 
  geom_point(E85.B6.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E85.B6.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color[1:4])) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

#E95
E95.bin50.celltype.df.sub.EC <- E95.bin50.celltype.df.sub[E95.bin50.celltype.df.sub$cell_type %in% cells, ]
E95.bin50.celltype.df.sub.EC$cell_type <- factor(E95.bin50.celltype.df.sub.EC$cell_type, levels = cells)
pdf("20210827-figures/figure7/Figure7-03-E95-bin50-EC.pdf", width = 9, height = 7)
ggplot() + 
  geom_point(E95.bin50.celltype.df.sub, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(E95.bin50.celltype.df.sub.EC, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = c(EC.color[1:5])) + 
  theme_classic() +
  theme(axis.text = element_blank(), legend.position = "none")
dev.off()

#看EC通路的基因表达
EC.interactions <- c("VEGF", "SEMA3", "CALCR", "CHEMERIN", "ANGPT", "HGF", "TRAIL", "APELIN")
UE.harmony.cellchat.EC <- subsetCellChat(UE.harmony.cellchat, idents.use = c(decidual.levels[c(1:5,7,6)], immune.levels[1:7], "FB-immune0", "FB-immune1", "FB-immune2", EC.levels[1:5]))
color2 <- c(decidual.color[c(1:5,7,6)], immune.color[c(4,5,3,2,7,6)], EC.color[1:5], rep(immune.color[10], 3))
pdf("20210827-figures/figure7/Figure7-04-EC-pathways.pdf", width = 15, height = 20)
plotGeneExpression(UE.harmony.cellchat.EC, signaling = EC.interactions, color.use = color2)
dev.off()

pdf("20210827-figures/figure7/Figure7-04-EC-pathways-1.pdf", width = 15, height = 5)
plotGeneExpression(UE.harmony.cellchat.EC, signaling = c("CD137", "CX3C"), color.use = color2)
dev.off()

#S100a8-D, NK, FB-immune给EC的signals
#与angiogenesis相关的pathways
df.net.LR <- subsetCommunication(UE.harmony.cellchat)
EC.interactions <- c("VEGF", "SEMA3", "CALCR", "CHEMERIN", "ANGPT", "HGF", "TRAIL", "APELIN")
EC.pairLR <- extractEnrichedLR(UE.harmony.cellchat, signaling = EC.interactions, geneLR.return = FALSE)

df.net.LR.interaction.EC <- df.net.LR[df.net.LR$interaction_name %in% EC.pairLR$interaction_name & 
                                        df.net.LR$source %in% c("FB-immune", "S100a8-D", "NK") & 
                                        df.net.LR$target %in% c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4"), c(1, 2, 5, 8)]
df.net.LR.interaction.EC$cell_LR <- paste(df.net.LR.interaction.EC$source, df.net.LR.interaction.EC$target, sep = "->")
df.net.LR.interaction.EC <- dcast(df.net.LR.interaction.EC[, 3:5], cell_LR ~ interaction_name_2, value.var = "prob")
rownames(df.net.LR.interaction.EC) <- df.net.LR.interaction.EC$cell_LR
df.net.LR.interaction.EC <- df.net.LR.interaction.EC[, -1]
df.net.LR.interaction.EC[is.na(df.net.LR.interaction.EC)] <- 0
color <- brewer.pal(11, "RdBu")
color <- color[c(1:5, 7:11)]
pdf("20210827-figures/figure7/figure7-02-EC-interactions-heatmap.pdf", width = 10, height = 5)
pheatmap(df.net.LR.interaction.EC, scale = "column", cluster_rows = F, cluster_cols = T, 
         color = colorRampPalette(colors = rev(color))(100))
dev.off()

#EC incoming pathways
UE.harmony.cellchat.EC.incoming <- subsetCellChat(UE.harmony.cellchat, idents.use = c("FB-immune", "S100a8-D", "NK", "EC-0", "EC-1", "EC-2", "EC-3", "EC-4"))
color2 <- c("#FDB462", "#CAB2D6", "#6A3D9A", EC.color[1:5])
pdf("20210827-figures/figure7/figure7-02-EC-incoming.pdf", width = 8, height = 15)
plotGeneExpression(UE.harmony.cellchat.EC.incoming, signaling = EC.interactions, color.use = color2)
dev.off()

#monocle
library("monocle")
load("../20210201-integrate/UE-endothelial-monocle.RData")
pdf("20210827-figures/figure7/Figure7-EC-cells-cluster-monocle-color-clusters.pdf", width = 10, height = 6)
plot_cell_trajectory(EC.cds, color_by = "seurat_clusters", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_manual(values = EC.color)
dev.off()
time.color <- c(brewer.pal(n = 11, "RdBu"))
time.color <- time.color[c(2:4, 8:11)]
pdf("20210827-figures/figure7/Figure7-EC-cells-cluster-monocle-color-time.pdf", width = 10, height = 6)
plot_cell_trajectory(EC.cds, color_by = "time", show_state_number = F, cell_size = 1, show_branch_points = F) + 
  scale_color_manual(values = time.color)
dev.off()
pdf("20210827-figures/figure7/Figure7-EC-cells-cluster-monocle-color-pseudotime.pdf", width = 10, height = 6)
plot_cell_trajectory(EC.cds, color_by = "Pseudotime", show_state_number = F, cell_size = 1, show_branch_points = T) + 
  scale_color_gradient2(low = "blue", high = "red", mid = "orange", midpoint = 15)
dev.off()
#branch2的差异基因分析
EC.BEAM_res <- BEAM(EC.cds, branch_point = 2, cores = 20)
EC.BEAM_res.order <- EC.BEAM_res[EC.BEAM_res$use_for_ordering == TRUE & EC.BEAM_res$num_cells_expressed>10 & EC.BEAM_res$qval<0.00001, ]
EC.branch2.heatmap <- plot_genes_branched_heatmap(EC.cds[row.names(EC.BEAM_res.order),],
                            branch_point = 2,
                            num_clusters = 6,
                            cores = 5,
                            show_rownames = F, 
                            hmcols = c(colorRampPalette(colors = rev(brewer.pal(11, "RdBu")))(60)), 
                            return_heatmap = T)

EC.branch2.heatmap.cluster <- data.frame(cluster = EC.branch2.heatmap$annotation_row, 
                                         gene = rownames(EC.branch2.heatmap$annotation_row))
EC.branch2.heatmap.cluster$Cluster <- gsub("^1$", "C6", EC.branch2.heatmap.cluster$Cluster)
EC.branch2.heatmap.cluster$Cluster <- gsub("^2$", "C2", EC.branch2.heatmap.cluster$Cluster)
EC.branch2.heatmap.cluster$Cluster <- gsub("^3$", "C5", EC.branch2.heatmap.cluster$Cluster)
EC.branch2.heatmap.cluster$Cluster <- gsub("^4$", "C1", EC.branch2.heatmap.cluster$Cluster)
EC.branch2.heatmap.cluster$Cluster <- gsub("^5$", "C3", EC.branch2.heatmap.cluster$Cluster)
EC.branch2.heatmap.cluster$Cluster <- gsub("^6$", "C4", EC.branch2.heatmap.cluster$Cluster)

EC.pseudotime.GO.df <- data.frame()
for (c in unique(EC.branch2.heatmap.cluster$Cluster)) {
  print(c)
  geneset <- EC.branch2.heatmap.cluster[EC.branch2.heatmap.cluster$Cluster == c, ]$gene
  enrich.GO <- enrich.GO.function(geneset)
  #  enrich.GO <- enrich.GO[1:10, ]
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    gene = enrich.GO$geneID)
  EC.pseudotime.GO.df <- rbind(EC.pseudotime.GO.df, tmp)
}
write.table(EC.pseudotime.GO.df, file = "20210827-figures/files/EC-pseudotime-DEG-GOs.csv", sep = ",", quote = F)

###20230126 DSC and iDSC to ECs
#decidual cells to immune cells
load("../20210201-integrate/UE-harmony-cellchat-subname-FB-subclusters-results-20211007.RData")
UE.harmony.cellchat
DSC.EC.pathway <- c("CALCR", "ANGPT", "ANGPTL", "CHEMERIN", "SEMA3", "VEGF")
UE.harmony.cellchat.df <- subsetCommunication(UE.harmony.cellchat, slot.name = "net")
UE.harmony.cellchat.df.DSC.EC <- UE.harmony.cellchat.df[UE.harmony.cellchat.df$source %in% c(decidual.cell[c(1:3,5:6)], "FB-immune0", "FB-immune1", "FB-immune2") & 
                                                              UE.harmony.cellchat.df$target %in% EC.cell[c(1:5)] & 
                                                              UE.harmony.cellchat.df$pathway_name %in% DSC.EC.pathway, ]
UE.harmony.cellchat.df.DSC.EC$source_target <- paste(UE.harmony.cellchat.df.DSC.EC$source, UE.harmony.cellchat.df.DSC.EC$target, sep = "_")
UE.harmony.cellchat.df.DSC.EC.d <- dcast(UE.harmony.cellchat.df.DSC.EC, source_target ~ interaction_name_2, value.var = "prob")
rownames(UE.harmony.cellchat.df.DSC.EC.d) <- UE.harmony.cellchat.df.DSC.EC.d$source_target
UE.harmony.cellchat.df.DSC.EC.d <- UE.harmony.cellchat.df.DSC.EC.d[, -1]
UE.harmony.cellchat.df.DSC.EC.d[is.na(UE.harmony.cellchat.df.DSC.EC.d)] <- 0
UE.harmony.cellchat.df.DSC.EC.d <- apply(UE.harmony.cellchat.df.DSC.EC.d, 2, function(x){(x-min(x))/(max(x)-min(x))})
UE.harmony.cellchat.df.DSC.EC.d <- UE.harmony.cellchat.df.DSC.EC.d[c(3,1,5,2,4,8,6,10,7,9,13,11,15,12,14,23,21,25,22,24,33,31,35,32,34,38,36,40,37,39,18,16,20,17,19,28,26,30,27,29), 
                                                                   c(13,1:12,15, 16,17,18,20,19,21,28,23,25,22,24,26,27, 29:35)]
pdf("20210827-figures/figure6/Figure6-01-DSC-to-EC-cellchat-heatmap-20230126.pdf", width = 10, height = 10)
pheatmap(UE.harmony.cellchat.df.DSC.EC.d, scale = "none", cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(colors = c("white", "#20B2AA"))(100))
dev.off()

load("../20210201-integrate/UE-EC-harmony-pc30-rename.RData")
UE.EC <- RenameIdents(UE.EC, "EC-2" = "Cycling EC", "EC-0" = "Angiogeneic EC", "EC-4" = "Arterial EC", "EC-1" = "Venous EC1", "EC-3" = "Venous EC2", "EC-5" = "Lympatic EC")
UE.EC.markers <- read.table("20210827-figures/files/UE-EC-markers.csv", sep = ",", check.names = F, header = T)
UE.EC.markers.top20 <- UE.EC.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
pdf("20220622-Comments/10.Adding/03.UE-EC-DEGs-heatmap.pdf", width = 10, height = 10)
DoHeatmap(UE.EC, features = as.character(UE.EC.markers.top20$gene)) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()

EC.genes <- c("Flt1", "Kdr", "Dll4", "Notch1")
pdf("20220622-Comments/10.Adding/06.EC-featureplot.pdf", width = 10, height = 10)
FeaturePlot(UE.EC, features = EC.genes, pt.size = 0.3) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()




