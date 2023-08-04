decidual.cell <- c("Lum-D", "Sfrp4-D", "Top2a-D", "Ptn-D", "Gatm-D", "S100a8-D", "Prl8a2-D")
immune.cell <- c("Mono", "Mac", "Pro-Mac", "DC-1", "DC-2", "Neutrophil", "NK", "B", "T", "FB-immune")
EC.cell <- c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "EC-5")

CalcSlideClosedCells <- function(df, bin = 50){
  df$x <- abs(as.numeric(df$x))
  df$y <- abs(as.numeric(df$y))
  celltype <- as.character(unique(df$cell_type))
  result.df <- data.frame(target = unique(df$cell_type), check.names = F)
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
    names <- c(names1, names2, names3, names4, names5, names6, names7, names8)
    closed.df <- df[names, ]#包含重复的点
    closed.df <- na.omit(closed.df)
    closed.df.tb <- data.frame(value = table(closed.df$cell_type), check.names = F)
    colnames(closed.df.tb) <- c("target", celltype[c])
    result.df <- merge(result.df, closed.df.tb, all = T)
    result.df[is.na(result.df)] <- 0
  }
  rownames(result.df) <- result.df$target
  result.df <- t(result.df[, -1])
  return(result.df)
}
##############################################看bin20下免疫细胞周围的细胞类型##########################################################
#E65.2.bin20.celltype.df.sub这个数据在FB-immune这个脚本中
E65.2.bin20.closed <- CalcSlideClosedCells(E65.2.bin20.celltype.df.sub, bin = 20)
E75.2.bin20.closed <- CalcSlideClosedCells(E75.2.bin20.celltype.df.sub, bin = 20)
E85.B6.bin20.closed <- CalcSlideClosedCells(E85.B6.bin20.celltype.df.sub, bin = 20)
E95.bin20.closed <- CalcSlideClosedCells(E95.bin20.celltype.df.sub, bin = 20)
E65.1.bin20.closed <- CalcSlideClosedCells(E65.1.bin20.celltype.df.sub, bin = 20)
E75.1.bin20.closed <- CalcSlideClosedCells(E75.1.bin20.celltype.df.sub, bin = 20)
E85.C1.bin20.closed <- CalcSlideClosedCells(E85.C1.bin20.celltype.df.sub, bin = 20)

E65.2.bin20.closed.r <- apply(E65.2.bin20.closed, 2, function(x){x/rowSums(E65.2.bin20.closed)})
E65.2.bin20.closed.r <- melt(E65.2.bin20.closed.r)
colnames(E65.2.bin20.closed.r) <- c("source", "bin20.closed", "value")
#E65.2.bin20.closed.r <- E65.2.bin20.closed.r %>% group_by(source) %>% top_n(n = 5, value)
E65.2.bin20.closed.r$sample <- "E65.2.bin20"

E65.1.bin20.closed.r <- apply(E65.1.bin20.closed, 2, function(x){x/rowSums(E65.1.bin20.closed)})
E65.1.bin20.closed.r <- melt(E65.1.bin20.closed.r)
colnames(E65.1.bin20.closed.r) <- c("source", "bin20.closed", "value")
#E65.1.bin20.closed.r <- E65.1.bin20.closed.r %>% group_by(source) %>% top_n(n = 5, value)
E65.1.bin20.closed.r$sample <- "E65.1.bin20"

E75.2.bin20.closed.r <- apply(E75.2.bin20.closed, 2, function(x){x/rowSums(E75.2.bin20.closed)})
E75.2.bin20.closed.r <- melt(E75.2.bin20.closed.r)
colnames(E75.2.bin20.closed.r) <- c("source", "bin20.closed", "value")
#E75.2.bin20.closed.r <- E75.2.bin20.closed.r %>% group_by(source) %>% top_n(n = 5, value)
E75.2.bin20.closed.r$sample <- "E75.2.bin20"

E75.1.bin20.closed.r <- apply(E75.1.bin20.closed, 2, function(x){x/rowSums(E75.1.bin20.closed)})
E75.1.bin20.closed.r <- melt(E75.1.bin20.closed.r)
colnames(E75.1.bin20.closed.r) <- c("source", "bin20.closed", "value")
#E75.1.bin20.closed.r <- E75.1.bin20.closed.r %>% group_by(source) %>% top_n(n = 5, value)
E75.1.bin20.closed.r$sample <- "E75.1.bin20"

E85.B6.bin20.closed.r <- apply(E85.B6.bin20.closed, 2, function(x){x/rowSums(E85.B6.bin20.closed)})
E85.B6.bin20.closed.r <- melt(E85.B6.bin20.closed.r)
colnames(E85.B6.bin20.closed.r) <- c("source", "bin20.closed", "value")
#E85.B6.bin20.closed.r <- E85.B6.bin20.closed.r %>% group_by(source) %>% top_n(n = 5, value)
E85.B6.bin20.closed.r$sample <- "E85.B6.bin20"

E85.C1.bin20.closed.r <- apply(E85.C1.bin20.closed, 2, function(x){x/rowSums(E85.C1.bin20.closed)})
E85.C1.bin20.closed.r <- melt(E85.C1.bin20.closed.r)
colnames(E85.C1.bin20.closed.r) <- c("source", "bin20.closed", "value")
#E85.C1.bin20.closed.r <- E85.C1.bin20.closed.r %>% group_by(source) %>% top_n(n = 5, value)
E85.C1.bin20.closed.r$sample <- "E85.C1.bin20"

E95.bin20.closed.r <- apply(E95.bin20.closed, 2, function(x){x/rowSums(E95.bin20.closed)})
E95.bin20.closed.r <- melt(E95.bin20.closed.r)
colnames(E95.bin20.closed.r) <- c("source", "bin20.closed", "value")
#E95.bin20.closed.r <- E95.bin20.closed.r %>% group_by(source) %>% top_n(n = 5, value)
E95.bin20.closed.r$sample <- "E95.bin20"

all.bin20.closed <- Reduce(rbind, list(E65.1.bin20.closed.r, E65.2.bin20.closed.r, E75.1.bin20.closed.r, E75.2.bin20.closed.r, E85.C1.bin20.closed.r, E85.B6.bin20.closed.r, E95.bin20.closed.r))
all.bin20.closed.immune <- all.bin20.closed[all.bin20.closed$source %in% immune.levels[1:7] & 
                                       all.bin20.closed$bin20.closed %in% decidual.levels[c(1:3,5:7)], ]
all.bin20.closed.immune$source <- factor(all.bin20.closed.immune$source, levels = immune.levels[1:7])
all.bin20.closed.immune$bin20.closed <- factor(all.bin20.closed.immune$bin20.closed, levels = decidual.levels[c(1:3,5:7)])
all.bin20.closed.immune[all.bin20.closed.immune==0] <- NA
all.bin20.closed.immune <- all.bin20.closed.immune[all.bin20.closed.immune$sample %in% c("E65.2.bin20", "E75.1.bin20", "E85.B6.bin20", "E95.bin20"), ]
pdf("20210827-figures/figure6/Figure6-01-immune-DSC-closed-bin20.pdf", width = 7, height = 5)
ggplot() + 
  geom_point(all.bin20.closed.immune, mapping = aes(x = source, y = bin20.closed, color = source, size = value), stat = "identity") + 
  scale_color_manual(values = c(immune.color[c(1:7)])) + 
  facet_wrap( ~ sample) + 
  scale_size(range = c(0, 10)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", color = "black"))
dev.off()

#decidual cells to immune cells
load("../20210201-integrate/UE-harmony-cellchat-subname-FB-subclusters-results-20211007.RData")
UE.harmony.cellchat
DSC.immune.pathway <- c("CCL","CXCL","MK","GALECTIN")
UE.harmony.cellchat.df <- subsetCommunication(UE.harmony.cellchat, slot.name = "net")
UE.harmony.cellchat.df.DSC.immune <- UE.harmony.cellchat.df[UE.harmony.cellchat.df$source %in% c(decidual.cell[c(1:3,5:6)], "FB-immune0", "FB-immune1", "FB-immune2") & 
                                                              UE.harmony.cellchat.df$target %in% immune.cell[c(2,3,4,5,6,7)] & 
                                                              UE.harmony.cellchat.df$pathway_name %in% DSC.immune.pathway, ]
UE.harmony.cellchat.df.DSC.immune$source_target <- paste(UE.harmony.cellchat.df.DSC.immune$source, UE.harmony.cellchat.df.DSC.immune$target, sep = "_")
UE.harmony.cellchat.df.DSC.immune.d <- dcast(UE.harmony.cellchat.df.DSC.immune, source_target ~ interaction_name_2, value.var = "prob")
rownames(UE.harmony.cellchat.df.DSC.immune.d) <- UE.harmony.cellchat.df.DSC.immune.d$source_target
UE.harmony.cellchat.df.DSC.immune.d <- UE.harmony.cellchat.df.DSC.immune.d[, -1]
UE.harmony.cellchat.df.DSC.immune.d[is.na(UE.harmony.cellchat.df.DSC.immune.d)] <- 0
UE.harmony.cellchat.df.DSC.immune.d <- apply(UE.harmony.cellchat.df.DSC.immune.d, 2, function(x){(x-min(x))/(max(x)-min(x))})

UE.harmony.cellchat.df.DSC.immune.d <- UE.harmony.cellchat.df.DSC.immune.d[c(3,6,1,2,5,4,9,12,7,8,11,10,15,18,13,14,17,16,27,30,25,26,29,28,39,42,37,38,41,40,45,48,43,44,47,46,21,24,19,20,23,22,33,36,31,32,35,34), 
                                                                           c(20,21,23,19,22,11,12,10,13,14,24,5,2,9,3,4,6,1,8,7,17,15,16,18)]
pdf("20210827-figures/figure6/Figure6-01-DSC-to-immune-cellchat-heatmap-20230126.pdf", width = 10, height = 10)
pheatmap(UE.harmony.cellchat.df.DSC.immune.d, scale = "none", cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(colors = c("white", "#FF8C00"))(100))
dev.off()
netVisual_chord_gene(UE.harmony.cellchat, signaling = c("CXCL"), sources.use = c(2,3,4,6,7,8,10:13,16,17), targets.use = c(10:13,16,17))
netVisual_chord_gene(UE.harmony.cellchat, signaling = c("CCL"), sources.use = c(2,3,4,6,7,8,10:13,16,17), targets.use = c(10:13,16,17))


UE.harmony.cellchat.subset <- subsetCellChat(UE.harmony.cellchat, 
                                             idents.use = c(decidual.levels[c(1:3, 5:6)], immune.levels[c(1:7)], "FB-immune0", "FB-immune1", "FB-immune2"))
netVisual_aggregate(UE.harmony.cellchat.subset, signaling = c("CXCL", "CCL"), layout = "chord")
netVisual_aggregate(UE.harmony.cellchat.subset, signaling = "CCL", layout = "chord")
pdf("20210827-figures/figure6/Figure6-08-CXCL-chord-cell.pdf", width = 10, height = 10)
netVisual_chord_cell(UE.harmony.cellchat.subset, signaling = c("CXCL"), transparency = 0.2,
                     sources.use = c(decidual.levels[c(1:3, 5:6)], "FB-immune0", "FB-immune1", "FB-immune2"), 
                     targets.use = c(immune.levels[c(1:7)]), 
                     color.use = c(decidual.color[c(1:3,5,6)], immune.color[c(4,5,3,2,7,6)], "#F0746F", "#39B54E", "#6292CC"))
dev.off()
pdf("20210827-figures/figure6/Figure6-08-CCL-chord-cell.pdf", width = 10, height = 10)
netVisual_chord_cell(UE.harmony.cellchat.subset, signaling = c("CCL"), transparency = 0.2,
                     sources.use = c(decidual.levels[c(1:3, 5:6)], "FB-immune0", "FB-immune1", "FB-immune2"), 
                     targets.use = c(immune.levels[c(1:7)]), 
                     color.use = c(decidual.color[c(1:3,5,6)], immune.color[c(4,5,3,2,7,6)], "#F0746F", "#39B54E", "#6292CC"))
dev.off()
pdf("20210827-figures/figure6/Figure6-08-MK-chord-cell.pdf", width = 10, height = 10)
netVisual_chord_cell(UE.harmony.cellchat.subset, signaling = c("MK"), transparency = 0.2,
                     sources.use = c(decidual.levels[c(1:3, 5:6)], "FB-immune0", "FB-immune1", "FB-immune2"), 
                     targets.use = c(immune.levels[c(1:7)]), 
                     color.use = c(decidual.color[c(1:3,5,6)], immune.color[c(4,5,3,2,7,6)], "#F0746F", "#39B54E", "#6292CC"))
dev.off()

UE.harmony.cellchat.subset <- subsetCellChat(UE.harmony.cellchat, 
                                             idents.use = c(decidual.levels[c(1:3, 5:7)], immune.levels[c(1:7)], "FB-immune0", "FB-immune1", "FB-immune2"))
pdf("20210827-figures/figure6/Figure6-08-GALECTIN-chord-cell.pdf", width = 10, height = 10)
netVisual_chord_cell(UE.harmony.cellchat.subset, signaling = c("GALECTIN"), transparency = 0.2,
                     sources.use = c(decidual.levels[c(1:3, 5:7)], "FB-immune0", "FB-immune1", "FB-immune2"), 
                     targets.use = c(immune.levels[c(1:7)]), 
                     color.use = c(decidual.color[c(1:3,5,7,6)], immune.color[c(4,5,3,2,7,6)], "#F0746F", "#39B54E", "#6292CC"))
dev.off()
#以每个片子的中间为界，统计免疫细胞在每个片子中上/下细胞比例值
#E65.1
E65.1.bin20.celltype.df.immune <- E65.1.bin20.celltype.df[E65.1.bin20.celltype.df$cell_type %in% immune.cell, ]
E65.1.bin20.celltype.df.immune$cell_type <- factor(E65.1.bin20.celltype.df.immune$cell_type, levels = immune.cell)
median <- (min(E65.1.bin20.celltype.df$y) + max(E65.1.bin20.celltype.df$y))/2
E65.1.bin20.celltype.df.immune$region <- "up"
E65.1.bin20.celltype.df.immune[E65.1.bin20.celltype.df.immune$y<10100, ]$region <- "down" #三分之一分位点
E65.1.bin20.celltype.df.immune.df <- table(E65.1.bin20.celltype.df.immune$region, E65.1.bin20.celltype.df.immune$cell_type)
E65.1.bin20.celltype.df.immune.df <- data.frame(apply(E65.1.bin20.celltype.df.immune.df, 1, function(x) x/colSums(E65.1.bin20.celltype.df.immune.df)), check.names = F)
E65.1.bin20.celltype.df.immune.df$sample <- "E65.1.bin20"
E65.1.bin20.celltype.df.immune.df$celltype <- rownames(E65.1.bin20.celltype.df.immune.df)
#E65.2
E65.2.bin20.celltype.df.immune <- E65.2.bin20.celltype.df[E65.2.bin20.celltype.df$cell_type %in% immune.cell, ]
E65.2.bin20.celltype.df.immune$cell_type <- factor(E65.2.bin20.celltype.df.immune$cell_type, levels = immune.cell)
median <- (min(E65.2.bin20.celltype.df$y) + max(E65.2.bin20.celltype.df$y))/2
E65.2.bin20.celltype.df.immune$region <- "up"
E65.2.bin20.celltype.df.immune[E65.2.bin20.celltype.df.immune$y<13900, ]$region <- "down"
E65.2.bin20.celltype.df.immune.df <- table(E65.2.bin20.celltype.df.immune$region, E65.2.bin20.celltype.df.immune$cell_type)
E65.2.bin20.celltype.df.immune.df <- data.frame(apply(E65.2.bin20.celltype.df.immune.df, 1, function(x) x/colSums(E65.2.bin20.celltype.df.immune.df)), check.names = F)
E65.2.bin20.celltype.df.immune.df$sample <- "E65.2.bin20"
E65.2.bin20.celltype.df.immune.df$celltype <- rownames(E65.2.bin20.celltype.df.immune.df)

#E75.1
E75.1.bin20.celltype.df.immune <- E75.1.bin20.celltype.df[E75.1.bin20.celltype.df$cell_type %in% immune.cell, ]
E75.1.bin20.celltype.df.immune$cell_type <- factor(E75.1.bin20.celltype.df.immune$cell_type, levels = immune.cell)
E75.1.bin20.celltype.df.immune$region <- "up"
E75.1.bin20.celltype.df.immune[E75.1.bin20.celltype.df.immune$y> -10600, ]$region <- "down"
E75.1.bin20.celltype.df.immune.df <- table(E75.1.bin20.celltype.df.immune$region, E75.1.bin20.celltype.df.immune$cell_type)
E75.1.bin20.celltype.df.immune.df <- data.frame(apply(E75.1.bin20.celltype.df.immune.df, 1, function(x) x/colSums(E75.1.bin20.celltype.df.immune.df)), check.names = F)
E75.1.bin20.celltype.df.immune.df$sample <- "E75.1.bin20"
E75.1.bin20.celltype.df.immune.df$celltype <- rownames(E75.1.bin20.celltype.df.immune.df)

#E75.2
E75.2.bin20.celltype.df.immune <- E75.2.bin20.celltype.df[E75.2.bin20.celltype.df$cell_type %in% immune.cell, ]
E75.2.bin20.celltype.df.immune$cell_type <- factor(E75.2.bin20.celltype.df.immune$cell_type, levels = immune.cell)
median <- (min(E75.2.bin20.celltype.df$y) + max(E75.2.bin20.celltype.df$y))/2
E75.2.bin20.celltype.df.immune$region <- "up"
E75.2.bin20.celltype.df.immune[E75.2.bin20.celltype.df.immune$y< -13300, ]$region <- "down"
E75.2.bin20.celltype.df.immune.df <- table(E75.2.bin20.celltype.df.immune$region, E75.2.bin20.celltype.df.immune$cell_type)
E75.2.bin20.celltype.df.immune.df <- data.frame(apply(E75.2.bin20.celltype.df.immune.df, 1, function(x) x/colSums(E75.2.bin20.celltype.df.immune.df)), check.names = F)
E75.2.bin20.celltype.df.immune.df$sample <- "E75.2.bin20"
E75.2.bin20.celltype.df.immune.df$celltype <- rownames(E75.2.bin20.celltype.df.immune.df)

#E85.C1
E85.C1.bin20.celltype.df.immune <- E85.C1.bin20.celltype.df[E85.C1.bin20.celltype.df$cell_type %in% immune.cell, ]
E85.C1.bin20.celltype.df.immune$cell_type <- factor(E85.C1.bin20.celltype.df.immune$cell_type, levels = immune.cell)
median <- (min(E85.C1.bin20.celltype.df$y) + max(E85.C1.bin20.celltype.df$y))/2
E85.C1.bin20.celltype.df.immune$region <- "up"
E85.C1.bin20.celltype.df.immune[E85.C1.bin20.celltype.df.immune$y>8100, ]$region <- "down"
E85.C1.bin20.celltype.df.immune.df <- table(E85.C1.bin20.celltype.df.immune$region, E85.C1.bin20.celltype.df.immune$cell_type)
E85.C1.bin20.celltype.df.immune.df <- data.frame(apply(E85.C1.bin20.celltype.df.immune.df, 1, function(x) x/colSums(E85.C1.bin20.celltype.df.immune.df)), check.names = F)
E85.C1.bin20.celltype.df.immune.df$sample <- "E85.C1.bin20"
E85.C1.bin20.celltype.df.immune.df <- E85.C1.bin20.celltype.df.immune.df[-5,]
E85.C1.bin20.celltype.df.immune.df$celltype <- rownames(E85.C1.bin20.celltype.df.immune.df)

#E85.B6
E85.B6.bin20.celltype.df.immune <- E85.B6.bin20.celltype.df[E85.B6.bin20.celltype.df$cell_type %in% immune.cell, ]
E85.B6.bin20.celltype.df.immune$cell_type <- factor(E85.B6.bin20.celltype.df.immune$cell_type, levels = immune.cell)
median <- (min(E85.B6.bin20.celltype.df$y) + max(E85.B6.bin20.celltype.df$y))/2
E85.B6.bin20.celltype.df.immune$region <- "up"
E85.B6.bin20.celltype.df.immune[E85.B6.bin20.celltype.df.immune$y< -5700, ]$region <- "down"
E85.B6.bin20.celltype.df.immune.df <- table(E85.B6.bin20.celltype.df.immune$region, E85.B6.bin20.celltype.df.immune$cell_type)
E85.B6.bin20.celltype.df.immune.df <- data.frame(apply(E85.B6.bin20.celltype.df.immune.df, 1, function(x) x/colSums(E85.B6.bin20.celltype.df.immune.df)), check.names = F)
E85.B6.bin20.celltype.df.immune.df$sample <- "E85.B6.bin20"
E85.B6.bin20.celltype.df.immune.df$celltype <- rownames(E85.B6.bin20.celltype.df.immune.df)

#E95
E95.bin20.celltype.df.immune <- E95.bin20.celltype.df[E95.bin20.celltype.df$cell_type %in% immune.cell, ]
E95.bin20.celltype.df.immune$cell_type <- factor(E95.bin20.celltype.df.immune$cell_type, levels = immune.cell)
E95.bin20.celltype.df.immune$region <- "up"
E95.bin20.celltype.df.immune[E95.bin20.celltype.df.immune$x>16600, ]$region <- "down"
E95.bin20.celltype.df.immune.df <- table(E95.bin20.celltype.df.immune$region, E95.bin20.celltype.df.immune$cell_type)
E95.bin20.celltype.df.immune.df <- data.frame(apply(E95.bin20.celltype.df.immune.df, 1, function(x) x/colSums(E95.bin20.celltype.df.immune.df)), check.names = F)
E95.bin20.celltype.df.immune.df$sample <- "E95.bin20"
E95.bin20.celltype.df.immune.df$celltype <- rownames(E95.bin20.celltype.df.immune.df)

all.slide.immune.region <- Reduce(rbind, list(E65.1.bin20.celltype.df.immune.df, E65.2.bin20.celltype.df.immune.df, 
                                              E75.1.bin20.celltype.df.immune.df, E75.2.bin20.celltype.df.immune.df, 
                                              E85.C1.bin20.celltype.df.immune.df, E85.B6.bin20.celltype.df.immune.df, 
                                              E95.bin20.celltype.df.immune.df))
all.slide.immune.region <- melt(all.slide.immune.region)
all.slide.immune.region$celltype <- factor(all.slide.immune.region$celltype, levels = immune.levels[1:7])
pdf("20210827-figures/figure6/Figure6-01-immune-slide-region-ratio-time.pdf", width = 7, height = 6)
ggplot(all.slide.immune.region, mapping = aes(x = variable, y = value, group = sample, color = sample)) +
  geom_point(size = 2, shape = 1) + 
  geom_line(size = 0.5) + 
  scale_color_manual(values = c(rep("#D6604D",2), rep("#F4A582",2), rep("#92C5DE",2), "#4393C3")) +
  facet_wrap(~ celltype) + 
  theme_classic()
dev.off()

all.slide.immune.region <- all.slide.immune.region[all.slide.immune.region$sample %in% c("E65.2.bin20", "E75.1.bin20", "E85.B6.bin20", "E95.bin20"), ]
pdf("20210827-figures/figure6/Figure6-01-immune-slide-region-ratio.pdf", width = 7, height = 4)
ggplot(all.slide.immune.region, mapping = aes(x = variable, y = value, group = celltype, color = celltype)) +
  geom_point(size = 2, shape = 1) + 
  geom_line(size = 0.5) + 
  scale_color_manual(values = immune.color[1:7]) +
  facet_wrap(~ sample) + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", color = "black"))
dev.off()

#挑选通路，看CXCL和CCL的ligand在空间上的表达
UE.harmony.cellchat.df.chemokine <- UE.harmony.cellchat.df[UE.harmony.cellchat.df$pathway_name %in% c("CXCL", "CCL", "MK"), ]
chemokine.ligands <- unique(UE.harmony.cellchat.df.chemokine$ligand)
chemokine.receptors <- c(unique(UE.harmony.cellchat.df.chemokine$receptor)[c(1:13,16,17)], "Itga4", "Itgb1", "Itga6")
pdf("20210827-figures/figure6/Figure6-02-E65-chemokines-expr.pdf", width = 10, height = 10)
my_spatial_plot(E65.2.bin50, features = chemokine.ligands, ncol = 4, color = brewer.pal(9, "PuBu"), ncol = 4)
dev.off()
pdf("20210827-figures/figure6/Figure6-02-E75-chemokines-expr.pdf", width = 12, height = 12)
my_spatial_plot(E75.2.bin50, features = chemokine.ligands, color = brewer.pal(9, "PuBu"), ncol = 4)
dev.off()
pdf("20210827-figures/figure6/Figure6-02-E75-1-chemokines-expr.pdf", width = 10, height = 10)
my_spatial_plot(E75.1.bin50, features = chemokine.ligands, color = brewer.pal(9, "PuBu"), ncol = 4)
dev.off()
load("20210818-spatial/B6-bin50-spatial-phago-map.RData")
pdf("20210827-figures/figure6/Figure6-02-E85-chemokines-expr.pdf", width = 12, height = 12)
my_spatial_plot(B6.bin50, features = chemokine.ligands, color = brewer.pal(9, "PuBu"), ncol = 4)
dev.off()
load("20210818-spatial/E95-2-bin50-spatial-phago-map.RData")
pdf("20210827-figures/figure6/Figure6-02-E95-chemokines-expr.pdf", width = 14, height = 12)
my_spatial_plot(E95.bin50, features = chemokine.ligands, color = brewer.pal(9, "PuBu"), ncol = 4)
dev.off()
pdf("20210827-figures/figure6/Figure6-02-UMAP-chemokines-expr.pdf", width = 25, height = 20)
FeaturePlot(UE.decidual, features = chemokine.ligands, pt.size = 0.5, ncol = 5) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()
#所有免疫细胞，DSCs 和 iDSC的0、1、2三群看chemokines的表达情况
load("../20210201-integrate/UE-immune-hamrony-rename.RData")
load("../20210201-integrate/UE-decidual-harmony-pc30-rename.RData")
UE.decidual$cellname <- gsub("Postn-D", "Lum-D", UE.decidual$cellname)
UE.decidual$cellname <- gsub("Ifit1-D", "S100a8-D", UE.decidual$cellname)
UE.decidual$cellname <- gsub("Hist1h2ap-D", "Top2a-D", UE.decidual$cellname)
UE.decidual <- subset(UE.decidual, subset = cellname != "Prg-" & cellname != "cluster9-1")
load("20210827-figures/figure5/FB-immune-sc-subclusters-SCT.RData")
FB.immune.SCT$cellname <- paste("iDSC", FB.immune.SCT@active.ident, sep = "")
UE.immune.harmony <- subset(UE.immune.harmony, subset = cellname != "FB-immune")
UE.immune.DSC <- merge(UE.immune.harmony, list(UE.decidual, FB.immune.SCT))
UE.immune.DSC@active.ident <- factor(UE.immune.DSC$cellname)
UE.immune.DSC <- subset(UE.immune.DSC, subset = cellname != "B" & cellname != "T" & cellname != "Ptn-D")
UE.immune.DSC <- NormalizeData(UE.immune.DSC, normalization.method = "RC", scale.factor = 1000000)#CPM normalize
UE.immune.DSC <- ScaleData(UE.immune.DSC, features = rownames(UE.immune.DSC))
UE.immune.DSC.avg <- AverageExpression(UE.immune.DSC, slot = "data", assays = "RNA", use.counts = T)
UE.immune.DSC.avg <- UE.immune.DSC.avg$RNA
UE.immune.DSC.avg.chemokine <- UE.immune.DSC.avg[c(chemokine.ligands, chemokine.receptors), c(7,15,16,3,14,12,4:6,8,13,1,2,10,11)]
UE.immune.DSC.avg.chemokine.L <- UE.immune.DSC.avg.chemokine[c(1:19), c(1:9)]
UE.immune.DSC.avg.chemokine.R <- UE.immune.DSC.avg.chemokine[c(20:37), c(10:15)]
color <- rev(brewer.pal(11, "Spectral"))
pdf("20210827-figures/figure6/Figure6-02-scRNA-DSC-chemokines-L-heatmap.pdf", width = 10, height = 10)
pheatmap(log2(UE.immune.DSC.avg.chemokine.L+1), scale = "none", cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(colors = color)(100))
dev.off()
pdf("20210827-figures/figure6/Figure6-02-scRNA-immune-chemokines-R-heatmap.pdf", width = 10, height = 10)
pheatmap(log2(UE.immune.DSC.avg.chemokine.R+1), scale = "none", cluster_rows = F, cluster_cols = F, 
         color = colorRampPalette(colors = color)(100))
dev.off()
#尝试画CXCL相关的network plot
library(ggraph)
library(igraph)
library(tidyverse)
UE.harmony.cellchat.df.CXCL <- UE.harmony.cellchat.df[UE.harmony.cellchat.df$pathway_name == "CXCL", ]
UE.harmony.cellchat.df.CXCL <- UE.harmony.cellchat.df.CXCL[UE.harmony.cellchat.df.CXCL$source %in% c(decidual.cell, immune.cell[c(1:5,7)], "FB-immune1") &
                                                             UE.harmony.cellchat.df.CXCL$target %in% c(immune.cell[c(1:5,7)], "FB-immune1"), ]
UE.harmony.cellchat.df.CXCL <- UE.harmony.cellchat.df.CXCL[, c(1,2,5,8,9)]
node <- data.frame(celltype = c(decidual.cell, immune.cell, "FB-immune0", "FB-immune1", "FB-immune2"), 
                   majortype = c(rep("Decidual", length(decidual.cell)), rep("Immune", length(immune.cell[c(1:5,7)])), "FB-immune"))
#使用igraph包创建网络数据
CXCL.g <- graph_from_data_frame(UE.harmony.cellchat.df.CXCL, directed = TRUE, vertices = node)
CXCL.g.gt <- as_tbl_graph(CXCL.g)
edge.color <- colorRampPalette(colors = brewer.pal(9, "Greys"))(20)
edge.color <- edge.color[seq(7,20,2)]
node.color <- c("#33A02C", "#FB9A99", "#6A3D9A", "#08519C", "#6A5ACD", "#1F78B4", "#A6CEE3", "#2E8B57", 
                "#CAB2D6", "#228B22", "#1F78B4", "#FDB462", "#9ACD32", "#FB8072")
pdf("20210827-figures/figure6/Figure6-02-CXCL-circle-plot-legend.pdf", width = 10, height = 10)
ggraph(CXCL.g.gt, layout = "circle") + 
  geom_edge_fan(aes(edge_width = prob, color = interaction_name_2), #fan可以画两个点之间多个线
                strength = 1, edge_alpha = 0.8, 
                 arrow = arrow(length = unit(3, 'mm'), type = "closed"), 
                 end_cap = circle(5, 'mm')) + 
  scale_edge_width(range=c(1, 2)) +
  scale_edge_color_manual(values = edge.color) + 
  geom_node_point(aes(color = name), size = 10) + 
  #geom_node_text(aes(label = name)) + 
  scale_color_manual(values = node.color) + 
  theme_void()# +
  #theme(legend.position = "none")
dev.off()

#CCL
UE.harmony.cellchat.df.CCL <- UE.harmony.cellchat.df[UE.harmony.cellchat.df$pathway_name == "CCL", ]
UE.harmony.cellchat.df.CCL <- UE.harmony.cellchat.df.CCL[UE.harmony.cellchat.df.CCL$source %in% c(decidual.cell, immune.cell[c(1:5,7)], "FB-immune1") &
                                                             UE.harmony.cellchat.df.CCL$target %in% c(immune.cell[c(1:5,7)], "FB-immune1"), ]
UE.harmony.cellchat.df.CCL <- UE.harmony.cellchat.df.CCL[, c(1,2,5,8,9)]
CCL.g <- graph_from_data_frame(UE.harmony.cellchat.df.CCL, directed = TRUE, vertices = node)
CCL.g.gt <- as_tbl_graph(CCL.g)
edge.color <- colorRampPalette(colors = brewer.pal(9, "Greys"))(20)
edge.color <- edge.color[seq(7,20)]

pdf("20210827-figures/figure6/Figure6-02-CCL-circle-plot.pdf", width = 5, height = 5)
ggraph(CCL.g.gt, layout = "circle") + 
  geom_edge_fan(aes(edge_width = prob, color = interaction_name_2), #fan可以画两个点之间多个线
                strength = 1, edge_alpha = 0.8, 
                arrow = arrow(length = unit(3, 'mm'), type = "closed"), 
                end_cap = circle(5, 'mm')) + 
  scale_edge_width(range=c(1, 2)) +
  scale_edge_color_manual(values = edge.color) + 
  geom_node_point(aes(color = name), size = 10) + 
  #geom_node_text(aes(label = name)) + 
  scale_color_manual(values = node.color) + 
  theme_void() + 
  theme(legend.position = "none")
dev.off()

#挑选CCL6，CCL9，CXCL12三个L-R
UE.harmony.cellchat.df.select <- UE.harmony.cellchat.df[UE.harmony.cellchat.df$ligand %in% c("Cxcl12", "Ccl6", "Ccl9") & 
                                                          UE.harmony.cellchat.df$source %in% c(c(decidual.cell, "FB-immune0", "FB-immune1", "FB-immune2")) &
                                                          UE.harmony.cellchat.df$target %in% c(immune.cell[c(1:5,7)]), ]
UE.harmony.cellchat.df.select <- UE.harmony.cellchat.df.select[, c(1,2,5,8,9)]
node <- data.frame(celltype = c(decidual.cell, immune.cell[c(1:5,7)], "FB-immune0", "FB-immune1", "FB-immune2"), 
                   majortype = c(rep("Decidual", length(decidual.cell)), rep("Immune", length(immune.cell[c(1:5,7)])), rep("FB-immune", 3)))
select.g <- graph_from_data_frame(UE.harmony.cellchat.df.select, directed = TRUE, vertices = node)
select.g.gt <- as_tbl_graph(select.g)
edge.color <- colorRampPalette(colors = brewer.pal(9, "Greys"))(20)
edge.color <- edge.color[seq(7,20,4)]
pdf("20210827-figures/figure6/Figure6-03-CCL6-CCL9-CXCL12-circleplot-legend.pdf", width = 10, height = 10)
ggraph(select.g.gt, layout = "circle") + 
  geom_edge_fan(aes(edge_width = prob, color = interaction_name_2), #fan可以画两个点之间多个线
                strength = 1, edge_alpha = 0.8, 
                arrow = arrow(length = unit(3, 'mm'), type = "closed"), 
                end_cap = circle(5, 'mm')) + 
  scale_edge_width(range=c(1, 2)) +
  scale_edge_color_manual(values = edge.color) + 
  geom_node_point(aes(color = name), size = 10) + 
  #geom_node_text(aes(label = name)) + 
  #scale_color_manual(values = node.color) + 
  theme_void() #+ 
  #theme(legend.position = "none")
dev.off()
#CCL6, CCL9, CXCL12和GALACTIN通路相关的基因表达 L-R
UE.harmony.cellchat.subset <- subsetCellChat(UE.harmony.cellchat, idents.use = c(decidual.cell, immune.cell, "FB-immune0", "FB-immune1", "FB-immune2"))
color <- c(decidual.color[c(1:3,5,7,6)], immune.color[c(4,5,3,2,7,6)], "#F0746F", "#39B54E", "#6292CC")
pdf("20210827-figures/figure6/Figure6-05-CXCL-violin.pdf", width = 10, height = 10)
plotGeneExpression(UE.harmony.cellchat.subset, signaling = "CXCL", color.use = color)
dev.off()
pdf("20210827-figures/figure6/Figure6-05-CCL-violin.pdf", width = 10, height = 10)
plotGeneExpression(UE.harmony.cellchat.subset, signaling = "CCL", color.use = color)
dev.off()
pdf("20210827-figures/figure6/Figure6-05-GALECTIN-violin.pdf", width = 10, height = 10)
plotGeneExpression(UE.harmony.cellchat.subset, signaling = "GALECTIN", color.use = color)
dev.off()
pdf("20210827-figures/figure6/Figure6-05-GALECTIN-violin.pdf", width = 10, height = 10)
plotGeneExpression(UE.harmony.cellchat.subset, signaling = "GALECTIN", color.use = color)
dev.off()
pdf("20210827-figures/figure6/Figure6-05-MK-violin.pdf", width = 10, height = 10)
plotGeneExpression(UE.harmony.cellchat.subset, signaling = "MK", color.use = color)
dev.off()
#增殖相关
UE.harmony.cellchat.df.pro <- UE.harmony.cellchat.df[UE.harmony.cellchat.df$source %in% c(immune.cell, "FB-immune0", "FB-immune1", "FB-immune2") &
                                                             UE.harmony.cellchat.df$target %in% c("Top2a-D"), ]
pdf("20210827-figures/figure6/Figure6-04-immune-cycling-DSC-dotplot.pdf", width = 4, height = 8)
netVisual_bubble(UE.harmony.cellchat, sources.use = c(13,12,17,10,11,16,32:34), targets.use = 4)
dev.off()

UE.harmony.cellchat.immune <- subsetCellChat(UE.harmony.cellchat, idents.use = c(decidual.levels[c(1:5,7,6)], immune.levels[1:7], "FB-immune0", "FB-immune1", "FB-immune2"))
color1 <- c(decidual.color[c(1:5,7,6)], immune.color[c(4,5,3,2,7,6)], rep(immune.color[10], 3))
pdf("20210827-figures/figure6/Figure6-07-immune-NK-releasing.pdf", width = 15, height = 10)
plotGeneExpressioncc(UE.harmony.cellchat.immune, signaling = c("SPP1", "PARs", "XCR", "FLT3"), color.use = color1)#NK releasing
dev.off()
pdf("20210827-figures/figure6/Figure6-07-immune-Neutrophil-releasing.pdf", width = 15, height = 10)
plotGeneExpression(UE.harmony.cellchat.immune, signaling = c("OSM", "CSF", "IL1"), color.use = color1)#Neutrophil releasing
dev.off()

pdf("20210827-figures/figure6/Figure6-immune-CXCL-CCL-dotplot.pdf", width = 6, height = 4.5)
netVisual_bubble(UE.harmony.cellchat, sources.use = c(decidual.levels[c(1:2)], "FB-immune0","FB-immune1"), 
                 targets.use = c(immune.levels[1:7]), 
                 signaling = c("CCL","CXCL","MK"), remove.isolate = FALSE)
dev.off()

#计算配体-受体对的基因的表达值
load("20210818-spatial/UE-harmony-pc30-UMAP-subname-iDSC.RData")
UE.harmony.SCT.iDSC@active.ident <- factor(UE.harmony.SCT.iDSC$subname)
UE.harmony.subname.avg <- AverageExpression(UE.harmony.SCT.iDSC, assay = "RNA")
UE.harmony.subname.avg <- UE.harmony.subname.avg$RNA
selected.path <- c("CXCL", "CCL")
UE.harmony.cellchat.CXCL.CCL <- UE.harmony.cellchat.df[UE.harmony.cellchat.df$source %in% c(decidual.levels[c(1:3, 5:7)], "FB-immune0","FB-immune1","FB-immune2") & 
                                                       UE.harmony.cellchat.df$target %in% c(immune.levels[1:7]) & 
                                                       UE.harmony.cellchat.df$pathway_name %in% selected.path, ]
value.c <- c()
for (i in seq(1, dim(UE.harmony.cellchat.CXCL.CCL)[1], 1)) {
  print(i)
  source <- as.character(UE.harmony.cellchat.CXCL.CCL[i, 1])
  target <- as.character(UE.harmony.cellchat.CXCL.CCL[i, 2])
  ligand <- as.character(UE.harmony.cellchat.CXCL.CCL[i, 3])
  receptor <- as.character(UE.harmony.cellchat.CXCL.CCL[i, 4])
  value <- UE.harmony.subname.avg[ligand, source] * UE.harmony.subname.avg[receptor, target]
  value.c <- c(value.c, value)
}
UE.harmony.cellchat.CXCL.CCL$express <- value.c
UE.harmony.cellchat.CXCL.CCL <- UE.harmony.cellchat.CXCL.CCL[UE.harmony.cellchat.CXCL.CCL$express>2, ]
UE.harmony.cellchat.CXCL.CCL <- UE.harmony.cellchat.CXCL.CCL[, c(1,2,5,8,9,12)]
UE.harmony.cellchat.CXCL.CCL$source <- factor(UE.harmony.cellchat.CXCL.CCL$source, 
                                              levels = c(decidual.levels[c(1:3,5:7)], immune.levels[c(1:7)], "FB-immune0", "FB-immune1", "FB-immune2"))
node <- data.frame(celltype = c(decidual.levels[c(1:3,5:7)], immune.levels[c(1:7)], "FB-immune0", "FB-immune1", "FB-immune2"), 
                   majortype = c(rep("Decidual", 6), rep("Immune", 7), rep("FB-immune", 3)))
node.color <- c(decidual.color[c(1:3,5:7)], immune.color[c(1:7)], rep("purple", 3))
select.g <- graph_from_data_frame(UE.harmony.cellchat.CXCL.CCL, directed = TR)
UE.harmony.cellchat.df.pattern$source_target <- paste(UE.harmony.cellchat.df.pattern$source, UE.harmony.cellchat.df.pattern$target, sep = "_")
UE.harmony.cellchat.df.pattern1 <- UE.harmony.cellchat.df.pattern[UE.harmony.cellchat.df.pattern$interaction_name_2 %in% pattern1, ]
UE.harmony.cellchat.df.pattern1 <- aggregate(UE.harmony.cellchat.df.pattern1$prob, by = list(UE.harmony.cellchat.df.pattern1$source_target), sum)
pattern1.add <- data.frame(Group.1 = c("FB-immune0_Mac"), x = c(0.000001))
UE.harmony.cellchat.df.pattern1 <- rbind(UE.harmony.cellchat.df.pattern1, pattern1.add)
UE.harmony.cellchat.df.pattern1 <- separate(UE.harmony.cellchat.df.pattern1, "Group.1", into = c("source", "target"), sep = "_")
chordDiagram(UE.harmony.cellchat.df.pattern1)
circos.clear()
grid.col <- c("Lum-D" = "#6A5ACD", "Sfrp4-D" = "#9ACD32", "FB-immune0cF0746F", "FB-immune1" = "#39B54E", 
              "Mac" = "#1F78B4", "Pro-Mac" = "#B2DF8A", "DC-1" = "#33A02C", "DC-2" = "#FB9A99", "Neutrophil" = "#2E8B57", "NK" = "#CAB2D6")
pdf("20210827-figures/figure6/Figure6-pattern1-circos.pdf", width = 6, height = 4.5)
chordDiagram(UE.harmony.cellchat.df.pattern1, 
             order = c(decidual.levels[1:2], "FB-immune0", "FB-immune1", immune.levels[2:7]), 
             grid.col = grid.col, directional = 1, 
             direction.type = c("diffHeight", "arrows"), 
             link.arr.type = "big.arrow", 
             annotationTrack = c("name", "grid"), 
             scale = "TRUE", link.visible = UE.harmony.cellchat.df.pattern1[[3]] >= 0.05) 
dev.off()

UE.harmony.cellchat.df.pattern2 <- UE.harmony.cellchat.df.pattern[UE.harmony.cellchat.df.pattern$interaction_name_2 %in% pattern2, ]
UE.harmony.cellchat.df.pattern2 <- aggregate(UE.harmony.cellchat.df.pattern2$prob, by = list(UE.harmony.cellchat.df.pattern2$source_target), sum)
pattern2.add <- data.frame(Group.1 = c("Sfrp4-D_Mac", "Lum-D_Mac"), x = c(0.000001, 0.000001))
UE.harmony.cellchat.df.pattern2 <- rbind(UE.harmony.cellchat.df.pattern2, pattern2.add)
UE.harmony.cellchat.df.pattern2 <- separate(UE.harmony.cellchat.df.pattern2, "Group.1", into = c("source", "target"), sep = "_")
pdf("20210827-figures/figure6/Figure6-pattern2-circos.pdf", width = 6, height = 4.5)
chordDiagram(UE.harmony.cellchat.df.pattern2, 
             order = c(decidual.levels[1:2], "FB-immune0", "FB-immune1", immune.levels[2:7]), 
             grid.col = grid.col, directional = 1, 
             direction.type = c("diffHeight", "arrows"), 
             link.arr.type = "big.arrow", 
             annotationTrack = c("name", "grid"), 
             scale = TRUE, link.visible = UE.harmony.cellchat.df.pattern2[[3]] >= 0.00001) 
dev.off()

UE.harmony.cellchat.df.pattern3 <- UE.harmony.cellchat.df.pattern[UE.harmony.cellchat.df.pattern$interaction_name_2 %in% pattern3, ]
UE.harmony.cellchat.df.pattern3 <- aggregate(UE.harmony.cellchat.df.pattern3$prob, by = list(UE.harmony.cellchat.df.pattern3$source_target), sum)
pattern3.add <- data.frame(Group.1 = c("Sfrp4-D_Pro-Mac", "Sfrp4-D_DC-1", "Sfrp4-D_NK"), x = c(rep(0.000001, 3)))
UE.harmony.cellchat.df.pattern3 <- rbind(UE.harmony.cellchat.df.pattern3, pattern3.add)
UE.harmony.cellchat.df.pattern3 <- separate(UE.harmony.cellchat.df.pattern3, "Group.1", into = c("source", "target"), sep = "_")
circos.par(gap.after = c(rep(1, length(unique(UE.harmony.cellchat.df.pattern3$source))-1), 10,
                         rep(1, length(unique(UE.harmony.cellchat.df.pattern3$target))-1), 10)) 
pdf("20210827-figures/figure6/Figure6-pattern3-circos.pdf", width = 6, height = 4.5)
chordDiagram(UE.harmony.cellchat.df.pattern3, 
             order = c(decidual.levels[1:2], "FB-immune0", "FB-immune1", immune.levels[2:7]), 
             grid.col = grid.col, directional = 1, 
             direction.type = c("diffHeight", "arrows"), 
             link.arr.type = "big.arrow", 
             annotationTrack = c("name", "grid"),
             scale = TRUE, link.visible = UE.harmony.cellchat.df.pattern3[[3]] >= 0.00001) 
dev.off()

UE.harmony.cellchat.df.pattern4 <- UE.harmony.cellchat.df.pattern[UE.harmony.cellchat.df.pattern$interaction_name_2 %in% pattern4, ]
UE.harmony.cellchat.df.pattern4 <- aggregate(UE.harmony.cellchat.df.pattern4$prob, by = list(UE.harmony.cellchat.df.pattern4$source_target), sum)
pattern4.add <- data.frame(Group.1 = c("Sfrp4-D_Pro-Mac"), x = c(rep(0.000001, 1)))
UE.harmony.cellchat.df.pattern4 <- rbind(UE.harmony.cellchat.df.pattern4, pattern4.add)
UE.harmony.cellchat.df.pattern4 <- separate(UE.harmony.cellchat.df.pattern4, "Group.1", into = c("source", "target"), sep = "_")
circos.par(gap.after = c(rep(1, length(unique(UE.harmony.cellchat.df.pattern4$source))-1), 10,
                         rep(1, length(unique(UE.harmony.cellchat.df.pattern4$target))-1), 10)) 
pdf("20210827-figures/figure6/Figure6-pattern4-circos.pdf", width = 6, height = 4.5)
chordDiagram(UE.harmony.cellchat.df.pattern4, 
             order = c(decidual.levels[1:2], "FB-immune0", "FB-immune1", immune.levels[2:7]), 
             grid.col = grid.col, directional = 1, 
             direction.type = c("diffHeight", "arrows"), 
             link.arr.type = "big.arrow", 
             annotationTrack = c("name", "grid"),
             scale = TRUE, link.visible = UE.harmony.cellchat.df.pattern4[[3]] >= 0.00001) 
dev.off()


#immune closed的细胞（bin50）
all.bin50.closed <- read.table("20210827-figures/files/all-slides-bin50-closed-cells-calc.csv", sep = ",", header = T, check.names = F)
all.bin50.closed.immune <- all.bin50.closed[all.bin50.closed$source %in% immune.levels[1:7], ]
all.bin50.closed.immune <- all.bin50.closed.immune[all.bin50.closed.immune$bin50.closed %in% 
                                                     c(decidual.levels[c(1:3,5:7)]), ]
all.bin50.closed.immune$bin50.closed <- factor(all.bin50.closed.immune$bin50.closed, 
                                               levels = c(decidual.levels[c(1:3,5:7)]))
all.bin50.closed.immune$source <- factor(all.bin50.closed.immune$source, levels = immune.levels[1:7])
all.bin50.closed.immune[all.bin50.closed.immune==0] <- NA
all.bin50.closed.immune <- all.bin50.closed.immune[all.bin50.closed.immune$sample %in% c("E65.2.bin50", "E75.1.bin50", "E85.B6.bin50", "E95.bin50"), ]
pdf("20210827-figures/figure6/Figure6-01-immune-DSC-closed-bin50.pdf", width = 7, height = 5)
ggplot() + 
  geom_point(all.bin50.closed.immune, mapping = aes(x = source, y = bin50.closed, size = value, color = source)) + 
  facet_wrap(~ sample) +
  scale_size(range = c(0, 10)) + 
  scale_color_manual(values = c(immune.color[1:7])) +
  theme(axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", color = "black"))
dev.off()

#20230201
load("../20210201-integrate/UE-immune-hamrony-rename.RData")
DimPlot(UE.immune.harmony, label = T)
UE.immune.per.df <- data.frame(table(UE.immune.harmony$time, UE.immune.harmony@active.ident), check.names = F)
UE.immune.per.df$Var2 <- factor(UE.immune.per.df$Var2, levels = immune.cell)
pdf("20220622-Comments/10.Adding/01.UE-immune-proportion-barplot.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(UE.immune.per.df, mapping = aes(x = Var1, y = Freq, fill = Var2), stat = "identity", position = "fill") + 
  scale_fill_manual(values = immune.color) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "black", fill = "white"))
dev.off()

UE.immune.markers <- read.table("20210827-figures/files/UE-immune-harmony-markers.csv", sep = ",", check.names = F)
UE.immune.markers.top20 <- UE.immune.markers %>% group_by(cluster) %>% top_n(n = 20, avg_logFC)
UE.immune.harmony <- NormalizeData(UE.immune.harmony)
UE.immune.harmony <- ScaleData(UE.immune.harmony, features = rownames(UE.immune.harmony))
pdf("20220622-Comments/10.Adding/02.UE-immune-DEGs-barplot.pdf", width = 10, height = 10)
DoHeatmap(UE.immune.harmony, features = as.character(UE.immune.markers.top20$gene)) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))
dev.off()


