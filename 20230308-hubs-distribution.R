#Display R default color
library(scales)
show_col(hue_pal()(3))
#Figure5: immune cells + eSF + iDSC0 + iDSC1 + iDSC2
#read files
E65.2.bin50.celltype.df <- read.table("20210827-figures/figure5/E65-2-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
E75.1.bin50.celltype.df <- read.table("20210827-figures/figure5/E75-1-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
E85.B6.bin50.celltype.df <- read.table("20210827-figures/figure5/E85-B6-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
E95.bin50.celltype.df <- read.table("20210827-figures/figure5/E95-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)

##################################################### Normal - bin50 immune cell distribution and zoom in#################################################
#A. Immune subclusters + iDSC subclusters + eSF and zoom in pics
function(){
  #E65.2.bin50
  E65.2.bin50.celltype.immune.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% 
                                                              c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E65.2.bin50.celltype.immune.df$cell_type <- factor(E65.2.bin50.celltype.immune.df$cell_type, 
                                                     levels= c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E65.2.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin50-immune-rect.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin50-immune.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin50.celltype.immune.df.annotation <- E65.2.bin50.celltype.immune.df[E65.2.bin50.celltype.immune.df$x >= 7000 & E65.2.bin50.celltype.immune.df$x <= 8000 & 
                                                                                E65.2.bin50.celltype.immune.df$y >= 16700 & E65.2.bin50.celltype.immune.df$y <= 17700, ]
  E65.2.bin50.celltype.df.annotation <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$x >= 7000 & E65.2.bin50.celltype.df$x <= 8000 & 
                                                                  E65.2.bin50.celltype.df$y >= 16700 & E65.2.bin50.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin50.celltype.immune.df$cell_type), unique(E65.2.bin50.celltype.immune.df.annotation$cell_type)), 
                                       levels(E65.2.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin50-immune-annotation.pdf", width = 1.2, height = 1.1)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E75.1.bin50
  E75.1.bin50.celltype.immune.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% 
                                                              c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E75.1.bin50.celltype.immune.df$cell_type <- factor(E75.1.bin50.celltype.immune.df$cell_type, 
                                                     levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E75.1.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin50-immune-rect.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin50-immune.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E75.1.bin50.celltype.immune.df.annotation <- E75.1.bin50.celltype.immune.df[E75.1.bin50.celltype.immune.df$x >= 7000 & E75.1.bin50.celltype.immune.df$x <= 9000 & 
                                                                                E75.1.bin50.celltype.immune.df$y <= 15000 & E75.1.bin50.celltype.immune.df$y >= 13000, ]
  E75.1.bin50.celltype.df.annotation <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$x >= 7000 & E75.1.bin50.celltype.df$x <= 9000 & 
                                                                  E75.1.bin50.celltype.df$y <= 15000 & E75.1.bin50.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin50.celltype.immune.df$cell_type), unique(E75.1.bin50.celltype.immune.df.annotation$cell_type)), 
                                       levels(E75.1.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin50-immune-annotation.pdf", width = 2, height = 2)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin50
  E85.B6.bin50.celltype.immune.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% 
                                                                c("Mono","Mac","Pro-Mac","DC-1","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E85.B6.bin50.celltype.immune.df$cell_type <- factor(E85.B6.bin50.celltype.immune.df$cell_type, 
                                                      levels = c("Mono","Mac","Pro-Mac","DC-1","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E85.B6.color <- c(immune.color[c(1:4,7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin50-immune-rect.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin50-immune.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin50.celltype.immune.df.annotation <- E85.B6.bin50.celltype.immune.df[E85.B6.bin50.celltype.immune.df$x >= 4000 & E85.B6.bin50.celltype.immune.df$x <= 6500 & 
                                                                                  E85.B6.bin50.celltype.immune.df$y >= -4000 & E85.B6.bin50.celltype.immune.df$y <= -2000, ]
  E85.B6.bin50.celltype.df.annotation <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$x >= 4000 & E85.B6.bin50.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin50.celltype.df$y >= -4000 & E85.B6.bin50.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin50.celltype.immune.df$cell_type), unique(E85.B6.bin50.celltype.immune.df.annotation$cell_type)), 
                                         levels(E85.B6.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin50-immune-annotation.pdf", width = 2.3, height = 1.9)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin50
  E95.bin50.celltype.immune.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% 
                                                          c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E95.bin50.celltype.immune.df$cell_type <- factor(E95.bin50.celltype.immune.df$cell_type, 
                                                   levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E95.bin50.celltype.immune.eSF.df <- E95.bin50.celltype.immune.df[E95.bin50.celltype.immune.df$cell_type %in% "Lum-D", ]
  E95.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin50-immune-rect.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin50-immune.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    geom_point(E95.bin50.celltype.immune.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 1.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin50.celltype.immune.df.annotation <- E95.bin50.celltype.immune.df[E95.bin50.celltype.immune.df$x >= 11300 & E95.bin50.celltype.immune.df$x <= 13800 & 
                                                                            E95.bin50.celltype.immune.df$y >= 17000 & E95.bin50.celltype.immune.df$y <= 20000, ]
  E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df[E95.bin50.celltype.df$x >= 11300 & E95.bin50.celltype.df$x <= 13800 & 
                                                              E95.bin50.celltype.df$y >= 17000 & E95.bin50.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin50.celltype.immune.df$cell_type), unique(E95.bin50.celltype.immune.df.annotation$cell_type)), 
                                   levels(E95.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin50-immune-annotation.pdf", width = 2.3, height = 2.7)
  ggplot() + 
    geom_point(E95.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
}
#B. Immune subclusters and zoom in pics
function(){
  #E65.2.bin50
  E65.2.bin50.celltype.immune.B.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% 
                                                              c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E65.2.bin50.celltype.immune.B.df$cell_type <- factor(E65.2.bin50.celltype.immune.B.df$cell_type, 
                                                     levels= c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"))
  E65.2.color <- c(immune.color[c(1:7)])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin50-immune-B-rect.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.B.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin50-immune-B.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.B.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin50.celltype.immune.B.df.annotation <- E65.2.bin50.celltype.immune.B.df[E65.2.bin50.celltype.immune.B.df$x >= 7000 & E65.2.bin50.celltype.immune.B.df$x <= 8000 & 
                                                                                E65.2.bin50.celltype.immune.B.df$y >= 16700 & E65.2.bin50.celltype.immune.B.df$y <= 17700, ]
  E65.2.bin50.celltype.df.annotation <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$x >= 7000 & E65.2.bin50.celltype.df$x <= 8000 & 
                                                                  E65.2.bin50.celltype.df$y >= 16700 & E65.2.bin50.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin50.celltype.immune.B.df$cell_type), unique(E65.2.bin50.celltype.immune.B.df.annotation$cell_type)), 
                                       levels(E65.2.bin50.celltype.immune.B.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin50-immune-B-annotation.pdf", width = 1.2, height = 1.1)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.B.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E75.1.bin50
  E75.1.bin50.celltype.immune.B.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% 
                                                              c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E75.1.bin50.celltype.immune.B.df$cell_type <- factor(E75.1.bin50.celltype.immune.B.df$cell_type, 
                                                     levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"))
  E75.1.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin50-immune-B-rect.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.B.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin50-immune-B.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.B.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E75.1.bin50.celltype.immune.B.df.annotation <- E75.1.bin50.celltype.immune.B.df[E75.1.bin50.celltype.immune.B.df$x >= 7000 & E75.1.bin50.celltype.immune.B.df$x <= 9000 & 
                                                                                E75.1.bin50.celltype.immune.B.df$y <= 15000 & E75.1.bin50.celltype.immune.B.df$y >= 13000, ]
  E75.1.bin50.celltype.df.annotation <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$x >= 7000 & E75.1.bin50.celltype.df$x <= 9000 & 
                                                                  E75.1.bin50.celltype.df$y <= 15000 & E75.1.bin50.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin50.celltype.immune.B.df$cell_type), unique(E75.1.bin50.celltype.immune.B.df.annotation$cell_type)), 
                                       levels(E75.1.bin50.celltype.immune.B.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin50-immune-B-annotation.pdf", width = 2, height = 2)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.B.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin50
  E85.B6.bin50.celltype.immune.B.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% 
                                                                c("Mono","Mac","Pro-Mac","DC-1","NK"), ]
  E85.B6.bin50.celltype.immune.B.df$cell_type <- factor(E85.B6.bin50.celltype.immune.B.df$cell_type, 
                                                      levels = c("Mono","Mac","Pro-Mac","DC-1","NK"))
  E85.B6.color <- c(immune.color[c(1:4,7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin50-immune-B-rect.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.B.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin50-immune-B.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.B.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin50.celltype.immune.B.df.annotation <- E85.B6.bin50.celltype.immune.B.df[E85.B6.bin50.celltype.immune.B.df$x >= 4000 & E85.B6.bin50.celltype.immune.B.df$x <= 6500 & 
                                                                                  E85.B6.bin50.celltype.immune.B.df$y >= -4000 & E85.B6.bin50.celltype.immune.B.df$y <= -2000, ]
  E85.B6.bin50.celltype.df.annotation <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$x >= 4000 & E85.B6.bin50.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin50.celltype.df$y >= -4000 & E85.B6.bin50.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin50.celltype.immune.B.df$cell_type), unique(E85.B6.bin50.celltype.immune.B.df.annotation$cell_type)), 
                                         levels(E85.B6.bin50.celltype.immune.B.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin50-immune-B-annotation.pdf", width = 2.3, height = 1.9)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.B.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin50
  E95.bin50.celltype.immune.B.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% 
                                                          c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E95.bin50.celltype.immune.B.df$cell_type <- factor(E95.bin50.celltype.immune.B.df$cell_type, 
                                                   levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"))
  E95.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin50-immune-B-rect.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.B.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin50-immune-B.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.B.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin50.celltype.immune.B.df.annotation <- E95.bin50.celltype.immune.B.df[E95.bin50.celltype.immune.B.df$x >= 11300 & E95.bin50.celltype.immune.B.df$x <= 13800 & 
                                                                            E95.bin50.celltype.immune.B.df$y >= 17000 & E95.bin50.celltype.immune.B.df$y <= 20000, ]
  E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df[E95.bin50.celltype.df$x >= 11300 & E95.bin50.celltype.df$x <= 13800 & 
                                                              E95.bin50.celltype.df$y >= 17000 & E95.bin50.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin50.celltype.immune.B.df$cell_type), unique(E95.bin50.celltype.immune.B.df.annotation$cell_type)), 
                                   levels(E95.bin50.celltype.immune.B.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin50-immune-B-annotation.pdf", width = 2.3, height = 2.7)
  ggplot() + 
    geom_point(E95.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.B.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
}
#C. Immune major cluster + iDSC subclusters + eSF and zoom in pics
function(){
  #E65.2.bin50
  E65.2.bin50.celltype.immune.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% 
                                                              c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E65.2.bin50.celltype.immune.df$cell_type <- factor(E65.2.bin50.celltype.immune.df$cell_type, 
                                                     levels= c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E65.2.color <- c(rep("#D2691E", 7), decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin50-immune-C-rect.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin50-immune-C.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin50.celltype.immune.df.annotation <- E65.2.bin50.celltype.immune.df[E65.2.bin50.celltype.immune.df$x >= 7000 & E65.2.bin50.celltype.immune.df$x <= 8000 & 
                                                                                E65.2.bin50.celltype.immune.df$y >= 16700 & E65.2.bin50.celltype.immune.df$y <= 17700, ]
  E65.2.bin50.celltype.df.annotation <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$x >= 7000 & E65.2.bin50.celltype.df$x <= 8000 & 
                                                                  E65.2.bin50.celltype.df$y >= 16700 & E65.2.bin50.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin50.celltype.immune.df$cell_type), unique(E65.2.bin50.celltype.immune.df.annotation$cell_type)), 
                                       levels(E65.2.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin50-immune-C-annotation.pdf", width = 1.2, height = 1.1)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E75.1.bin50
  E75.1.bin50.celltype.immune.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% 
                                                              c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E75.1.bin50.celltype.immune.df$cell_type <- factor(E75.1.bin50.celltype.immune.df$cell_type, 
                                                     levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E75.1.color <- c(rep("#D2691E", 7), decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin50-immune-C-rect.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin50-immune-C.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E75.1.bin50.celltype.immune.df.annotation <- E75.1.bin50.celltype.immune.df[E75.1.bin50.celltype.immune.df$x >= 7000 & E75.1.bin50.celltype.immune.df$x <= 9000 & 
                                                                                E75.1.bin50.celltype.immune.df$y <= 15000 & E75.1.bin50.celltype.immune.df$y >= 13000, ]
  E75.1.bin50.celltype.df.annotation <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$x >= 7000 & E75.1.bin50.celltype.df$x <= 9000 & 
                                                                  E75.1.bin50.celltype.df$y <= 15000 & E75.1.bin50.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin50.celltype.immune.df$cell_type), unique(E75.1.bin50.celltype.immune.df.annotation$cell_type)), 
                                       levels(E75.1.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin50-immune-C-annotation.pdf", width = 2, height = 2)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin50
  E85.B6.bin50.celltype.immune.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% 
                                                                c("Mono","Mac","Pro-Mac","DC-1","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E85.B6.bin50.celltype.immune.df$cell_type <- factor(E85.B6.bin50.celltype.immune.df$cell_type, 
                                                      levels = c("Mono","Mac","Pro-Mac","DC-1","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E85.B6.color <- c(rep("#D2691E", 5), decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin50-immune-C-rect.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin50-immune-C.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin50.celltype.immune.df.annotation <- E85.B6.bin50.celltype.immune.df[E85.B6.bin50.celltype.immune.df$x >= 4000 & E85.B6.bin50.celltype.immune.df$x <= 6500 & 
                                                                                  E85.B6.bin50.celltype.immune.df$y >= -4000 & E85.B6.bin50.celltype.immune.df$y <= -2000, ]
  E85.B6.bin50.celltype.df.annotation <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$x >= 4000 & E85.B6.bin50.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin50.celltype.df$y >= -4000 & E85.B6.bin50.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin50.celltype.immune.df$cell_type), unique(E85.B6.bin50.celltype.immune.df.annotation$cell_type)), 
                                         levels(E85.B6.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin50-immune-C-annotation.pdf", width = 2.3, height = 1.9)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin50
  E95.bin50.celltype.immune.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% 
                                                          c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E95.bin50.celltype.immune.df$cell_type <- factor(E95.bin50.celltype.immune.df$cell_type, 
                                                   levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E95.color <- c(rep("#D2691E", 7), decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  E95.bin50.celltype.immune.eSF.df <- E95.bin50.celltype.immune.df[E95.bin50.celltype.immune.df$cell_type %in% "Lum-D", ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin50-immune-C-rect.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin50-immune-C.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    geom_point(E95.bin50.celltype.immune.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 1.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin50.celltype.immune.df.annotation <- E95.bin50.celltype.immune.df[E95.bin50.celltype.immune.df$x >= 11300 & E95.bin50.celltype.immune.df$x <= 13800 & 
                                                                            E95.bin50.celltype.immune.df$y >= 17000 & E95.bin50.celltype.immune.df$y <= 20000, ]
  E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df[E95.bin50.celltype.df$x >= 11300 & E95.bin50.celltype.df$x <= 13800 & 
                                                              E95.bin50.celltype.df$y >= 17000 & E95.bin50.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin50.celltype.immune.df$cell_type), unique(E95.bin50.celltype.immune.df.annotation$cell_type)), 
                                   levels(E95.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin50-immune-C-annotation.pdf", width = 2.3, height = 2.7)
  ggplot() + 
    geom_point(E95.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
}
#D. Immune subclusters + eSF
function(){
  #E65.2.bin50
  E65.2.bin50.celltype.immune.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% 
                                                              c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"), ]
  E65.2.bin50.celltype.immune.df$cell_type <- factor(E65.2.bin50.celltype.immune.df$cell_type, 
                                                     levels= c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"))
  E65.2.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E65-2-bin50-immune-eSF-D-rect.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E65-2-bin50-immune-eSF-D.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin50.celltype.immune.df.annotation <- E65.2.bin50.celltype.immune.df[E65.2.bin50.celltype.immune.df$x >= 7000 & E65.2.bin50.celltype.immune.df$x <= 8000 & 
                                                                                E65.2.bin50.celltype.immune.df$y >= 16700 & E65.2.bin50.celltype.immune.df$y <= 17700, ]
  E65.2.bin50.celltype.df.annotation <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$x >= 7000 & E65.2.bin50.celltype.df$x <= 8000 & 
                                                                  E65.2.bin50.celltype.df$y >= 16700 & E65.2.bin50.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin50.celltype.immune.df$cell_type), unique(E65.2.bin50.celltype.immune.df.annotation$cell_type)), 
                                       levels(E65.2.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E65-2-bin50-immune-eSF-D-annotation.pdf", width = 1.2, height = 1.1)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E75.1.bin50
  E75.1.bin50.celltype.immune.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% 
                                                              c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"), ]
  E75.1.bin50.celltype.immune.df$cell_type <- factor(E75.1.bin50.celltype.immune.df$cell_type, 
                                                     levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"))
  E75.1.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E75-1-bin50-immune-eSF-D-rect.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E75-1-bin50-immune-eSF-D.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E75.1.bin50.celltype.immune.df.annotation <- E75.1.bin50.celltype.immune.df[E75.1.bin50.celltype.immune.df$x >= 7000 & E75.1.bin50.celltype.immune.df$x <= 9000 & 
                                                                                E75.1.bin50.celltype.immune.df$y <= 15000 & E75.1.bin50.celltype.immune.df$y >= 13000, ]
  E75.1.bin50.celltype.df.annotation <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$x >= 7000 & E75.1.bin50.celltype.df$x <= 9000 & 
                                                                  E75.1.bin50.celltype.df$y <= 15000 & E75.1.bin50.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin50.celltype.immune.df$cell_type), unique(E75.1.bin50.celltype.immune.df.annotation$cell_type)), 
                                       levels(E75.1.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E75-1-bin50-immune-eSF-D-annotation.pdf", width = 2, height = 2)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin50
  E85.B6.bin50.celltype.immune.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% 
                                                                c("Mono","Mac","Pro-Mac","DC-1","NK", "Lum-D"), ]
  E85.B6.bin50.celltype.immune.df$cell_type <- factor(E85.B6.bin50.celltype.immune.df$cell_type, 
                                                      levels = c("Mono","Mac","Pro-Mac","DC-1","NK", "Lum-D"))
  E85.B6.color <- c(immune.color[c(1:4,7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E85.B6-2-bin50-immune-eSF-D-rect.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E85.B6-2-bin50-immune-eSF-D.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin50.celltype.immune.df.annotation <- E85.B6.bin50.celltype.immune.df[E85.B6.bin50.celltype.immune.df$x >= 4000 & E85.B6.bin50.celltype.immune.df$x <= 6500 & 
                                                                                  E85.B6.bin50.celltype.immune.df$y >= -4000 & E85.B6.bin50.celltype.immune.df$y <= -2000, ]
  E85.B6.bin50.celltype.df.annotation <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$x >= 4000 & E85.B6.bin50.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin50.celltype.df$y >= -4000 & E85.B6.bin50.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin50.celltype.immune.df$cell_type), unique(E85.B6.bin50.celltype.immune.df.annotation$cell_type)), 
                                         levels(E85.B6.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E85.B6-2-bin50-immune-eSF-D-annotation.pdf", width = 2.3, height = 1.9)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin50
  E95.bin50.celltype.immune.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% 
                                                          c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"), ]
  E95.bin50.celltype.immune.df$cell_type <- factor(E95.bin50.celltype.immune.df$cell_type, 
                                                   levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"))
  E95.bin50.celltype.immune.eSF.df <- E95.bin50.celltype.immune.df[E95.bin50.celltype.immune.df$cell_type %in% "Lum-D", ]
  E95.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E95-2-bin50-immune-eSF-D-rect.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E95-2-bin50-immune-eSF-D.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    geom_point(E95.bin50.celltype.immune.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 1.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin50.celltype.immune.df.annotation <- E95.bin50.celltype.immune.df[E95.bin50.celltype.immune.df$x >= 11300 & E95.bin50.celltype.immune.df$x <= 13800 & 
                                                                            E95.bin50.celltype.immune.df$y >= 17000 & E95.bin50.celltype.immune.df$y <= 20000, ]
  E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df[E95.bin50.celltype.df$x >= 11300 & E95.bin50.celltype.df$x <= 13800 & 
                                                              E95.bin50.celltype.df$y >= 17000 & E95.bin50.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin50.celltype.immune.df$cell_type), unique(E95.bin50.celltype.immune.df.annotation$cell_type)), 
                                   levels(E95.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/05.E95-2-bin50-immune-eSF-D-annotation.pdf", width = 2.3, height = 2.7)
  ggplot() + 
    geom_point(E95.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
}
#E. iDSC subclusters + eSF
function(){
  #E65.2.bin50
  E65.2.bin50.celltype.immune.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% 
                                                              c("Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E65.2.bin50.celltype.immune.df$cell_type <- factor(E65.2.bin50.celltype.immune.df$cell_type, 
                                                     levels= c("Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E65.2.color <- c(decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E65-2-bin50-iDSC-eSF-rect.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E65-2-bin50-iDSC-eSF-immune.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin50.celltype.immune.df.annotation <- E65.2.bin50.celltype.immune.df[E65.2.bin50.celltype.immune.df$x >= 7000 & E65.2.bin50.celltype.immune.df$x <= 8000 & 
                                                                                E65.2.bin50.celltype.immune.df$y >= 16700 & E65.2.bin50.celltype.immune.df$y <= 17700, ]
  E65.2.bin50.celltype.df.annotation <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$x >= 7000 & E65.2.bin50.celltype.df$x <= 8000 & 
                                                                  E65.2.bin50.celltype.df$y >= 16700 & E65.2.bin50.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin50.celltype.immune.df$cell_type), unique(E65.2.bin50.celltype.immune.df.annotation$cell_type)), 
                                       levels(E65.2.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E65-2-bin50-immune-iDSC-eSF-annotation.pdf", width = 1.2, height = 1.1)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E75.1.bin50
  E75.1.bin50.celltype.immune.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% 
                                                              c("Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E75.1.bin50.celltype.immune.df$cell_type <- factor(E75.1.bin50.celltype.immune.df$cell_type, 
                                                     levels = c("Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E75.1.color <- c(decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E75-1-bin50-immune-iDSC-eSF-rect.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E75-1-bin50-immune-iDSC-eSF.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E75.1.bin50.celltype.immune.df.annotation <- E75.1.bin50.celltype.immune.df[E75.1.bin50.celltype.immune.df$x >= 7000 & E75.1.bin50.celltype.immune.df$x <= 9000 & 
                                                                                E75.1.bin50.celltype.immune.df$y <= 15000 & E75.1.bin50.celltype.immune.df$y >= 13000, ]
  E75.1.bin50.celltype.df.annotation <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$x >= 7000 & E75.1.bin50.celltype.df$x <= 9000 & 
                                                                  E75.1.bin50.celltype.df$y <= 15000 & E75.1.bin50.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin50.celltype.immune.df$cell_type), unique(E75.1.bin50.celltype.immune.df.annotation$cell_type)), 
                                       levels(E75.1.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E75-1-bin50-immune-iDSC-eSF-annotation.pdf", width = 2, height = 2)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin50
  E85.B6.bin50.celltype.immune.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% 
                                                                c("Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E85.B6.bin50.celltype.immune.df$cell_type <- factor(E85.B6.bin50.celltype.immune.df$cell_type, 
                                                      levels = c("Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E85.B6.color <- c(decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E85.B6-2-bin50-immune-iDSC-eSF-rect.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E85.B6-2-bin50-immune-iDSC-eSF.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin50.celltype.immune.df.annotation <- E85.B6.bin50.celltype.immune.df[E85.B6.bin50.celltype.immune.df$x >= 4000 & E85.B6.bin50.celltype.immune.df$x <= 6500 & 
                                                                                  E85.B6.bin50.celltype.immune.df$y >= -4000 & E85.B6.bin50.celltype.immune.df$y <= -2000, ]
  E85.B6.bin50.celltype.df.annotation <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$x >= 4000 & E85.B6.bin50.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin50.celltype.df$y >= -4000 & E85.B6.bin50.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin50.celltype.immune.df$cell_type), unique(E85.B6.bin50.celltype.immune.df.annotation$cell_type)), 
                                         levels(E85.B6.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E85.B6-2-bin50-immune-iDSC-eSF-annotation.pdf", width = 2.3, height = 1.9)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin50
  E95.bin50.celltype.immune.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% 
                                                          c("Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E95.bin50.celltype.immune.df$cell_type <- factor(E95.bin50.celltype.immune.df$cell_type, 
                                                   levels = c("Lum-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E95.bin50.celltype.immune.eSF.df <- E95.bin50.celltype.immune.df[E95.bin50.celltype.immune.df$cell_type %in% "Lum-D", ]
  E95.color <- c(decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E95-2-bin50-immune-iDSC-eSF-rect.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E95-2-bin50-immune-iDSC-eSF.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    geom_point(E95.bin50.celltype.immune.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 1.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin50.celltype.immune.df.annotation <- E95.bin50.celltype.immune.df[E95.bin50.celltype.immune.df$x >= 11300 & E95.bin50.celltype.immune.df$x <= 13800 & 
                                                                            E95.bin50.celltype.immune.df$y >= 17000 & E95.bin50.celltype.immune.df$y <= 20000, ]
  E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df[E95.bin50.celltype.df$x >= 11300 & E95.bin50.celltype.df$x <= 13800 & 
                                                              E95.bin50.celltype.df$y >= 17000 & E95.bin50.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin50.celltype.immune.df$cell_type), unique(E95.bin50.celltype.immune.df.annotation$cell_type)), 
                                   levels(E95.bin50.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E95-2-bin50-immune-iDSC-eSF-annotation.pdf", width = 2.3, height = 2.7)
  ggplot() + 
    geom_point(E95.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
}

##################################################### Normal - bin50 EC subcluster distribution and zoom in#################################################
#A. EC suclusters + D6 + iDSC subclusters
function(){
  #E65.2.bin50
  E65.2.bin50.celltype.EC.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% 
                                                          c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E65.2.bin50.celltype.EC.df$cell_type <- factor(E65.2.bin50.celltype.EC.df$cell_type, 
                                                 levels= c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E65.2.color <- c(EC.color[1:5], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E65-2-bin50-iDSC-D6-rect.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/07.E65-2-bin50-iDSC-D6-EC.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin50.celltype.EC.df.annotation <- E65.2.bin50.celltype.EC.df[E65.2.bin50.celltype.EC.df$x >= 7000 & E65.2.bin50.celltype.EC.df$x <= 8000 & 
                                                                        E65.2.bin50.celltype.EC.df$y >= 16700 & E65.2.bin50.celltype.EC.df$y <= 17700, ]
  E65.2.bin50.celltype.df.annotation <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$x >= 7000 & E65.2.bin50.celltype.df$x <= 8000 & 
                                                                  E65.2.bin50.celltype.df$y >= 16700 & E65.2.bin50.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin50.celltype.EC.df$cell_type), unique(E65.2.bin50.celltype.EC.df.annotation$cell_type)), 
                                       levels(E65.2.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E65-2-bin50-EC-iDSC-eSF-annotation.pdf", width = 1.2, height = 1.1)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E75.1.bin50
  E75.1.bin50.celltype.EC.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% 
                                                          c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E75.1.bin50.celltype.EC.df$cell_type <- factor(E75.1.bin50.celltype.EC.df$cell_type, 
                                                 levels = c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E75.1.color <- c(EC.color[c(1,3:5)], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E75-1-bin50-EC-iDSC-eSF-rect.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/07.E75-1-bin50-iDSC-D6-EC.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E75.1.bin50.celltype.EC.df.annotation <- E75.1.bin50.celltype.EC.df[E75.1.bin50.celltype.EC.df$x >= 7000 & E75.1.bin50.celltype.EC.df$x <= 9000 & 
                                                                        E75.1.bin50.celltype.EC.df$y <= 15000 & E75.1.bin50.celltype.EC.df$y >= 13000, ]
  E75.1.bin50.celltype.df.annotation <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$x >= 7000 & E75.1.bin50.celltype.df$x <= 9000 & 
                                                                  E75.1.bin50.celltype.df$y <= 15000 & E75.1.bin50.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin50.celltype.EC.df$cell_type), unique(E75.1.bin50.celltype.EC.df.annotation$cell_type)), 
                                       levels(E75.1.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E75-1-bin50-EC-iDSC-eSF-annotation.pdf", width = 2, height = 2)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin50
  E85.B6.bin50.celltype.EC.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% 
                                                            c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E85.B6.bin50.celltype.EC.df$cell_type <- factor(E85.B6.bin50.celltype.EC.df$cell_type, 
                                                  levels = c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E85.B6.color <- c(EC.color[1:4], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E85.B6-2-bin50-EC-iDSC-eSF-rect.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/07.E85.B6-2-bin50-iDSC-D6-EC.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin50.celltype.EC.df.annotation <- E85.B6.bin50.celltype.EC.df[E85.B6.bin50.celltype.EC.df$x >= 4000 & E85.B6.bin50.celltype.EC.df$x <= 6500 & 
                                                                          E85.B6.bin50.celltype.EC.df$y >= -4000 & E85.B6.bin50.celltype.EC.df$y <= -2000, ]
  E85.B6.bin50.celltype.df.annotation <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$x >= 4000 & E85.B6.bin50.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin50.celltype.df$y >= -4000 & E85.B6.bin50.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin50.celltype.EC.df$cell_type), unique(E85.B6.bin50.celltype.EC.df.annotation$cell_type)), 
                                         levels(E85.B6.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E85.B6-2-bin50-EC-iDSC-eSF-annotation.pdf", width = 2.3, height = 1.9)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin50
  E95.bin50.celltype.EC.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% 
                                                      c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E95.bin50.celltype.EC.df$cell_type <- factor(E95.bin50.celltype.EC.df$cell_type, 
                                               levels = c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E95.color <- c(EC.color[1:5], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E95-2-bin50-EC-iDSC-eSF-rect.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/07.E95-2-bin50-iDSC-D6-EC.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin50.celltype.EC.df.annotation <- E95.bin50.celltype.EC.df[E95.bin50.celltype.EC.df$x >= 11300 & E95.bin50.celltype.EC.df$x <= 13800 & 
                                                                    E95.bin50.celltype.EC.df$y >= 17000 & E95.bin50.celltype.EC.df$y <= 20000, ]
  E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df[E95.bin50.celltype.df$x >= 11300 & E95.bin50.celltype.df$x <= 13800 & 
                                                              E95.bin50.celltype.df$y >= 17000 & E95.bin50.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin50.celltype.EC.df$cell_type), unique(E95.bin50.celltype.EC.df.annotation$cell_type)), 
                                   levels(E95.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E95-2-bin50-EC-iDSC-eSF-annotation.pdf", width = 2.3, height = 2.7)
  ggplot() + 
    geom_point(E95.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
}
#B. EC subclusters + D6
function(){
  #E65.2.bin50
  E65.2.bin50.celltype.EC.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% 
                                                          c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D"), ]
  E65.2.bin50.celltype.EC.df$cell_type <- factor(E65.2.bin50.celltype.EC.df$cell_type, 
                                                 levels= c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D"))
  E65.2.color <- c(EC.color[1:5], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E65-2-bin50-D6-rect.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/08.E65-2-bin50-D6-EC.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin50.celltype.EC.df.annotation <- E65.2.bin50.celltype.EC.df[E65.2.bin50.celltype.EC.df$x >= 7000 & E65.2.bin50.celltype.EC.df$x <= 8000 & 
                                                                        E65.2.bin50.celltype.EC.df$y >= 16700 & E65.2.bin50.celltype.EC.df$y <= 17700, ]
  E65.2.bin50.celltype.df.annotation <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$x >= 7000 & E65.2.bin50.celltype.df$x <= 8000 & 
                                                                  E65.2.bin50.celltype.df$y >= 16700 & E65.2.bin50.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin50.celltype.EC.df$cell_type), unique(E65.2.bin50.celltype.EC.df.annotation$cell_type)), 
                                       levels(E65.2.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E65-2-bin50-EC-eSF-annotation.pdf", width = 1.2, height = 1.1)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E75.1.bin50
  E75.1.bin50.celltype.EC.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% 
                                                          c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D"), ]
  E75.1.bin50.celltype.EC.df$cell_type <- factor(E75.1.bin50.celltype.EC.df$cell_type, 
                                                 levels = c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D"))
  E75.1.color <- c(EC.color[c(1,3:5)], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E75-1-bin50-EC-eSF-rect.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/08.E75-1-bin50-D6-EC.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E75.1.bin50.celltype.EC.df.annotation <- E75.1.bin50.celltype.EC.df[E75.1.bin50.celltype.EC.df$x >= 7000 & E75.1.bin50.celltype.EC.df$x <= 9000 & 
                                                                        E75.1.bin50.celltype.EC.df$y <= 15000 & E75.1.bin50.celltype.EC.df$y >= 13000, ]
  E75.1.bin50.celltype.df.annotation <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$x >= 7000 & E75.1.bin50.celltype.df$x <= 9000 & 
                                                                  E75.1.bin50.celltype.df$y <= 15000 & E75.1.bin50.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin50.celltype.EC.df$cell_type), unique(E75.1.bin50.celltype.EC.df.annotation$cell_type)), 
                                       levels(E75.1.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E75-1-bin50-EC-eSF-annotation.pdf", width = 2, height = 2)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin50
  E85.B6.bin50.celltype.EC.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% 
                                                            c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D"), ]
  E85.B6.bin50.celltype.EC.df$cell_type <- factor(E85.B6.bin50.celltype.EC.df$cell_type, 
                                                  levels = c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D"))
  E85.B6.color <- c(EC.color[1:4], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E85.B6-2-bin50-EC-eSF-rect.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/08.E85.B6-2-bin50-D6-EC.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin50.celltype.EC.df.annotation <- E85.B6.bin50.celltype.EC.df[E85.B6.bin50.celltype.EC.df$x >= 4000 & E85.B6.bin50.celltype.EC.df$x <= 6500 & 
                                                                          E85.B6.bin50.celltype.EC.df$y >= -4000 & E85.B6.bin50.celltype.EC.df$y <= -2000, ]
  E85.B6.bin50.celltype.df.annotation <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$x >= 4000 & E85.B6.bin50.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin50.celltype.df$y >= -4000 & E85.B6.bin50.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin50.celltype.EC.df$cell_type), unique(E85.B6.bin50.celltype.EC.df.annotation$cell_type)), 
                                         levels(E85.B6.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E85.B6-2-bin50-EC-eSF-annotation.pdf", width = 2.3, height = 1.9)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin50
  E95.bin50.celltype.EC.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% 
                                                      c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D"), ]
  E95.bin50.celltype.EC.df$cell_type <- factor(E95.bin50.celltype.EC.df$cell_type, 
                                               levels = c("EC-0", "EC-1", "EC-2", "EC-3", "EC-4", "S100a8-D"))
  E95.color <- c(EC.color[1:5], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E95-2-bin50-EC-eSF-rect.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/08.E95-2-bin50-D6-EC.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin50.celltype.EC.df.annotation <- E95.bin50.celltype.EC.df[E95.bin50.celltype.EC.df$x >= 11300 & E95.bin50.celltype.EC.df$x <= 13800 & 
                                                                    E95.bin50.celltype.EC.df$y >= 17000 & E95.bin50.celltype.EC.df$y <= 20000, ]
  E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df[E95.bin50.celltype.df$x >= 11300 & E95.bin50.celltype.df$x <= 13800 & 
                                                              E95.bin50.celltype.df$y >= 17000 & E95.bin50.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin50.celltype.EC.df$cell_type), unique(E95.bin50.celltype.EC.df.annotation$cell_type)), 
                                   levels(E95.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E95-2-bin50-EC-eSF-annotation.pdf", width = 2.3, height = 2.7)
  ggplot() + 
    geom_point(E95.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
}
#C. D6 + iDSC subclusters
function(){
  #E65.2.bin50
  E65.2.bin50.celltype.EC.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% 
                                                          c("S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E65.2.bin50.celltype.EC.df$cell_type <- factor(E65.2.bin50.celltype.EC.df$cell_type, 
                                                 levels= c("S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E65.2.color <- c(decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E65-2-bin50-iDSC-D6-rect.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/09.E65-2-bin50-iDSC-D6.pdf", width = 4.5, height = 4)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin50.celltype.EC.df.annotation <- E65.2.bin50.celltype.EC.df[E65.2.bin50.celltype.EC.df$x >= 7000 & E65.2.bin50.celltype.EC.df$x <= 8000 & 
                                                                        E65.2.bin50.celltype.EC.df$y >= 16700 & E65.2.bin50.celltype.EC.df$y <= 17700, ]
  E65.2.bin50.celltype.df.annotation <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$x >= 7000 & E65.2.bin50.celltype.df$x <= 8000 & 
                                                                  E65.2.bin50.celltype.df$y >= 16700 & E65.2.bin50.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin50.celltype.EC.df$cell_type), unique(E65.2.bin50.celltype.EC.df.annotation$cell_type)), 
                                       levels(E65.2.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E65-2-bin50-iDSC-eSF-annotation.pdf", width = 1.2, height = 1.1)
  ggplot() + 
    geom_point(E65.2.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E75.1.bin50
  E75.1.bin50.celltype.EC.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% 
                                                          c("S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E75.1.bin50.celltype.EC.df$cell_type <- factor(E75.1.bin50.celltype.EC.df$cell_type, 
                                                 levels = c("S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E75.1.color <- c(decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E75-1-bin50-iDSC-eSF-rect.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/09.E75-1-bin50-iDSC-D6.pdf", width = 5.6, height = 5)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E75.1.bin50.celltype.EC.df.annotation <- E75.1.bin50.celltype.EC.df[E75.1.bin50.celltype.EC.df$x >= 7000 & E75.1.bin50.celltype.EC.df$x <= 9000 & 
                                                                        E75.1.bin50.celltype.EC.df$y <= 15000 & E75.1.bin50.celltype.EC.df$y >= 13000, ]
  E75.1.bin50.celltype.df.annotation <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$x >= 7000 & E75.1.bin50.celltype.df$x <= 9000 & 
                                                                  E75.1.bin50.celltype.df$y <= 15000 & E75.1.bin50.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin50.celltype.EC.df$cell_type), unique(E75.1.bin50.celltype.EC.df.annotation$cell_type)), 
                                       levels(E75.1.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E75-1-bin50-iDSC-eSF-annotation.pdf", width = 2, height = 2)
  ggplot() + 
    geom_point(E75.1.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin50
  E85.B6.bin50.celltype.EC.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% 
                                                            c("S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E85.B6.bin50.celltype.EC.df$cell_type <- factor(E85.B6.bin50.celltype.EC.df$cell_type, 
                                                  levels = c("S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E85.B6.color <- c(decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E85.B6-2-bin50-iDSC-eSF-rect.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/09.E85.B6-2-bin50-iDSC-D6.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin50.celltype.EC.df.annotation <- E85.B6.bin50.celltype.EC.df[E85.B6.bin50.celltype.EC.df$x >= 4000 & E85.B6.bin50.celltype.EC.df$x <= 6500 & 
                                                                          E85.B6.bin50.celltype.EC.df$y >= -4000 & E85.B6.bin50.celltype.EC.df$y <= -2000, ]
  E85.B6.bin50.celltype.df.annotation <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$x >= 4000 & E85.B6.bin50.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin50.celltype.df$y >= -4000 & E85.B6.bin50.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin50.celltype.EC.df$cell_type), unique(E85.B6.bin50.celltype.EC.df.annotation$cell_type)), 
                                         levels(E85.B6.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E85.B6-2-bin50-iDSC-eSF-annotation.pdf", width = 2.3, height = 1.9)
  ggplot() + 
    geom_point(E85.B6.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin50
  E95.bin50.celltype.EC.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% 
                                                      c("S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"), ]
  E95.bin50.celltype.EC.df$cell_type <- factor(E95.bin50.celltype.EC.df$cell_type, 
                                               levels = c("S100a8-D", "FB_immune_0","FB_immune_1", "FB_immune_2"))
  E95.color <- c(decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E95-2-bin50-iDSC-eSF-rect.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/09.E95-2-bin50-iDSC-D6.pdf", width = 10, height = 8)
  ggplot() + 
    geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin50.celltype.EC.df.annotation <- E95.bin50.celltype.EC.df[E95.bin50.celltype.EC.df$x >= 11300 & E95.bin50.celltype.EC.df$x <= 13800 & 
                                                                    E95.bin50.celltype.EC.df$y >= 17000 & E95.bin50.celltype.EC.df$y <= 20000, ]
  E95.bin50.celltype.df.annotation <- E95.bin50.celltype.df[E95.bin50.celltype.df$x >= 11300 & E95.bin50.celltype.df$x <= 13800 & 
                                                              E95.bin50.celltype.df$y >= 17000 & E95.bin50.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin50.celltype.EC.df$cell_type), unique(E95.bin50.celltype.EC.df.annotation$cell_type)), 
                                   levels(E95.bin50.celltype.EC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/06.E95-2-bin50-iDSC-eSF-annotation.pdf", width = 2.3, height = 2.7)
  ggplot() + 
    geom_point(E95.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin50.celltype.EC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
}


##################################################### CBA - bin50 immune subcluster distribution and zoom in################################################
#read files
DBA.E3.bin50.iDSC.subcluster.df <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-iDSC-subclusters-celltype-20230306.csv", sep = ",", check.names = F)
DBA.E6.bin50.iDSC.subcluster.df <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-iDSC-subclusters-celltype-20230306.csv", sep = ",", check.names = F)
DBA.F1.bin50.iDSC.subcluster.df <- read.table("20220622-Comments/05.DBA/DBA-F1-bin50-iDSC-subclusters-celltype-20230306.csv", sep = ",", check.names = F)

#A. Immune subclusters + eSF
function(){
  #E3
  DBA.E3.bin50.celltype.immune.df <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% 
                                                              c("Monocyte","Mac","Proliferating Mac","Neutrophil","NK", "D1-eSF"), ]
  DBA.E3.bin50.celltype.immune.df$cell_type <- factor(DBA.E3.bin50.celltype.immune.df$cell_type, 
                                                     levels= c("Monocyte","Mac","Proliferating Mac","Neutrophil","NK", "D1-eSF"))
  DBA.E3.color <- c(immune.color[c(1:3,6:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  DBA.E3.bin50.celltype.eSF.df <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/10.E3-bin50-eSF-immune.pdf", width = 6.5, height = 6.5)
  ggplot() + 
    geom_point(DBA.E3.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E3.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    geom_point(DBA.E3.bin50.celltype.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 1) + 
    scale_color_manual(values = DBA.E3.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  #E6
  DBA.E6.bin50.celltype.immune.df <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("Monocyte","Mac","Proliferating Mac","Neutrophil","NK", "D1-eSF"), ]
  DBA.E6.bin50.celltype.immune.df$cell_type <- factor(DBA.E6.bin50.celltype.immune.df$cell_type, 
                                                      levels= c("Monocyte","Mac","Proliferating Mac","Neutrophil","NK", "D1-eSF"))
  DBA.E6.color <- c(immune.color[c(1:3,6:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  DBA.E6.bin50.celltype.eSF.df <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/10.E6-bin50-eSF-immune.pdf", width = 9, height = 7)
  ggplot() + 
    geom_point(DBA.E6.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    geom_point(DBA.E6.bin50.celltype.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 1) + 
    scale_color_manual(values = DBA.E6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  #F1
  DBA.F1.bin50.celltype.immune.df <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("Monocyte","Mac","Proliferating Mac","Neutrophil","NK", "D1-eSF"), ]
  DBA.F1.bin50.celltype.immune.df$cell_type <- factor(DBA.F1.bin50.celltype.immune.df$cell_type, 
                                                      levels= c("Monocyte","Mac","Proliferating Mac","Neutrophil","NK", "D1-eSF"))
  DBA.F1.color <- c(immune.color[c(1:3,6:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  DBA.F1.bin50.celltype.eSF.df <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/10.F1-bin50-eSF-immune.pdf", width = 6.5, height = 6)
  ggplot() + 
    geom_point(DBA.F1.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.F1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = DBA.F1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
}
#B. iDSC subclusters + eSF
function(){
  #E3
  DBA.E3.bin50.celltype.immune.df <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D1-eSF"), ]
  DBA.E3.bin50.celltype.immune.df$cell_type <- factor(DBA.E3.bin50.celltype.immune.df$cell_type, 
                                                      levels= c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D1-eSF"))
  DBA.E3.color <- c("#F8766D", "#00BA38", "#619CFF", decidual.color[1])
  DBA.E3.bin50.celltype.eSF.df <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/11.E3-bin50-eSF-iDSC.pdf", width = 6.5, height = 6.5)
  ggplot() + 
    geom_point(DBA.E3.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E3.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    geom_point(DBA.E3.bin50.celltype.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 1) + 
    scale_color_manual(values = DBA.E3.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  #E6
  DBA.E6.bin50.celltype.immune.df <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D1-eSF"), ]
  DBA.E6.bin50.celltype.immune.df$cell_type <- factor(DBA.E6.bin50.celltype.immune.df$cell_type, 
                                                      levels= c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D1-eSF"))
  DBA.E6.color <- c("#F8766D", "#00BA38", "#619CFF", decidual.color[1])
  DBA.E6.bin50.celltype.eSF.df <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/11.E6-bin50-eSF-iDSC.pdf", width = 9, height = 7)
  ggplot() + 
    geom_point(DBA.E6.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E6.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    geom_point(DBA.E6.bin50.celltype.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 1) + 
    scale_color_manual(values = DBA.E6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  #F1
  DBA.F1.bin50.celltype.immune.df <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D1-eSF"), ]
  DBA.F1.bin50.celltype.immune.df$cell_type <- factor(DBA.F1.bin50.celltype.immune.df$cell_type, 
                                                      levels= c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D1-eSF"))
  DBA.F1.color <- c("#F8766D", "#00BA38", "#619CFF", decidual.color[1])
  DBA.F1.bin50.celltype.eSF.df <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/11.F1-bin50-eSF-iDSC.pdf", width = 6.5, height = 6)
  ggplot() + 
    geom_point(DBA.F1.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.F1.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = DBA.F1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
}

#####################################################CBA - bin50 EC subcluster distribution and zoom in#####################################################
#A. EC subclusters + D6
function(){
  #E3
  DBA.E3.bin50.celltype.EC.df <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("Angiogenic EC", "Proliterating EC", "Venous EC1", "Venous EC2", "Arterial EC", "D6-S100a8"), ]
  DBA.E3.bin50.celltype.EC.df$cell_type <- factor(DBA.E3.bin50.celltype.EC.df$cell_type, 
                                                      levels= c("Angiogenic EC", "Proliterating EC", "Venous EC1", "Venous EC2", "Arterial EC", "D6-S100a8"))
  DBA.E3.color <- c(EC.color[c(1,3,2,4)], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  DBA.E3.bin50.celltype.eSF.df <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% c("D6-S100a8"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/12.E3-bin50-D6-EC-subclusters.pdf", width = 6.5, height = 6.5)
  ggplot() + 
    geom_point(DBA.E3.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E3.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = DBA.E3.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  #E6
  DBA.E6.bin50.celltype.EC.df <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("Angiogenic EC", "Proliterating EC", "Venous EC1", "Venous EC2", "Arterial EC", "D6-S100a8"), ]
  DBA.E6.bin50.celltype.EC.df$cell_type <- factor(DBA.E6.bin50.celltype.EC.df$cell_type, 
                                                      levels= c("Angiogenic EC", "Proliterating EC", "Venous EC1", "Venous EC2", "Arterial EC", "D6-S100a8"))
  DBA.E6.color <- c(EC.color[c(1,3,2,4)], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  DBA.E6.bin50.celltype.eSF.df <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% c("D6-S100a8"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/12.E6-bin50-D6-EC-subclusters.pdf", width = 9, height = 7)
  ggplot() + 
    geom_point(DBA.E6.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E6.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) + 
    scale_color_manual(values = DBA.E6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  #F1
  DBA.F1.bin50.celltype.EC.df <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("Angiogenic EC", "Proliterating EC", "Venous EC1", "Venous EC2", "Arterial EC", "D6-S100a8"), ]
  DBA.F1.bin50.celltype.EC.df$cell_type <- factor(DBA.F1.bin50.celltype.EC.df$cell_type, 
                                                      levels= c("Angiogenic EC", "Proliterating EC", "Venous EC1", "Venous EC2", "Arterial EC", "D6-S100a8"))
  DBA.F1.color <- c(EC.color[c(1,3,2,4)], decidual.color[6], "#F8766D", "#00BA38", "#619CFF")
  DBA.F1.bin50.celltype.eSF.df <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% c("D6-S100a8"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/12.F1-bin50-D6-EC-subclusters.pdf", width = 6.5, height = 6)
  ggplot() + 
    geom_point(DBA.F1.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.F1.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = DBA.F1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
}
#B. iDSC subclusters + D6
function(){
  #E3
  DBA.E3.bin50.celltype.EC.df <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D6-S100a8"), ]
  DBA.E3.bin50.celltype.EC.df$cell_type <- factor(DBA.E3.bin50.celltype.EC.df$cell_type, 
                                                      levels= c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D6-S100a8"))
  DBA.E3.color <- c("#F8766D", "#00BA38", "#619CFF", decidual.color[6])
  DBA.E3.bin50.celltype.S100a8.df <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% c("D6-S100a8"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/13.E3-bin50-S100a8-iDSC.pdf", width = 6.5, height = 6.5)
  ggplot() + 
    geom_point(DBA.E3.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E3.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) + 
    scale_color_manual(values = DBA.E3.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  #E6
  DBA.E6.bin50.celltype.EC.df <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D6-S100a8"), ]
  DBA.E6.bin50.celltype.EC.df$cell_type <- factor(DBA.E6.bin50.celltype.EC.df$cell_type, 
                                                      levels= c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D6-S100a8"))
  DBA.E6.color <- c("#F8766D", "#00BA38", "#619CFF", decidual.color[6])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/13.E6-bin50-S100a8-iDSC.pdf", width = 9, height = 7)
  ggplot() + 
    geom_point(DBA.E6.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E6.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) + 
    scale_color_manual(values = DBA.E6.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  #F1
  DBA.F1.bin50.celltype.EC.df <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% 
                                                                       c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D6-S100a8"), ]
  DBA.F1.bin50.celltype.EC.df$cell_type <- factor(DBA.F1.bin50.celltype.EC.df$cell_type, 
                                                      levels= c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2", "D6-S100a8"))
  DBA.F1.color <- c("#F8766D", "#00BA38", "#619CFF", decidual.color[6])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/13.F1-bin50-S100a8-iDSC.pdf", width = 6.5, height = 6)
  ggplot() + 
    geom_point(DBA.F1.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.F1.bin50.celltype.EC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = DBA.F1.color) + 
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
}

#####################################################bin20 immune cell distribution and zoom in#################################################
#A. Immune subclusters + iDSC subclusters + eSF and zoom in pics
function(){
  #E65.2.bin20
  E65.2.bin20.celltype.df <- read.table("20210827-figures/figure5/E65-2-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E65.2.bin50.celltype.df <- read.table("20210827-figures/figure5/E65-2-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E65.2.bin20.celltype.immune.df <- E65.2.bin20.celltype.df[E65.2.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E65.2.bin50.eSF.iDSC.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% c("Lum-D", "FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E65.2.bin20.celltype.immune.eSF.iDSC.df <- rbind(E65.2.bin20.celltype.immune.df, E65.2.bin50.eSF.iDSC.df)
  E65.2.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(E65.2.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                              levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E65.2.color <- c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E65-2-bin20-immune-eSF-iDSC.pdf", width =7.5, height = 6)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin20-immune-eSF-iDSC-rect.pdf", width = 7.5, height = 6)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin20.celltype.immune.eSF.iDSC.df.annotation <- E65.2.bin20.celltype.immune.eSF.iDSC.df[E65.2.bin20.celltype.immune.eSF.iDSC.df$x >= 7000 & E65.2.bin20.celltype.immune.eSF.iDSC.df$x <= 8000 & 
                                                                                E65.2.bin20.celltype.immune.eSF.iDSC.df$y >= 16700 & E65.2.bin20.celltype.immune.eSF.iDSC.df$y <= 17700, ]
  E65.2.bin20.celltype.df.annotation <- E65.2.bin20.celltype.df[E65.2.bin20.celltype.df$x >= 7000 & E65.2.bin20.celltype.df$x <= 8000 & 
                                                                  E65.2.bin20.celltype.df$y >= 16700 & E65.2.bin20.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin20.celltype.immune.eSF.iDSC.df$cell_type), unique(E65.2.bin20.celltype.immune.eSF.iDSC.df.annotation$cell_type)), 
                                       levels(E65.2.bin20.celltype.immune.eSF.iDSC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E65-2-bin20-immune-eSF-iDSC-annotation.pdf", width = 2.2, height = 2.)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.eSF.iDSC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()

  E75.1.bin20.celltype.df <- read.table("20210827-figures/figure5/E75-1-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E75.1.bin50.celltype.df <- read.table("20210827-figures/figure5/E75-1-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E75.1.bin20.celltype.immune.df <- E75.1.bin20.celltype.df[E75.1.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E75.1.bin50.celltype.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% c("Lum-D", "FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E75.1.bin20.celltype.immune.eSF.iDSC.df <- rbind(E75.1.bin20.celltype.immune.df, E75.1.bin50.celltype.df)
  E75.1.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(E75.1.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                              levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E75-1-bin20-immune-eSF-iDSC.pdf", width = 9, height = 8.5)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  E75.1.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin20-immune-eSF-iDSC-rect.pdf", width = 9, height = 8.5)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  
  E75.1.bin20.celltype.immune.eSF.iDSC.df.annotation <- E75.1.bin20.celltype.immune.eSF.iDSC.df[E75.1.bin20.celltype.immune.eSF.iDSC.df$x >= 7000 & E75.1.bin20.celltype.immune.eSF.iDSC.df$x <= 9000 & 
                                                                                E75.1.bin20.celltype.immune.eSF.iDSC.df$y <= 15000 & E75.1.bin20.celltype.immune.eSF.iDSC.df$y >= 13000, ]
  E75.1.bin20.celltype.df.annotation <- E75.1.bin20.celltype.df[E75.1.bin20.celltype.df$x >= 7000 & E75.1.bin20.celltype.df$x <= 9000 & 
                                                                  E75.1.bin20.celltype.df$y <= 15000 & E75.1.bin20.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin20.celltype.immune.eSF.iDSC.df$cell_type), unique(E75.1.bin20.celltype.immune.eSF.iDSC.df.annotation$cell_type)), 
                                       levels(E75.1.bin20.celltype.immune.eSF.iDSC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E75-1-bin20-immune-eSF-iDSC-annotation.pdf", width = 5.5, height = 4.1)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.eSF.iDSC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'right', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin20
  E85.B6.bin20.celltype.df <- read.table("20210827-figures/figure5/E85-B6-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E85.B6.bin50.celltype.df <- read.table("20210827-figures/figure5/E85-B6-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E85.B6.bin20.celltype.immune.df <- E85.B6.bin20.celltype.df[E85.B6.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E85.B6.bin50.eSF.iDSC.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% c("Lum-D", "FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E85.B6.bin20.celltype.immune.eSF.iDSC.df <- rbind(E85.B6.bin20.celltype.immune.df, E85.B6.bin50.eSF.iDSC.df)
  E85.B6.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(E85.B6.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                               levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E85-B6-bin20-immune-eSF-iDSC.pdf", width =10.5, height =10)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.color <- c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin20-immune-eSF-iDSC-rect.pdf", width = 10.5, height = 10)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin20.celltype.immune.eSF.iDSC.df.annotation <- E85.B6.bin20.celltype.immune.eSF.iDSC.df[E85.B6.bin20.celltype.immune.eSF.iDSC.df$x >= 4000 & E85.B6.bin20.celltype.immune.eSF.iDSC.df$x <= 6500 & 
                                                                                  E85.B6.bin20.celltype.immune.eSF.iDSC.df$y >= -4000 & E85.B6.bin20.celltype.immune.eSF.iDSC.df$y <= -2000, ]
  E85.B6.bin20.celltype.df.annotation <- E85.B6.bin20.celltype.df[E85.B6.bin20.celltype.df$x >= 4000 & E85.B6.bin20.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin20.celltype.df$y >= -4000 & E85.B6.bin20.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin20.celltype.immune.eSF.iDSC.df$cell_type), unique(E85.B6.bin20.celltype.immune.eSF.iDSC.df.annotation$cell_type)), 
                                         levels(E85.B6.bin20.celltype.immune.eSF.iDSC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E85.B6-2-bin20-immune-eSF-iDSC-annotation.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.eSF.iDSC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin20
  E95.bin20.celltype.df <- read.table("20210827-figures/figure5/E95-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E95.bin50.celltype.df <- read.table("20210827-figures/figure5/E95-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E95.bin20.celltype.immune.df <- E95.bin20.celltype.df[E95.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E95.bin50.eSF.iDSC.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% c("Lum-D", "FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E95.bin20.celltype.immune.eSF.iDSC.df <- rbind(E95.bin20.celltype.immune.df, E95.bin50.eSF.iDSC.df)
  E95.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(E95.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                            levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E95.bin20.celltype.immune.eSF.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% c("Lum-D"), ]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E95-bin20-immune-eSF-iDSC.pdf", width = 15, height = 13)
  ggplot() + 
    geom_point(E95.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    geom_point(E95.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 1.5) + 
    scale_color_manual(values = c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin20-immune-eSF-iDSC-rect.pdf", width = 15, height = 13)
  ggplot() + 
    geom_point(E95.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin20.celltype.immune.eSF.iDSC.df.annotation <- E95.bin20.celltype.immune.eSF.iDSC.df[E95.bin20.celltype.immune.eSF.iDSC.df$x >= 11300 & E95.bin20.celltype.immune.eSF.iDSC.df$x <= 13800 & 
                                                                            E95.bin20.celltype.immune.eSF.iDSC.df$y >= 17000 & E95.bin20.celltype.immune.eSF.iDSC.df$y <= 20000, ]
  E95.bin20.celltype.df.annotation <- E95.bin20.celltype.df[E95.bin20.celltype.df$x >= 11300 & E95.bin20.celltype.df$x <= 13800 & 
                                                              E95.bin20.celltype.df$y >= 17000 & E95.bin20.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin20.celltype.immune.eSF.iDSC.df$cell_type), unique(E95.bin20.celltype.immune.eSF.iDSC.df.annotation$cell_type)), 
                                   levels(E95.bin20.celltype.immune.eSF.iDSC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.E95-2-bin20-immune-eSF-iDSC-annotation.pdf", width = 5.2, height = 6.2)
  ggplot() + 
    geom_point(E95.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.eSF.iDSC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
}
#B. Immune subclusters and zoom in pics
function(){
  #E65.2.bin20
  E65.2.bin20.celltype.df <- read.table("20210827-figures/figure5/E65-2-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E65.2.bin20.celltype.immune.df <- E65.2.bin20.celltype.df[E65.2.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E65.2.bin20.celltype.immune.df$cell_type <- factor(E65.2.bin20.celltype.immune.df$cell_type, 
                                                              levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"))
  E65.2.color <- c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E65-2-bin20-immune-B.pdf", width =7.5, height = 6)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E65-2-bin20-immune-B-rect.pdf", width = 7.5, height = 6)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin20.celltype.immune.df.annotation <- E65.2.bin20.celltype.immune.df[E65.2.bin20.celltype.immune.df$x >= 7000 & E65.2.bin20.celltype.immune.df$x <= 8000 & 
                                                                                                  E65.2.bin20.celltype.immune.df$y >= 16700 & E65.2.bin20.celltype.immune.df$y <= 17700, ]
  E65.2.bin20.celltype.df.annotation <- E65.2.bin20.celltype.df[E65.2.bin20.celltype.df$x >= 7000 & E65.2.bin20.celltype.df$x <= 8000 & 
                                                                  E65.2.bin20.celltype.df$y >= 16700 & E65.2.bin20.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin20.celltype.immune.df$cell_type), unique(E65.2.bin20.celltype.immune.df.annotation$cell_type)), 
                                       levels(E65.2.bin20.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E65-2-bin20-immune-B-annotation.pdf", width = 2.2, height = 2.)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  E75.1.bin20.celltype.df <- read.table("20210827-figures/figure5/E75-1-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E75.1.bin20.celltype.immune.df <- E75.1.bin20.celltype.df[E75.1.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E75.1.bin20.celltype.immune.df$cell_type <- factor(E75.1.bin20.celltype.immune.df$cell_type, 
                                                              levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"))
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E75-1-bin20-immune-B.pdf", width = 9, height = 8.5)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  E75.1.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E75-1-bin20-immune-B-rect.pdf", width = 9, height = 8.5)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  
  E75.1.bin20.celltype.immune.df.annotation <- E75.1.bin20.celltype.immune.df[E75.1.bin20.celltype.immune.df$x >= 7000 & E75.1.bin20.celltype.immune.df$x <= 9000 & 
                                                                                                  E75.1.bin20.celltype.immune.df$y <= 15000 & E75.1.bin20.celltype.immune.df$y >= 13000, ]
  E75.1.bin20.celltype.df.annotation <- E75.1.bin20.celltype.df[E75.1.bin20.celltype.df$x >= 7000 & E75.1.bin20.celltype.df$x <= 9000 & 
                                                                  E75.1.bin20.celltype.df$y <= 15000 & E75.1.bin20.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin20.celltype.immune.df$cell_type), unique(E75.1.bin20.celltype.immune.df.annotation$cell_type)), 
                                       levels(E75.1.bin20.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E75-1-bin20-immune-B-annotation.pdf", width = 5.5, height = 4.1)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'right', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin20
  E85.B6.bin20.celltype.df <- read.table("20210827-figures/figure5/E85-B6-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E85.B6.bin20.celltype.immune.df <- E85.B6.bin20.celltype.df[E85.B6.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E85.B6.bin20.celltype.immune.df$cell_type <- factor(E85.B6.bin20.celltype.immune.df$cell_type, 
                                                               levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"))
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E85-B6-bin20-immune-B.pdf", width =10.5, height =10)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.color <- c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E85.B6-2-bin20-immune-B-rect.pdf", width = 10.5, height = 10)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin20.celltype.immune.df.annotation <- E85.B6.bin20.celltype.immune.df[E85.B6.bin20.celltype.immune.df$x >= 4000 & E85.B6.bin20.celltype.immune.df$x <= 6500 & 
                                                                                                    E85.B6.bin20.celltype.immune.df$y >= -4000 & E85.B6.bin20.celltype.immune.df$y <= -2000, ]
  E85.B6.bin20.celltype.df.annotation <- E85.B6.bin20.celltype.df[E85.B6.bin20.celltype.df$x >= 4000 & E85.B6.bin20.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin20.celltype.df$y >= -4000 & E85.B6.bin20.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin20.celltype.immune.df$cell_type), unique(E85.B6.bin20.celltype.immune.df.annotation$cell_type)), 
                                         levels(E85.B6.bin20.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E85.B6-2-bin20-immune-B-annotation.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin20
  E95.bin20.celltype.df <- read.table("20210827-figures/figure5/E95-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E95.bin20.celltype.immune.df <- E95.bin20.celltype.df[E95.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E95.bin20.celltype.immune.df$cell_type <- factor(E95.bin20.celltype.immune.df$cell_type, 
                                                            levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"))
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E95-bin20-immune-B.pdf", width = 15, height = 13)
  ggplot() + 
    geom_point(E95.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) + 
    scale_color_manual(values = c(immune.color[1:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  E95.color <- c(immune.color[c(1:7)], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E95-2-bin20-immune-B-rect.pdf", width = 15, height = 13)
  ggplot() + 
    geom_point(E95.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin20.celltype.immune.df.annotation <- E95.bin20.celltype.immune.df[E95.bin20.celltype.immune.df$x >= 11300 & E95.bin20.celltype.immune.df$x <= 13800 & 
                                                                                              E95.bin20.celltype.immune.df$y >= 17000 & E95.bin20.celltype.immune.df$y <= 20000, ]
  E95.bin20.celltype.df.annotation <- E95.bin20.celltype.df[E95.bin20.celltype.df$x >= 11300 & E95.bin20.celltype.df$x <= 13800 & 
                                                              E95.bin20.celltype.df$y >= 17000 & E95.bin20.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin20.celltype.immune.df$cell_type), unique(E95.bin20.celltype.immune.df.annotation$cell_type)), 
                                   levels(E95.bin20.celltype.immune.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E95-2-bin20-immune-B-annotation.pdf", width = 5.2, height = 6.2)
  ggplot() + 
    geom_point(E95.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
}
#C. Immune major cluster + iDSC subclusters + eSF and zoom in pics
function(){
  #E65.2.bin20
  E65.2.bin20.celltype.df <- read.table("20210827-figures/figure5/E65-2-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E65.2.bin50.celltype.df <- read.table("20210827-figures/figure5/E65-2-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E65.2.bin20.celltype.immune.df <- E65.2.bin20.celltype.df[E65.2.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E65.2.bin50.eSF.iDSC.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% c("Lum-D", "FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E65.2.bin20.celltype.immune.eSF.iDSC.df <- rbind(E65.2.bin20.celltype.immune.df, E65.2.bin50.eSF.iDSC.df)
  E65.2.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(E65.2.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                              levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E65.2.color <- c(rep("#D2691E", 7), decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E65-2-bin20-immune-eSF-iDSC-C.pdf", width =7.5, height = 6)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E65-2-bin20-immune-eSF-iDSC-rect-C.pdf", width = 7.5, height = 6)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 16700)) +
    geom_hline(aes(yintercept = 17700)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 8000)) +
    theme(axis.text = element_blank())
  dev.off()
  E65.2.bin20.celltype.immune.eSF.iDSC.df.annotation <- E65.2.bin20.celltype.immune.eSF.iDSC.df[E65.2.bin20.celltype.immune.eSF.iDSC.df$x >= 7000 & E65.2.bin20.celltype.immune.eSF.iDSC.df$x <= 8000 & 
                                                                                                  E65.2.bin20.celltype.immune.eSF.iDSC.df$y >= 16700 & E65.2.bin20.celltype.immune.eSF.iDSC.df$y <= 17700, ]
  E65.2.bin20.celltype.df.annotation <- E65.2.bin20.celltype.df[E65.2.bin20.celltype.df$x >= 7000 & E65.2.bin20.celltype.df$x <= 8000 & 
                                                                  E65.2.bin20.celltype.df$y >= 16700 & E65.2.bin20.celltype.df$y <= 17700, ]
  
  E65.2.ann.color <- E65.2.color[match(intersect(levels(E65.2.bin20.celltype.immune.eSF.iDSC.df$cell_type), unique(E65.2.bin20.celltype.immune.eSF.iDSC.df.annotation$cell_type)), 
                                       levels(E65.2.bin20.celltype.immune.eSF.iDSC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E65-2-bin20-immune-eSF-iDSC-annotation-C.pdf", width = 2.2, height = 2.)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.eSF.iDSC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  E75.1.bin20.celltype.df <- read.table("20210827-figures/figure5/E75-1-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E75.1.bin50.celltype.df <- read.table("20210827-figures/figure5/E75-1-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E75.1.bin20.celltype.immune.df <- E75.1.bin20.celltype.df[E75.1.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E75.1.bin50.celltype.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% c("Lum-D", "FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E75.1.bin20.celltype.immune.eSF.iDSC.df <- rbind(E75.1.bin20.celltype.immune.df, E75.1.bin50.celltype.df)
  E75.1.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(E75.1.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                              levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E75.1.color <- c(rep("#D2691E", 7), decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E75-1-bin20-immune-eSF-iDSC-C.pdf", width = 9, height = 8.5)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E75-1-bin20-immune-eSF-iDSC-rect-C.pdf", width = 9, height = 8.5)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) + 
    theme_classic() +
    geom_hline(aes(yintercept = 15000)) +
    geom_hline(aes(yintercept = 13000)) + 
    geom_vline(aes(xintercept = 7000)) +
    geom_vline(aes(xintercept = 9000)) +
    theme(axis.text = element_blank())
  dev.off()
  
  E75.1.bin20.celltype.immune.eSF.iDSC.df.annotation <- E75.1.bin20.celltype.immune.eSF.iDSC.df[E75.1.bin20.celltype.immune.eSF.iDSC.df$x >= 7000 & E75.1.bin20.celltype.immune.eSF.iDSC.df$x <= 9000 & 
                                                                                                  E75.1.bin20.celltype.immune.eSF.iDSC.df$y <= 15000 & E75.1.bin20.celltype.immune.eSF.iDSC.df$y >= 13000, ]
  E75.1.bin20.celltype.df.annotation <- E75.1.bin20.celltype.df[E75.1.bin20.celltype.df$x >= 7000 & E75.1.bin20.celltype.df$x <= 9000 & 
                                                                  E75.1.bin20.celltype.df$y <= 15000 & E75.1.bin20.celltype.df$y >= 13000, ]
  
  E75.1.ann.color <- E75.1.color[match(intersect(levels(E75.1.bin20.celltype.immune.eSF.iDSC.df$cell_type), unique(E75.1.bin20.celltype.immune.eSF.iDSC.df.annotation$cell_type)), 
                                       levels(E75.1.bin20.celltype.immune.eSF.iDSC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E75-1-bin20-immune-eSF-iDSC-annotation-C.pdf", width = 5.5, height = 4.1)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.eSF.iDSC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'right', axis.text = element_blank()) 
  dev.off()
  
  #E85.B6.bin20
  E85.B6.bin20.celltype.df <- read.table("20210827-figures/figure5/E85-B6-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E85.B6.bin50.celltype.df <- read.table("20210827-figures/figure5/E85-B6-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E85.B6.bin20.celltype.immune.df <- E85.B6.bin20.celltype.df[E85.B6.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E85.B6.bin50.eSF.iDSC.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% c("Lum-D", "FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E85.B6.bin20.celltype.immune.eSF.iDSC.df <- rbind(E85.B6.bin20.celltype.immune.df, E85.B6.bin50.eSF.iDSC.df)
  E85.B6.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(E85.B6.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                               levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E85.B6.color <- c(rep("#D2691E", 7), decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E85-B6-bin20-immune-eSF-iDSC-C.pdf", width =10.5, height =10)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E85.B6-2-bin20-immune-eSF-iDSC-rect-C.pdf", width = 10.5, height = 10)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 4000)) +
    geom_vline(aes(xintercept = 6500)) + 
    geom_hline(aes(yintercept = -2000)) +
    geom_hline(aes(yintercept = -4000)) +
    theme(axis.text = element_blank())
  dev.off()
  E85.B6.bin20.celltype.immune.eSF.iDSC.df.annotation <- E85.B6.bin20.celltype.immune.eSF.iDSC.df[E85.B6.bin20.celltype.immune.eSF.iDSC.df$x >= 4000 & E85.B6.bin20.celltype.immune.eSF.iDSC.df$x <= 6500 & 
                                                                                                    E85.B6.bin20.celltype.immune.eSF.iDSC.df$y >= -4000 & E85.B6.bin20.celltype.immune.eSF.iDSC.df$y <= -2000, ]
  E85.B6.bin20.celltype.df.annotation <- E85.B6.bin20.celltype.df[E85.B6.bin20.celltype.df$x >= 4000 & E85.B6.bin20.celltype.df$x <= 6500 & 
                                                                    E85.B6.bin20.celltype.df$y >= -4000 & E85.B6.bin20.celltype.df$y <= -2000, ]
  E85.B6.ann.color <- E85.B6.color[match(intersect(levels(E85.B6.bin20.celltype.immune.eSF.iDSC.df$cell_type), unique(E85.B6.bin20.celltype.immune.eSF.iDSC.df.annotation$cell_type)), 
                                         levels(E85.B6.bin20.celltype.immune.eSF.iDSC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E85.B6-2-bin20-immune-eSF-iDSC-annotation-C.pdf", width = 5, height = 4.5)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.eSF.iDSC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
  
  #E95.bin20
  E95.bin20.celltype.df <- read.table("20210827-figures/figure5/E95-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E95.bin50.celltype.df <- read.table("20210827-figures/figure5/E95-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E95.bin20.celltype.immune.df <- E95.bin20.celltype.df[E95.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E95.bin50.eSF.iDSC.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% c("Lum-D", "FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E95.bin20.celltype.immune.eSF.iDSC.df <- rbind(E95.bin20.celltype.immune.df, E95.bin50.eSF.iDSC.df)
  E95.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(E95.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                            levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E95.color <- c(rep("#D2691E", 7), decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E95-bin20-immune-eSF-iDSC-C.pdf", width = 15, height = 13)
  ggplot() + 
    geom_point(E95.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) + 
    scale_color_manual(values = c(rep("#D2691E", 7), decidual.color[1], "#F8766D", "#00BA38", "#619CFF")) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E95-2-bin20-immune-eSF-iDSC-rect-C.pdf", width = 15, height = 13)
  ggplot() + 
    geom_point(E95.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) + 
    theme_classic() +
    geom_vline(aes(xintercept = 11300)) +
    geom_vline(aes(xintercept = 13800)) + 
    geom_hline(aes(yintercept = 17000)) +
    geom_hline(aes(yintercept = 20000)) +
    theme(axis.text = element_blank())
  dev.off()
  E95.bin20.celltype.immune.eSF.iDSC.df.annotation <- E95.bin20.celltype.immune.eSF.iDSC.df[E95.bin20.celltype.immune.eSF.iDSC.df$x >= 11300 & E95.bin20.celltype.immune.eSF.iDSC.df$x <= 13800 & 
                                                                                              E95.bin20.celltype.immune.eSF.iDSC.df$y >= 17000 & E95.bin20.celltype.immune.eSF.iDSC.df$y <= 20000, ]
  E95.bin20.celltype.df.annotation <- E95.bin20.celltype.df[E95.bin20.celltype.df$x >= 11300 & E95.bin20.celltype.df$x <= 13800 & 
                                                              E95.bin20.celltype.df$y >= 17000 & E95.bin20.celltype.df$y <= 20000, ]
  E95.ann.color <- E95.color[match(intersect(levels(E95.bin20.celltype.immune.eSF.iDSC.df$cell_type), unique(E95.bin20.celltype.immune.eSF.iDSC.df.annotation$cell_type)), 
                                   levels(E95.bin20.celltype.immune.eSF.iDSC.df$cell_type))]
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E95-2-bin20-immune-eSF-iDSC-annotation-C.pdf", width = 5.2, height = 6.2)
  ggplot() + 
    geom_point(E95.bin20.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.eSF.iDSC.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.ann.color) + 
    theme_classic() + 
    theme(legend.position = 'none', axis.text = element_blank()) 
  dev.off()
}

#sptial distribution of iDSC subclusters and eSF in E9.5-normal and CBA samples
iDSC.eSF.color <- c(decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
#E3
DBA.E3.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin50.celltype.sub <- DBA.E3.bin50.celltype[DBA.E3.bin50.celltype$cell_type != "iDSC-DBA1", ]
DBA.E3.bin50.celltype.sub <- DBA.E3.bin50.celltype.sub[, c(1,2,3,4)]
load("20220622-Comments/05.DBA/DBA-E3-iDSC-subclusters-spatial.RData")
DBA.E3.bin50.iDSC.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-iDSC-celltype-20230306.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin50.iDSC.celltype <- DBA.E3.bin50.iDSC.celltype[, c(1,2,3,4)]
DBA.E3.bin50.iDSC.subcluster.df <- rbind(DBA.E3.bin50.celltype.sub, DBA.E3.bin50.iDSC.celltype)
DBA.E3.bin50.iDSC.eSF <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"), ]
DBA.E3.bin50.iDSC.eSF$cell_type <- factor(DBA.E3.bin50.iDSC.eSF$cell_type, levels = c("D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/03.DBA-E3-bin50-iDSC-subclusters.pdf", width = 5, height = 5)
ggplot() + 
  geom_point(DBA.E3.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.E3.bin50.iDSC.eSF, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = iDSC.eSF.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#E6
DBA.E6.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E6.bin50.celltype.sub <- DBA.E6.bin50.celltype[DBA.E6.bin50.celltype$cell_type != "iDSC-DBA1", ]
DBA.E6.bin50.celltype.sub <- DBA.E6.bin50.celltype.sub[, c(1,2,3,4)]
load("20220622-Comments/05.DBA/DBA-E6-iDSC-subclusters-spatial.RData")
DBA.E6.bin50.iDSC.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-iDSC-celltype-20230306.csv", sep = ",", check.names = F, header = T)
DBA.E6.bin50.iDSC.celltype <- DBA.E6.bin50.iDSC.celltype[, c(1,2,3,4)]
DBA.E6.bin50.iDSC.subcluster.df <- rbind(DBA.E6.bin50.celltype.sub, DBA.E6.bin50.iDSC.celltype)
DBA.E6.bin50.iDSC.eSF <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"), ]
DBA.E6.bin50.iDSC.eSF$cell_type <- factor(DBA.E6.bin50.iDSC.eSF$cell_type, levels = c("D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/03.DBA-E6-bin50-iDSC-subclusters.pdf", width = 8, height = 6)
ggplot() + 
  geom_point(DBA.E6.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.E6.bin50.iDSC.eSF, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = iDSC.eSF.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#F1
DBA.F1.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.F1.bin50.celltype.sub <- DBA.F1.bin50.celltype[DBA.F1.bin50.celltype$cell_type != "iDSC-DBA1", ]
DBA.F1.bin50.celltype.sub <- DBA.F1.bin50.celltype.sub[, c(1,2,3,4)]
load("20220622-Comments/05.DBA/DBA-F1-iDSC-subclusters-spatial.RData")
DBA.F1.bin50.iDSC.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin50-iDSC-celltype-20230306.csv", sep = ",", check.names = F, header = T)
DBA.F1.bin50.iDSC.celltype <- DBA.F1.bin50.iDSC.celltype[, c(1,2,3,4)]
DBA.F1.bin50.iDSC.subcluster.df <- rbind(DBA.F1.bin50.celltype.sub, DBA.F1.bin50.iDSC.celltype)
DBA.F1.bin50.iDSC.eSF <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"), ]
DBA.F1.bin50.iDSC.eSF$cell_type <- factor(DBA.F1.bin50.iDSC.eSF$cell_type, levels = c("D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/03.DBA-F1-bin50-iDSC-subclusters.pdf", width = 5.4, height = 5)
ggplot() + 
  geom_point(DBA.F1.bin50.iDSC.subcluster.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size= 0.5) + 
  geom_point(DBA.F1.bin50.iDSC.eSF, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = iDSC.eSF.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#Normal E95
E95.bin50.celltype.df <- read.table("20210827-figures/figure5/E95-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
E95.bin50.iDSC.eSF <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% c("Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"), ]
E95.bin50.iDSC.eSF$cell_type <- factor(E95.bin50.iDSC.eSF$cell_type, levels = c("Lum-D", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/03.E95-bin50-iDSC-subclusters.pdf", width = 10, height = 8)
ggplot() + 
  geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size= 0.5) + 
  geom_point(E95.bin50.iDSC.eSF, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = iDSC.eSF.color) + 
  theme_classic() + 
  theme(axis.text = element_blank())
dev.off()

#CBA immune cell distribution - bin20
CBA.immune.eSF.color <- c(immune.color[1:4], immune.color[6:7], decidual.color[1], "#F8766D", "#00BA38", "#619CFF")
#E3
DBA.E3.bin20.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin20-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin20.celltype <- DBA.E3.bin20.celltype[, 1:4]
DBA.E3.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin50.celltype <- DBA.E3.bin50.celltype[, 1:4]
DBA.E3.bin20.celltype.immune.df <- DBA.E3.bin20.celltype[DBA.E3.bin20.celltype$cell_type %in% c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK"), ]
DBA.E3.bin20.celltype.immune.eSF.iDSC.df <- rbind(DBA.E3.bin20.celltype.immune.df, DBA.E3.bin50.iDSC.eSF)
DBA.E3.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(DBA.E3.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                          levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))
DBA.E3.bin20.celltype.immune.eSF.df <- DBA.E3.bin20.celltype.immune.eSF.iDSC.df[DBA.E3.bin20.celltype.immune.eSF.iDSC.df$cell_type %in% c("D1-eSF"), ]
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/04.DBA-E3-bin20-immune-eSF-iDSC.pdf", width = 11, height = 11)
ggplot() + 
  geom_point(DBA.E3.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(DBA.E3.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  geom_point(DBA.E3.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 0.5) + 
  scale_color_manual(values = CBA.immune.eSF.color) +
  theme_classic() +
  theme(axis.text = element_blank())
dev.off()

#E6
DBA.E6.bin20.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin20-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E6.bin20.celltype <- DBA.E6.bin20.celltype[, 1:4]
DBA.E6.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E6.bin50.celltype <- DBA.E6.bin50.celltype[, 1:4]
DBA.E6.bin20.celltype.immune.df <- DBA.E6.bin20.celltype[DBA.E6.bin20.celltype$cell_type %in% c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK"), ]
DBA.E6.bin20.celltype.immune.eSF.iDSC.df <- rbind(DBA.E6.bin20.celltype.immune.df, DBA.E6.bin50.iDSC.eSF)
DBA.E6.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(DBA.E6.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                             levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))
DBA.E6.bin20.celltype.immune.eSF.df <- DBA.E6.bin20.celltype.immune.eSF.iDSC.df[DBA.E6.bin20.celltype.immune.eSF.iDSC.df$cell_type %in% c("D1-eSF"), ]
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/04.DBA-E6-bin20-immune-eSF-iDSC.pdf", width = 18, height = 15)
ggplot() + 
  geom_point(DBA.E6.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(DBA.E6.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  geom_point(DBA.E6.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 0.5) + 
  scale_color_manual(values = CBA.immune.eSF.color) +
  theme_classic() +
  theme(axis.text = element_blank())
dev.off()

#F1
DBA.F1.bin20.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin20-celltype.csv", sep = ",", check.names = F, header = T)
DBA.F1.bin20.celltype <- DBA.F1.bin20.celltype[, 1:4]
DBA.F1.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.F1.bin50.celltype <- DBA.F1.bin50.celltype[, 1:4]
DBA.F1.bin20.celltype.immune.df <- DBA.F1.bin20.celltype[DBA.F1.bin20.celltype$cell_type %in% c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK"), ]
DBA.F1.bin20.celltype.immune.eSF.iDSC.df <- rbind(DBA.F1.bin20.celltype.immune.df, DBA.F1.bin50.iDSC.eSF)
DBA.F1.bin20.celltype.immune.eSF.iDSC.df$cell_type <- factor(DBA.F1.bin20.celltype.immune.eSF.iDSC.df$cell_type, 
                                                             levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))
DBA.F1.bin20.celltype.immune.eSF.df <- DBA.F1.bin20.celltype.immune.eSF.iDSC.df[DBA.F1.bin20.celltype.immune.eSF.iDSC.df$cell_type %in% c("D1-eSF"), ]
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/04.DBA-F1-bin20-immune-eSF-iDSC.pdf", width = 14, height = 14)
ggplot() + 
  geom_point(DBA.F1.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(DBA.F1.bin20.celltype.immune.eSF.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  geom_point(DBA.F1.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y), color = decidual.color[1], size = 0.5) + 
  scale_color_manual(values = CBA.immune.eSF.color) +
  theme_classic() +
  theme(axis.text = element_blank())
dev.off()

#immune cell + eSF + iDSCs, bin50 distribution
#DBA.E3.bin50
DBA.E3.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin50.celltype <- DBA.E3.bin50.celltype[, 1:4]
DBA.E3.bin50.iDSC.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-iDSC-celltype-20230306.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin50.iDSC.celltype <- DBA.E3.bin50.iDSC.celltype[, 1:4]
DBA.E3.bin50.celltype.immune.df <- DBA.E3.bin50.celltype[DBA.E3.bin50.celltype$cell_type %in% 
                                                        c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "D1-eSF"), ]
DBA.E3.bin50.celltype.immune.df <- rbind(DBA.E3.bin50.celltype.immune.df, DBA.E3.bin50.iDSC.celltype)
DBA.E3.bin50.celltype.immune.df$cell_type <- factor(DBA.E3.bin50.celltype.immune.df$cell_type, 
                                                 levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "D1-eSF", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))

pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.DBA.E3-2-bin50-immune-rect.pdf", width = 10, height = 8)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(DBA.E3.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = CBA.immune.eSF.color) + 
  theme_classic() +
  geom_vline(aes(xintercept = 11300)) +
  geom_vline(aes(xintercept = 13800)) + 
  geom_hline(aes(yintercept = 17000)) +
  geom_hline(aes(yintercept = 20000)) +
  theme(axis.text = element_blank())
dev.off()
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.DBA.E3-2-bin50-immune.pdf", width = 10, height = 8)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(DBA.E3.bin50.celltype.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = CBA.immune.eSF.color) + 
  theme_classic() +
  theme(axis.text = element_blank())
dev.off()
DBA.E3.bin50.celltype.immune.df.annotation <- DBA.E3.bin50.celltype.immune.df[DBA.E3.bin50.celltype.immune.df$x >= 11300 & DBA.E3.bin50.celltype.immune.df$x <= 13800 & 
                                                                          DBA.E3.bin50.celltype.immune.df$y >= 17000 & DBA.E3.bin50.celltype.immune.df$y <= 20000, ]
DBA.E3.bin50.celltype.df.annotation <- DBA.E3.bin50.celltype.df[DBA.E3.bin50.celltype.df$x >= 11300 & DBA.E3.bin50.celltype.df$x <= 13800 & 
                                                            DBA.E3.bin50.celltype.df$y >= 17000 & DBA.E3.bin50.celltype.df$y <= 20000, ]
DBA.E3.ann.color <- DBA.E3.color[match(intersect(levels(DBA.E3.bin50.celltype.immune.df$cell_type), unique(DBA.E3.bin50.celltype.immune.df.annotation$cell_type)), 
                                 levels(DBA.E3.bin50.celltype.immune.df$cell_type))]
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/02.DBA.E3-2-bin50-immune-annotation.pdf", width = 2.3, height = 2.7)
ggplot() + 
  geom_point(DBA.E3.bin50.celltype.df.annotation, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
  geom_point(DBA.E3.bin50.celltype.immune.df.annotation, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
  scale_color_manual(values = DBA.E3.ann.color) + 
  theme_classic() + 
  theme(legend.position = 'none', axis.text = element_blank()) 
dev.off()

############################################bin20 immune - 20230501################################################
#Normal tissue Immune subclusters + iDSC subclusters + eSF and zoom in pics
function(){
  #E65.2.bin20
  E65.2.bin20.celltype.df <- read.table("20210827-figures/figure5/E65-2-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E65.2.bin50.celltype.df <- read.table("20210827-figures/figure5/E65-2-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E65.2.bin20.celltype.immune.df <- E65.2.bin20.celltype.df[E65.2.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E65.2.bin50.eSF.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% c("Lum-D"), ]
  E65.2.bin50.iDSC.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% c("FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E65.2.bin20.celltype.immune.eSF.df <- rbind(E65.2.bin20.celltype.immune.df, E65.2.bin50.eSF.df)
  E65.2.bin20.celltype.immune.eSF.df$cell_type <- factor(E65.2.bin20.celltype.immune.eSF.df$cell_type, 
                                                              levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"))
  E65.2.color <- c(immune.color[1:7], decidual.color[1])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E65-2-bin20-immune-eSF.pdf", width =7.5, height = 6)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  E65.2.bin50.iDSC.df <- E65.2.bin50.celltype.df[E65.2.bin50.celltype.df$cell_type %in% c("FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E65.2.bin20.celltype.immune.iDSC.df <- rbind(E65.2.bin20.celltype.immune.df, E65.2.bin50.iDSC.df)
  E65.2.bin20.celltype.immune.iDSC.df$cell_type <- factor(E65.2.bin20.celltype.immune.iDSC.df$cell_type, 
                                                         levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E65.2.color <- c(immune.color[1:7], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E65-2-bin20-immune-iDSC.pdf", width =7.5, height = 6)
  ggplot() + 
    geom_point(E65.2.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E65.2.bin20.celltype.immune.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E65.2.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  E75.1.bin20.celltype.df <- read.table("20210827-figures/figure5/E75-1-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E75.1.bin50.celltype.df <- read.table("20210827-figures/figure5/E75-1-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E75.1.bin20.celltype.immune.df <- E75.1.bin20.celltype.df[E75.1.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E75.1.bin50.eSF.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% c("Lum-D"), ]
  E75.1.bin50.iDSC.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% c("FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E75.1.bin20.celltype.immune.eSF.df <- rbind(E75.1.bin20.celltype.immune.df, E75.1.bin50.eSF.df)
  E75.1.bin20.celltype.immune.eSF.df$cell_type <- factor(E75.1.bin20.celltype.immune.eSF.df$cell_type, 
                                                         levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"))
  E75.1.color <- c(immune.color[1:7], decidual.color[1])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E75-1-bin20-immune-eSF.pdf", width = 7.5, height = 6)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  E75.1.bin50.iDSC.df <- E75.1.bin50.celltype.df[E75.1.bin50.celltype.df$cell_type %in% c("FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E75.1.bin20.celltype.immune.iDSC.df <- rbind(E75.1.bin20.celltype.immune.df, E75.1.bin50.iDSC.df)
  E75.1.bin20.celltype.immune.iDSC.df$cell_type <- factor(E75.1.bin20.celltype.immune.iDSC.df$cell_type, 
                                                          levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E75.1.color <- c(immune.color[1:7], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E75-1-bin20-immune-iDSC.pdf", width = 7.5, height = 6)
  ggplot() + 
    geom_point(E75.1.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E75.1.bin20.celltype.immune.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E75.1.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()

  #E85.B6.bin20
  E85.B6.bin20.celltype.df <- read.table("20210827-figures/figure5/E85-B6-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E85.B6.bin50.celltype.df <- read.table("20210827-figures/figure5/E85-B6-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E85.B6.bin20.celltype.immune.df <- E85.B6.bin20.celltype.df[E85.B6.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E85.B6.bin50.eSF.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% c("Lum-D"), ]
  E85.B6.bin50.iDSC.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% c("FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E85.B6.bin20.celltype.immune.eSF.df <- rbind(E85.B6.bin20.celltype.immune.df, E85.B6.bin50.eSF.df)
  E85.B6.bin20.celltype.immune.eSF.df$cell_type <- factor(E85.B6.bin20.celltype.immune.eSF.df$cell_type, 
                                                         levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"))
  E85.B6.color <- c(immune.color[1:7], decidual.color[1])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E85-B6-bin20-immune-eSF.pdf", width = 10.5, height = 10)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  E85.B6.bin50.iDSC.df <- E85.B6.bin50.celltype.df[E85.B6.bin50.celltype.df$cell_type %in% c("FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E85.B6.bin20.celltype.immune.iDSC.df <- rbind(E85.B6.bin20.celltype.immune.df, E85.B6.bin50.iDSC.df)
  E85.B6.bin20.celltype.immune.iDSC.df$cell_type <- factor(E85.B6.bin20.celltype.immune.iDSC.df$cell_type, 
                                                          levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E85.B6.color <- c(immune.color[1:7], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E85-B6-bin20-immune-iDSC.pdf", width = 10.5, height = 10)
  ggplot() + 
    geom_point(E85.B6.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E85.B6.bin20.celltype.immune.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E85.B6.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  #E95.bin20
  E95.bin20.celltype.df <- read.table("20210827-figures/figure5/E95-bin20-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E95.bin50.celltype.df <- read.table("20210827-figures/figure5/E95-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
  E95.bin20.celltype.immune.df <- E95.bin20.celltype.df[E95.bin20.celltype.df$cell_type %in% c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK"), ]
  E95.bin50.eSF.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% c("Lum-D"), ]
  E95.bin50.iDSC.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% c("FB_immune_0", "FB_immune_1", "FB_immune_2"), ]
  E95.bin20.celltype.immune.eSF.df <- rbind(E95.bin20.celltype.immune.df, E95.bin50.eSF.df)
  E95.bin20.celltype.immune.eSF.df$cell_type <- factor(E95.bin20.celltype.immune.eSF.df$cell_type, 
                                                          levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "Lum-D"))
  E95.color <- c(immune.color[1:7], decidual.color[1])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E95-bin20-immune-eSF.pdf", width = 15, height = 13)
  ggplot() + 
    geom_point(E95.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  E95.bin50.iDSC.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% c("FB_immune_1", "FB_immune_0", "FB_immune_2"), ]
  E95.bin20.celltype.immune.iDSC.df <- rbind(E95.bin20.celltype.immune.df, E95.bin50.iDSC.df)
  E95.bin20.celltype.immune.iDSC.df$cell_type <- factor(E95.bin20.celltype.immune.iDSC.df$cell_type, 
                                                           levels = c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK", "FB_immune_0", "FB_immune_1", "FB_immune_2"))
  E95.color <- c(immune.color[1:7], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/01.E95-bin20-immune-iDSC.pdf", width = 15, height = 13)
  ggplot() + 
    geom_point(E95.bin20.celltype.df, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(E95.bin20.celltype.immune.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5) +
    scale_color_manual(values = E95.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
}

#CBA tissue Immune subclusters + iDSC subclusters + eSF and zoom in pics
function(){
  DBA.E3.bin20.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin20-celltype.csv", sep = ",", check.names = F, header = T)
  DBA.E3.bin20.celltype <- DBA.E3.bin20.celltype[, 1:4]
  DBA.E3.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-celltype.csv", sep = ",", check.names = F, header = T)
  DBA.E3.bin50.celltype <- DBA.E3.bin50.celltype[, 1:4]
  DBA.E3.bin20.celltype.immune.df <- DBA.E3.bin20.celltype[DBA.E3.bin20.celltype$cell_type %in% c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK"), ]
  DBA.E3.bin50.eSF <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
  DBA.E3.bin50.iDSC <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"), ]
  DBA.E3.bin20.celltype.immune.eSF.df <- rbind(DBA.E3.bin20.celltype.immune.df, DBA.E3.bin50.eSF)
  DBA.E3.bin20.celltype.immune.eSF.df$cell_type <- factor(DBA.E3.bin20.celltype.immune.eSF.df$cell_type, 
                                                          levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "D1-eSF"))
  E3.color <- c(immune.color[1:4], immune.color[6:7], decidual.color[1])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/04.DBA-E3-bin20-immune-eSF.pdf", width = 11, height = 11)
  ggplot() + 
    geom_point(DBA.E3.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E3.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5)+ 
    scale_color_manual(values = E3.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  DBA.E3.bin50.iDSC <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"), ]
  DBA.E3.bin20.celltype.immune.iDSC.df <- rbind(DBA.E3.bin20.celltype.immune.df, DBA.E3.bin50.iDSC)
  DBA.E3.bin20.celltype.immune.iDSC.df$cell_type <- factor(DBA.E3.bin20.celltype.immune.iDSC.df$cell_type, 
                                                           levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))
  E3.color <- c(immune.color[1:4], immune.color[6:7], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/04.DBA-E3-bin20-immune-iDSC.pdf", width = 11, height = 11)
  ggplot() + 
    geom_point(DBA.E3.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E3.bin20.celltype.immune.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5)+ 
    scale_color_manual(values = E3.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  DBA.E6.bin20.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin20-celltype.csv", sep = ",", check.names = F, header = T)
  DBA.E6.bin20.celltype <- DBA.E6.bin20.celltype[, 1:4]
  DBA.E6.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-celltype.csv", sep = ",", check.names = F, header = T)
  DBA.E6.bin50.celltype <- DBA.E6.bin50.celltype[, 1:4]
  DBA.E6.bin20.celltype.immune.df <- DBA.E6.bin20.celltype[DBA.E6.bin20.celltype$cell_type %in% c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK"), ]
  DBA.E6.bin50.eSF <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
  DBA.E6.bin50.iDSC <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"), ]
  DBA.E6.bin20.celltype.immune.eSF.df <- rbind(DBA.E6.bin20.celltype.immune.df, DBA.E6.bin50.eSF)
  DBA.E6.bin20.celltype.immune.eSF.df$cell_type <- factor(DBA.E6.bin20.celltype.immune.eSF.df$cell_type, 
                                                          levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "D1-eSF"))
  E6.color <- c(immune.color[1:4], immune.color[6:7], decidual.color[1])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/04.DBA-E6-bin20-immune-eSF.pdf", width = 18, height = 15)
  ggplot() + 
    geom_point(DBA.E6.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E6.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5)+ 
    scale_color_manual(values = E6.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  DBA.E6.bin50.iDSC <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"), ]
  DBA.E6.bin20.celltype.immune.iDSC.df <- rbind(DBA.E6.bin20.celltype.immune.df, DBA.E6.bin50.iDSC)
  DBA.E6.bin20.celltype.immune.iDSC.df$cell_type <- factor(DBA.E6.bin20.celltype.immune.iDSC.df$cell_type, 
                                                           levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))
  E6.color <- c(immune.color[1:4], immune.color[6:7], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/04.DBA-E6-bin20-immune-iDSC.pdf", width = 18, height = 15)
  ggplot() + 
    geom_point(DBA.E6.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.E6.bin20.celltype.immune.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5)+ 
    scale_color_manual(values = E6.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  DBA.F1.bin20.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin20-celltype.csv", sep = ",", check.names = F, header = T)
  DBA.F1.bin20.celltype <- DBA.F1.bin20.celltype[, 1:4]
  DBA.F1.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-F1-bin50-celltype.csv", sep = ",", check.names = F, header = T)
  DBA.F1.bin50.celltype <- DBA.F1.bin50.celltype[, 1:4]
  
  DBA.F1.bin20.celltype.immune.df <- DBA.F1.bin20.celltype[DBA.F1.bin20.celltype$cell_type %in% c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK"), ]
  DBA.F1.bin50.eSF <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
  DBA.F1.bin50.iDSC <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"), ]
  DBA.F1.bin20.celltype.immune.eSF.df <- rbind(DBA.F1.bin20.celltype.immune.df, DBA.F1.bin50.eSF)
  DBA.F1.bin20.celltype.immune.eSF.df$cell_type <- factor(DBA.F1.bin20.celltype.immune.eSF.df$cell_type, 
                                                          levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "D1-eSF"))
  F1.color <- c(immune.color[1:4], immune.color[6:7], decidual.color[1])
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/04.DBA-F1-bin20-immune-eSF.pdf", width = 14, height = 14)
  ggplot() + 
    geom_point(DBA.F1.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.F1.bin20.celltype.immune.eSF.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5)+ 
    scale_color_manual(values = F1.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
  
  DBA.F1.bin50.iDSC <- DBA.F1.bin50.iDSC.subcluster.df[DBA.F1.bin50.iDSC.subcluster.df$cell_type %in% c("DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"), ]
  DBA.F1.bin20.celltype.immune.iDSC.df <- rbind(DBA.F1.bin20.celltype.immune.df, DBA.F1.bin50.iDSC)
  DBA.F1.bin20.celltype.immune.iDSC.df$cell_type <- factor(DBA.F1.bin20.celltype.immune.iDSC.df$cell_type, 
                                                           levels = c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK", "DBA.1.iDSC0", "DBA.1.iDSC1", "DBA.1.iDSC2"))
  F1.color <- c(immune.color[1:4], immune.color[6:7], "#F8766D", "#00BA38", "#619CFF")
  pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/04.DBA-F1-bin20-immune-iDSC.pdf", width = 14, height = 14)
  ggplot() + 
    geom_point(DBA.F1.bin20.celltype, mapping = aes(x = x, y = y), color = "#E7E7E7", size = 0.5) + 
    geom_point(DBA.F1.bin20.celltype.immune.iDSC.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.5)+ 
    scale_color_manual(values = F1.color) +
    theme_classic() +
    theme(axis.text = element_blank())
  dev.off()
}

#20230507 immune cell number in EC hub
#Normal - E9.5, 14000<x<17500, 16500<y<20000
E95.bin50.celltype.df <- read.table("20210827-figures/figure5/E95-bin50-FB-immune-subclusters-mapping-all-slide.csv", sep = ",", check.names = F)
E95.bin50.eSF.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% c("Lum-D"), ]
E95.bin50.eSF.sub.df <- E95.bin50.eSF.df[E95.bin50.eSF.df$x<15000 & E95.bin50.eSF.df$x>12500 & E95.bin50.eSF.df$y<20000 & E95.bin50.eSF.df$y>17500, ]
E95.bin50.tropho.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% c("tropho-0", "tropho-1", "tropho-2", "tropho-3", "tropho-4"), ]
E95.bin50.immune.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$cell_type %in% E95.immune, ]

ggplot() + 
  geom_point(E95.bin50.celltype.df, mapping = aes(x = x, y = y), color = "lightgray", size = 0.1) + 
  geom_point(E95.bin50.tropho.df, mapping = aes(x = x, y = y), size = 0.1) + 
  geom_point(E95.bin50.eSF.df, mapping = aes(x = x, y = y), size = 0.1) + 
  geom_vline(aes(xintercept = 14000)) +
  geom_vline(aes(xintercept = 17500)) +
  geom_hline(aes(yintercept = 20000)) +
  geom_hline(aes(yintercept = 16500)) +
  geom_point(E95.bin50.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.1) +
  scale_color_manual(values = E95.color)
  
E95.immune <- c("Mono","Mac","Pro-Mac","DC-1","DC-2","Neutrophil","NK")
E95.color <- c(immune.color[1:7])
E95.bin50.sub.df <- E95.bin50.celltype.df[E95.bin50.celltype.df$x > 14000 & E95.bin50.celltype.df$x < 17500 &
                                            E95.bin50.celltype.df$y > 16500 & E95.bin50.celltype.df$y < 20000, ]
E95.bin50.immune.sub.df <- E95.bin50.sub.df[E95.bin50.sub.df$cell_type %in% E95.immune, ]
E95.bin50.immune.sub.ration.df <- data.frame(table(E95.bin50.immune.sub.df$cell_type)/4761, check.names = F)
E95.bin50.immune.sub.ration.df$Var1 <- c("Monocyte", "NK")
E95.bin50.immune.sub.ration.df$sample <- "Normal_E95"

#CBA/J - E3, 12750<y<15000
DBA.E3.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E3-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E3.bin50.celltype <- DBA.E3.bin50.celltype[, 1:4]
DBA.E3.bin50.eSF <- DBA.E3.bin50.iDSC.subcluster.df[DBA.E3.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
DBA.E3.bin50.tropho.df <- DBA.E3.bin50.celltype[DBA.E3.bin50.celltype$cell_type %in% c("tropho-0", "tropho-1", "tropho-2", "tropho-3", "tropho-4"), ]
E3.bin50.immune.df <- DBA.E3.bin50.celltype[DBA.E3.bin50.celltype$cell_type %in% E3.immune, ]

ggplot() + 
  geom_point(DBA.E3.bin50.celltype, mapping = aes(x = x, y = y), color = "lightgray") + 
  geom_point(DBA.E3.bin50.tropho.df, mapping = aes(x = x, y = y)) + 
  geom_point(DBA.E3.bin50.eSF, mapping = aes(x = x, y = y)) + 
  geom_hline(aes(yintercept = 15000)) +
  geom_hline(aes(yintercept = 12750)) + 
  geom_vline(aes(xintercept = 8500)) +
  geom_vline(aes(xintercept = 13500)) +
  geom_point(E3.bin50.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.1) +
  scale_color_manual(values = E3.color)

E3.immune <- c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK")
E3.color <- c(immune.color[1:4], immune.color[6:7])
E3.bin50.sub.df <- DBA.E3.bin50.celltype[DBA.E3.bin50.celltype$x > 8500 & DBA.E3.bin50.celltype$x < 13500 &
                                                 DBA.E3.bin50.celltype$y > 12750 & DBA.E3.bin50.celltype$y < 15000, ]
E3.bin50.immune.sub.df <- E3.bin50.sub.df[E3.bin50.sub.df$cell_type %in% E3.immune, ]
E3.bin50.immune.sub.ratio.df <- data.frame(table(E3.bin50.immune.sub.df$cell_type)/4352, check.names = F)
E3.bin50.immune.sub.ratio.df$sample <- "CBA_E95-a"
  
#CBA/J - E6, 13400<y<16400
DBA.E6.bin50.celltype <- read.table("20220622-Comments/05.DBA/DBA-E6-bin50-celltype.csv", sep = ",", check.names = F, header = T)
DBA.E6.bin50.celltype <- DBA.E6.bin50.celltype[, 1:4]
DBA.E6.bin50.eSF <- DBA.E6.bin50.iDSC.subcluster.df[DBA.E6.bin50.iDSC.subcluster.df$cell_type %in% c("D1-eSF"), ]
DBA.E6.bin50.tropho.df <- DBA.E6.bin50.celltype[DBA.E6.bin50.celltype$cell_type %in% c("tropho-0", "tropho-1", "tropho-2", "tropho-3", "tropho-4"), ]
E6.bin50.immune.df <- DBA.E6.bin50.celltype[DBA.E6.bin50.celltype$cell_type %in% E6.immune, ]

ggplot() + 
  geom_point(DBA.E6.bin50.celltype, mapping = aes(x = x, y = y), color = "lightgray", size = 0.1) + 
  geom_point(DBA.E6.bin50.tropho.df, mapping = aes(x = x, y = y), size = 0.1) + 
  geom_point(DBA.E6.bin50.eSF, mapping = aes(x = x, y = y), size = 0.1) + 
  geom_vline(aes(xintercept = 13750)) + 
  geom_vline(aes(xintercept = 16450)) + 
  geom_hline(aes(yintercept = 15000)) + 
  geom_hline(aes(yintercept = 10000)) +
  geom_point(E6.bin50.immune.df, mapping = aes(x = x, y = y, color = cell_type), size = 0.1) +
  scale_color_manual(values = E6.color)
  
E6.immune <- c("Monocyte","Mac","Proliferating Mac","DC-1","Neutrophil","NK")
E6.color <- c(immune.color[1:4], immune.color[6:7])
E6.bin50.sub.df <- DBA.E6.bin50.celltype[DBA.E6.bin50.celltype$x > 13750 & DBA.E6.bin50.celltype$x < 16450 &
                                           DBA.E6.bin50.celltype$y > 10000 & DBA.E6.bin50.celltype$y < 15000, ]
E6.bin50.immune.sub.df <- E6.bin50.sub.df[E6.bin50.sub.df$cell_type %in% E6.immune, ]
E6.bin50.immune.sub.ratio.df <- data.frame(table(E6.bin50.immune.sub.df$cell_type)/5247, check.names = F)
E6.bin50.immune.sub.ratio.df$sample <- "CBA_E95-b"

Vascular.hub.immune.ratio <- rbind(E95.bin50.immune.sub.ration.df, E3.bin50.immune.sub.ratio.df, E6.bin50.immune.sub.ratio.df)
Vascular.hub.immune.ratio$sample <- factor(Vascular.hub.immune.ratio$sample, 
                                           levels = c("Normal_E95", "CBA_E95-a", "CBA_E95-b"))
Vascular.hub.immune.ratio$Var1 <- factor(Vascular.hub.immune.ratio$Var1, levels = E6.immune)
pdf("20220622-Comments/05.DBA/00.figures/20230308-hubs/14.Vascular-hub-immune-cell-ratio-bin50.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(Vascular.hub.immune.ratio, mapping = aes(x= sample, y = Freq, fill = Var1), stat = "identity") + 
  scale_fill_manual(values = immune.color[c(1,2,4,6,7)]) + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", color = "black"))
dev.off()

#(Pro) macrophage - VEGF signal pathways








