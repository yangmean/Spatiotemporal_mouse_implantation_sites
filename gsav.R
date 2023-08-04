library("msigdbr")
library("dplyr")
library("GSVA")
#查看支持的物种
msigdbr_show_species()
#mouse gene set
m_df = msigdbr(species = "Mus musculus")
#查看基因集
a <- m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
#挑选基因集
m_df = msigdbr(species = "Mus musculus", category = "H")
#构建参考数据集
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
#提取数据
Mac.data <- as.matrix(Mac@assays$RNA@data)
meta <- Mac@meta.data[,c("tdTomato")]
meta <- gsub("TRUE", "positive", meta)
meta <- gsub("FALSE", "negative", meta)
#运行GSVA，使用的count数据
mac.cp <- gsva(Mac.data, msigdbr_list, kcdf = "Gaussian", method = "gsva", parallel.sz = 20)
#limma做通路的差异分析
library(limma)
group <- factor(meta,levels = c("positive", "negative"),ordered = F)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(Mac.data)
fit <- lmFit(mac.cp, design)
con.matrx <- makeContrasts(contrasts = c("positive-negative"), levels = design)
fit2 <- contrasts.fit(fit, con.matrx)
fit2 <- eBayes(fit2)
mac.pathway <- topTable(fit2, coef = "positive-negative", number = 20)
mac.pathway$pathway <- rownames(mac.pathway)
write.table(mac.pathway, file = "20210827-figures/figure4/macrophage-postive-negtive.csv", sep = ",", quote = F)
mac.pathway$pathway <- factor(mac.pathway$pathway, levels = mac.pathway[order(mac.pathway$t), ]$pathway)
pdf("20210827-figures/figure4/Figure4-2-tdTomato-compare-GSVA.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(mac.pathway, mapping = aes(x = t, y = pathway), stat = "identity") +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "black", fill = "white"))
dev.off()

#做FB-immune的两个位置的GSVA分析
FB.immune.data <- as.matrix(B6.bin50.FB.immune@assays$SCT@data)
meta <- B6.bin50.FB.immune@meta.data[,c("cluster_celltype")]
meta <- gsub("3_FB-immune", "FB3", meta)
meta <- gsub("6_FB-immune", "FB6", meta)
#运行GSVA，使用的count数据
FB.immune.cp <- gsva(FB.immune.data, msigdbr_list, kcdf = "Gaussian", method = "gsva", parallel.sz = 20)
#limma做通路的差异分析
library(limma)
group <- factor(meta,levels = c("FB3", "FB6"),ordered = F)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(FB.immune.data)
fit <- lmFit(FB.immune.cp, design)
#构建比较矩阵
con.matrx <- makeContrasts(contrasts = c("FB3-FB6"), levels = design)
fit2 <- contrasts.fit(fit, con.matrx)
fit2 <- eBayes(fit2)
FB.immune.pathway  <- topTable(fit2, coef = "FB3-FB6", number = 15, lfc = 0)#number是拿出最大的几个结果
FB.immune.pathway$pathway <- rownames(FB.immune.pathway)
FB.immune.pathway$pathway <- factor(FB.immune.pathway$pathway, levels = FB.immune.pathway[order(FB.immune.pathway$t), ]$pathway)
write.table(FB.immune.BP, file = "20210827-figures/figure4/FB.immune-cluster3-6-H.csv", sep = ",", quote = F)
pdf("20210827-figures/figure4/Figure4-2-FB-immune-cluster3-6-GSVA.pdf", width = 10, height = 10)
ggplot() + 
  geom_bar(FB.immune.pathway, mapping = aes(x = t, y = pathway), stat = "identity") + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "black", fill = "white"))
dev.off()

#cluster1,3,6的cycling DSC的comparing
cycling.DSC.data <- as.matrix(B6.bin50.cycling.DSC@assays$SCT@data)
meta <- B6.bin50.cycling.DSC@meta.data[,c("cluster_celltype")]
meta <- gsub("1_Top2a-D", "cluster1", meta)
meta <- gsub("3_Top2a-D", "cluster3", meta)
meta <- gsub("6_Top2a-D", "cluster6", meta)
#运行GSVA，使用的count数据
cycling.DSC.cp <- gsva(cycling.DSC.data, msigdbr_list, kcdf = "Gaussian", method = "gsva", parallel.sz = 20)
#limma做通路的差异分析
library(limma)
group <- factor(meta,levels = c("cluster1", "cluster3", "cluster6"),ordered = F)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(cycling.DSC.data)
fit <- lmFit(cycling.DSC.cp, design)
#构建比较矩阵
con.matrx <- makeContrasts(contrasts = c("cluster1-cluster3", "cluster3-cluster6", "cluster6-cluster1"), levels = design)
fit2 <- contrasts.fit(fit, con.matrx)
fit2 <- eBayes(fit2)
cycling.DSC.pathway <- topTable(fit2, number = 20)
cycling.DSC.cp.pathway <- cycling.DSC.cp[rownames(cycling.DSC.pathway), ]
cycling.DSC.cp.pathway <- data.frame(t(cycling.DSC.cp.pathway), check.names = F)
cycling.DSC.cp.pathway$cellname <- rownames(cycling.DSC.cp.pathway)
cellcluster <- data.frame(cellname = colnames(B6.bin50.cycling.DSC), cluster = meta)
cycling.DSC.cp.pathway <- merge(cycling.DSC.cp.pathway, cellcluster, by = "cellname", all = F)
cycling.DSC.cp.pathway <- aggregate(cycling.DSC.cp.pathway[, 2:21], by = list(cycling.DSC.cp.pathway$cluster), FUN = mean)
rownames(cycling.DSC.cp.pathway) <- cycling.DSC.cp.pathway$Group.1
pdf("20210827-figures/figure4/Figure4-12-Cycling-DSC-cluster1-3-6-GSVA.pdf", width = 15, height = 10)
pheatmap(cycling.DSC.cp.pathway[, -1], scale = "column", cluster_rows = T, cluster_cols = T)
dev.off()

#FB-immune分亚群之后的分析
FB.immune.data <- as.matrix(FB.immune.SCT@assays$SCT@data)
meta <- FB.immune.SCT@meta.data[,c("seurat_clusters")]
meta <- paste("cluster", meta, sep = "")
#运行GSVA，使用的count数据
FB.immune.SCT.cp <- gsva(FB.immune.data, msigdbr_list, kcdf = "Gaussian", method = "gsva", parallel.sz = 20)
#limma做通路的差异分析
library(limma)
group <- factor(meta,levels = c("cluster0", "cluster1", "cluster2"),ordered = F)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(FB.immune.data)
fit.FB.sub <- lmFit(FB.immune.SCT.cp, design)
#构建比较矩阵
con.matrx <- makeContrasts(contrasts = c("cluster0-cluster1", "cluster2-cluster1", "cluster0-cluster2"), levels = design)
fit2.FB.sub <- contrasts.fit(fit.FB.sub, con.matrx)
fit2.FB.sub <- eBayes(fit2.FB.sub)
FB.immune.SCT.cp.pathway.t <- fit2.FB.sub$t
color <- brewer.pal(11, "Spectral")
color <- color[c(1:5, 7:11)]
#cut color, <0 的用蓝色绿色，>0 的用黄色红色
#breaks
bk <- c(seq(-20,0,by=0.1),seq(0,20,by=0.1))
bk <- unique(bk)
pdf("20210827-figures/figure5/Figure5-FB-imune-3clusters-gsva-heatmap.pdf")
pheatmap(FB.immune.SCT.cp.pathway.t,
         scale = "none",cluster_cols = T,
         color = rev(c(colorRampPalette(colors = c(color[1:3], "white"))(length(bk)/2),
                       colorRampPalette(colors = c("white", color[9:11]))(length(bk)/2))),
         legend_breaks=seq(-20,20,5),
         breaks=bk)
dev.off()

#每个cluster单独用GSVA的结果来show
FB.immune.SCT.pathway <- topTable(fit2, number = 50)#取全部结果
FB.immune.SCT.cp.pathway <- FB.immune.SCT.cp[rownames(FB.immune.SCT.pathway), ]
FB.immune.SCT.cp.pathway <- data.frame(t(FB.immune.SCT.cp.pathway), check.names = F)
FB.immune.SCT.cp.pathway$cellname <- rownames(FB.immune.SCT.cp.pathway)
cellcluster <- data.frame(cellname = colnames(FB.immune.SCT), cluster = FB.immune.SCT$seurat_clusters)
FB.immune.SCT.cp.pathway <- merge(FB.immune.SCT.cp.pathway, cellcluster, by = "cellname", all = F)
FB.immune.SCT.cp.pathway <- aggregate(FB.immune.SCT.cp.pathway[, 2:51], by = list(FB.immune.SCT.cp.pathway$cluster), FUN = mean)
rownames(FB.immune.SCT.cp.pathway) <- FB.immune.SCT.cp.pathway$Group.1
pdf("20210827-figures/figure5/Figure5-FB-imune-3clusters-gsva-heatmap-GSVA.pdf", width = 15, height = 10)
pheatmap(FB.immune.SCT.cp.pathway[, -1], scale = "column", cluster_rows = T, cluster_cols = T)
dev.off()







