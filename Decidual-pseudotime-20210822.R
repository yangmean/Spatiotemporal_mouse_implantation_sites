#monocle2
library(monocle)
library(Seurat)
load("../20210201-integrate/UE-decidual-harmony-pc30-rename.RData")
DimPlot(UE.decidual)
UE.decidual.susbet <- subset(UE.decidual, subset = cellname != "cluster9-1" & cellname != "Prg-" )
UE.decidual.susbet <- subset(UE.decidual, subset = cellname != "cluster9-1" & cellname != "Lum-D" & cellname != "Postn-D" & cellname != "Prg-" )
UE.decidual.susbet$cellname <- gsub("Ifit1-D", "S100a8-D", UE.decidual.susbet$cellname)
UE.decidual.susbet$cellname <- gsub("Postn-D", "Lum-D", UE.decidual.susbet$cellname)
UE.decidual.susbet$cellname <- gsub("Hist1h2ap-D", "Top2a-D", UE.decidual.susbet$cellname)
DimPlot(UE.decidual.susbet, group.by = "cellname")
UE.decidual.susbet@active.ident <- factor(UE.decidual.susbet$cellname)
UE.decidual.susbet.markers <- FindAllMarkers(UE.decidual.susbet, min.pct = 0.5, logfc.threshold = 0.5, only.pos = T)
#UE.decidual.susbet.monocle <- subset(UE.decidual.susbet, subset = time != "T105")
decidual.data <- as(as.matrix(UE.decidual.susbet@assays$RNA@counts), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", UE.decidual.susbet@meta.data)
fd <- data.frame(gene_short_name = row.names(decidual.data), row.names = row.names(decidual.data))
fd <- new("AnnotatedDataFrame", data = fd)
decidual.cds <- newCellDataSet(decidual.data,
                               phenoData = pd,
                               featureData = fd, 
                               lowerDetectionLimit = 0.5,
                               expressionFamily = negbinomial.size())
#过滤
decidual.cds <- estimateSizeFactors(decidual.cds)
decidual.cds <- estimateDispersions(decidual.cds)
decidual.cds <- detectGenes(decidual.cds, min_expr = 0.1)
expressed_genes <- rownames(subset(fData(decidual.cds), num_cells_expressed >= 10))
length(expressed_genes)
#using variable genes
#UE.decidual.susbet <- FindVariableFeatures(UE.decidual.susbet, selection.method = "vst", nfeatures = 2000)
order.genes <- differentialGeneTest(decidual.cds[expressed_genes, ],
                                    fullModelFormulaStr = "~time",
                                    cores = 20)
order.genes.top2000 <- rownames(order.genes)[order(order.genes$qval)][1:2000]
order.genes.top2000 <- unique(UE.decidual.susbet.markers$gene)
decidual.cds <- setOrderingFilter(decidual.cds, ordering_genes = order.genes.top2000)
decidual.cds <- reduceDimension(decidual.cds, method = "DDRTree")
decidual.cds <- orderCells(decidual.cds)
#decidual.cds <- orderCells(decidual.cds, root_state = 4)
decidual.cds$cellname <- factor(decidual.cds$cellname, levels = c("Lum-D", "Sfrp4-D", "Postn-D", "Top2a-D", "Hist1h2ap-D", "Gatm-D", 
                                                                  "Ptn-D", "Ifit1-D", "S100a8-D", "Prl8a2-D"))
plot_cell_trajectory(decidual.cds, color_by = "cellname", show_state_number = F, cell_size = 0.5, show_branch_points = F) + 
  scale_color_manual(values = decidual.cluster.color)
plot_cell_trajectory(decidual.cds, color_by = "Pseudotime", show_state_number = F, cell_size = 0.5, show_branch_points = F)
plot_cell_trajectory(decidual.cds, color_by = "State", show_state_number = F, cell_size = 0.5, show_branch_points = F)
save(decidual.cds, file = "20210821-decidual-pseudotime/Decidual-gene-cellname-monocle2.RData")
table(decidual.cds$cellname, decidual.cds$State)

plot_cell_trajectory(decidual.cds, color_by = "State", show_state_number = F, cell_size = 0.5, show_branch_points = F)
plot_complex_cell_trajectory(decidual.cds, color_by = "State", cell_size = 0.1)

plot_complex_cell_trajectory(decidual.cds, color_by = "cellname", cell_size = 0.1)

plot_cell_trajectory(decidual.cds, color_by = "time", cell_size = 0.1)

#随着时间的差异基因分析
diff_test_res <- differentialGeneTest(decidual.cds[expressed_genes, ], 
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)", 
                                      cores = 20)
sig_genes_names <- rownames(subset(diff_test_res, qval<0.00001 & num_cells_expressed >10 & use_for_ordering == TRUE))

decidual.pseudo.heatmap <- plot_pseudotime_heatmap(decidual.cds[sig_genes_names, ],
                        num_clusters = 8,
                        cores = 10,
                        show_rownames = T,
                        return_heatmap = T,
                        use_gene_short_name = T, 
                        hmcols = c(colorRampPalette(colors = rev(brewer.pal(11, "RdBu")))(60)))
pdf("20210827-figures/figure3/Figure3-6-pseudotime-heatmap-1960-genes.pdf", width = 10, height = 20)
decidual.pseudo.heatmap
dev.off()
decidual.pseudo.clusters <- cutree(decidual.pseudo.heatmap$tree_row, 8)
decidual.pseudo.clusters <- data.frame(genes = names(decidual.pseudo.clusters), cluster = decidual.pseudo.clusters)
decidual.pseudo.clusters$cluster <- factor(decidual.pseudo.clusters$cluster, levels = seq(1, 8, 1))
decidual.pseudo.clusters.GO <- data.frame()
for (c in levels(decidual.pseudo.clusters$cluster)) {
  print(c)
  geneset <- decidual.pseudo.clusters[decidual.pseudo.clusters$cluster == c, ]$genes
  enrich.GO <- enrich.GO.function(geneset)
  enrich.GO <- enrich.GO[1:5, ]
  tmp <- data.frame(cluster = c, 
                    GO_ID = enrich.GO$ID, 
                    GO_Term = enrich.GO$Description, 
                    pvalue = enrich.GO$pvalue, 
                    qvalue = enrich.GO$qvalue, 
                    gene = enrich.GO$geneID)
  decidual.pseudo.clusters.GO <- rbind(decidual.pseudo.clusters.GO, tmp)
}
decidual.pseudo.clusters.GO$GO_Term <- factor(decidual.pseudo.clusters.GO$GO_Term, levels = rev(unique(decidual.pseudo.clusters.GO$GO_Term)))
ggplot() + 
  geom_point(decidual.pseudo.clusters.GO, mapping = aes(x = cluster, y = GO_Term, color = cluster, size = -log10(qvalue))) + 
  scale_size(range = c(1, 8)) + 
  theme(panel.background = element_rect(color = "black", fill = "white", size = 1),
        panel.grid = element_blank(),
        axis.line = element_blank())
#URD
library(URD)
library(rgl)
load("../20210201-integrate/UE-decidual-harmony-pc30-subset-rename.RData")
UE.decidual.susbet$cellname <- gsub("Ifit1-D", "S100a8-D", UE.decidual.susbet$cellname)
UE.decidual.susbet$cellname <- gsub("Postn-D", "Lum-D", UE.decidual.susbet$cellname)
UE.decidual.susbet$cellname <- gsub("Hist1h2ap-D", "Top2a-D", UE.decidual.susbet$cellname)
#downsample up to 1000 of each cluster
UE.decidual.susbet@active.ident <- factor(UE.decidual.susbet$cellname)
UE.decidual.susbet <- subset(UE.decidual.susbet, downsample = 1000)
#UE.decidual.susbet <- subset(UE.decidual.susbet, cellname %in% c("Top2a-D", "S100a8-D", "Gatm-D"))
count.decidual <- as.matrix(UE.decidual.susbet@assays$RNA@counts)
meta.decidual <- UE.decidual.susbet@meta.data
meta.decidual <- meta.decidual[, c(1, 2, 3, 4, 10, 15, 16)]
decidual.URD <- createURD(count.data = count.decidual, meta = meta.decidual, min.cells = 3, min.counts = 3)
decidual.URD@group.ids$stage <- as.character(decidual.URD@meta[rownames(decidual.URD@group.ids), "cellname"])
decidual.URD@group.ids$time <- as.character(decidual.URD@meta[rownames(decidual.URD@group.ids), "time"])
decidual.URD <- findVariableGenes(decidual.URD, set.object.var.genes = T, diffCV.cutoff = 0.5, mean.min = 0.005, mean.max = 100, do.plot = T)

#var.genes <- sort(unique(unlist(var.by.stage)))
#decidual.URD@var.genes <- var.genes

#calculate PCA and tSNE
decidual.URD <- calcPCA(decidual.URD, mp.factor = 2)
pcSDPlot(decidual.URD)
set.seed(19)
decidual.URD <- calcTsne(decidual.URD)
plotDim(decidual.URD, "stage", plot.title = "tSNE: Stage")
plotDim(decidual.URD, "time", plot.title = "tSNE: Time")
plotDim(decidual.URD, "Sfrp4", plot.title = "tSNE: Sfrp4 expression")
#remove outliers
decidual.URD <- calcKNN(decidual.URD, nn = 120)
outliers <- knnOutliers(decidual.URD, nn.1 = 1, nn.2 = 20, x.max = 40, slope.r = 1.1, int.r = 2.9, 
                        slope.b = 0.85, int.b = 10, title = "Identifying Outliers by k-NN Distance")
#calculate diffusion map
decidual.URD.sigma <- calcDM(decidual.URD, knn = 120, sigma = 12)
#save(decidual.URD, file = "20210821-decidual-pseudotime/Decidual-URD.RData")
plotDimArray(decidual.URD.sigma, reduction.use = "dm", dims.to.plot = 1:18,  
             outer.title = "Diffusion Map (Sigma 8, 120 NNs): Stage", 
             label = "stage", plot.title = "", legend = T)
plotDim(decidual.URD, "stage", transitions.plot = 10000, plot.title = "Developmental stage (with transitions)")
root.cells <- cellsInCluster(decidual.URD.sigma, "stage", "Top2a-D")
decidual.floods.sigma <- floodPseudotime(decidual.URD.sigma, root.cells = root.cells, n = 120, minimum.cells.flooded = 0)

decidual.URD.sigma <- floodPseudotimeProcess(decidual.URD.sigma, decidual.floods.sigma, floods.name = "pseudotime")
pseudotimePlotStabilityOverall(decidual.URD.sigma)
plotDim(decidual.URD.sigma, "pseudotime")
#Investigate the distribution of pseudotime for each developmental stage.
plotDists(decidual.URD.sigma, "pseudotime", "stage", plot.title = "Pseudotime by stage (Sigma 12)")
#Find tips
tip.cells <- setdiff(cellsInCluster(decidual.URD.sigma, "stage", c("Gatm-D", "S100a8-D")), rownames(decidual.URD.sigma@pseudotime)[is.na(decidual.URD.sigma@pseudotime$pseudotime)])
decidual.URD.sigma.tip <- urdSubset(decidual.URD.sigma, cells.keep = cellsInCluster(decidual.URD.sigma, "stage", c("Gatm-D", "S100a8-D")))
decidual.URD.sigma.tip <- findVariableGenes(decidual.URD.sigma.tip, set.object.var.genes = T, diffCV.cutoff = 0.5, mean.min = 0.005, mean.max = 100, do.plot = T)
decidual.URD.sigma.tip <- calcPCA(decidual.URD.sigma.tip, mp.factor = 1.5)
pcSDPlot(decidual.URD.sigma.tip)
set.seed(20)
decidual.URD.sigma.tip <- calcTsne(decidual.URD.sigma.tip)
decidual.URD.sigma.tip <- graphClustering(decidual.URD.sigma.tip, num.nn = 50, do.jaccard = T, method = "Louvain")
plotDim(decidual.URD.sigma.tip, "Louvain-50", plot.title = "Louvain (50 NN) graph clustering", point.size=3, label.clusters = T)
#decidual.URD.sigma.tip@group.ids$`Louvain-50` <- gsub("[4,6,7]", 2, decidual.URD.sigma.tip@group.ids$`Louvain-50`)
#decidual.URD.sigma.tip@group.ids$`Louvain-50` <- gsub("[3,5,8]", 1, decidual.URD.sigma.tip@group.ids$`Louvain-50`)
#Based random walks
decidual.URD.sigma@group.ids[rownames(decidual.URD.sigma.tip@group.ids), "tip.clusters"] <- decidual.URD.sigma.tip@group.ids$`Louvain-50`
decidual.URD.sigma.ptlogistic <- pseudotimeDetermineLogistic(decidual.URD.sigma, "pseudotime", 
                                                               optimal.cells.forward = 20, max.cells.back = 40, do.plot = T)
decidual.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(decidual.URD.sigma, "pseudotime", logistic.params = decidual.URD.sigma.ptlogistic))
decidual.walks <- simulateRandomWalksFromTips(decidual.URD.sigma, tip.group.id = "tip.clusters", root.cells = root.cells, transition.matrix = decidual.biased.tm, 
                                              n.per.tip = 25000, root.visits = 1, max.steps = ncol(decidual.URD.sigma@logupx.data), verbose = T)
decidual.URD.sigma <- processRandomWalksFromTips(decidual.URD.sigma, decidual.walks, verbose = T)
plotDim(decidual.URD.sigma, "tip.clusters", plot.title = "Cells in each tip")
#Build tree
decidual.URD.sigma.tree <- loadTipCells(decidual.URD.sigma, "tip.clusters")
decidual.URD.sigma.tree <- buildTree(decidual.URD.sigma.tree, pseudotime = "pseudotime", tips.use = c(3,5,2), divergence.method = "preference", cells.per.pseudotime.bin = 25, 
                        bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh = 0.001)
save(decidual.URD.sigma.tree, file = "20210827-figures/figure4/Top2a-Gatm-S100a8-URD-tree.RData")
pdf("20210827-figures/figure4/Top2a-Gatm-S100a8-URD-tree.pdf", width = 5, height = 4)
plotTree(decidual.URD.sigma.tree, "stage", title = "Developmental Stage", cell.size = 1, discrete.colors = c("#08519C", "#FDB462", "#FB8072"))
dev.off()
#Force-directed layout
decidual.URD.sigma.tree <- treeForceDirectedLayout(decidual.URD.sigma.tree, num.nn = 100, cut.unconnected.segments = 2, verbose = T)
plotTreeForce(decidual.URD.sigma.tree, "stage", title = "stage", title.cex = 2, title.line = 2.5)
plotTreeForceStore3DView(decidual.URD.sigma.tree, "stage")
save(decidual.URD.sigma.tree, file = "20210821-decidual-pseudotime/Decidual-URD-tree.RData")

#只计算AM区域的pseudotime
library(monocle)
load("../20210201-integrate/UE-decidual-harmony-pc30-rename.RData")
DimPlot(UE.decidual)
UE.decidual.susbet <- subset(UE.decidual, subset = cellname != "cluster9-1" & cellname != "Prg-" )
UE.decidual.susbet$cellname <- gsub("Ifit1-D", "S100a8-D", UE.decidual.susbet$cellname)
UE.decidual.susbet$cellname <- gsub("Hist1h2ap-D", "Top2a-D", UE.decidual.susbet$cellname)
UE.decidual.susbet$cellname <- gsub("Postn-D", "Lum-D", UE.decidual.susbet$cellname)
DimPlot(UE.decidual.susbet, group.by = "cellname")
UE.decidual.susbet@active.ident <- factor(UE.decidual.susbet$cellname)
UE.decidual.AM <- subset(UE.decidual.susbet, subset = cellname %in% c("Lum-D", "Sfrp4-D", "Gatm-D", "Prl8a2-D"))
UE.decidual.MET <- subset(UE.decidual.susbet, subset = cellname %in% c("Lum-D", "Sfrp4-D", "Top2a-D", "S100a8-D", "Gatm-D"))
save(UE.decidual.AM, file = "../20210201-integrate/UE-decidual-AM-20210906.RData")
save(UE.decidual.MET, file = "../20210201-integrate/UE-decidual-MET-20210906.RData")
decidual.data <- as(as.matrix(UE.decidual.AM@assays$RNA@counts), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", UE.decidual.AM@meta.data)
fd <- data.frame(gene_short_name = row.names(decidual.data), row.names = row.names(decidual.data))
fd <- new("AnnotatedDataFrame", data = fd)
decidual.cds <- newCellDataSet(decidual.data,
                               phenoData = pd,
                               featureData = fd, 
                               lowerDetectionLimit = 0.5,
                               expressionFamily = negbinomial.size())
#过滤
decidual.cds <- estimateSizeFactors(decidual.cds)
decidual.cds <- estimateDispersions(decidual.cds)
decidual.cds <- detectGenes(decidual.cds, min_expr = 0.1)
expressed_genes <- rownames(subset(fData(decidual.cds), num_cells_expressed >= 10))
length(expressed_genes)
AM.DEGs <- FindAllMarkers(UE.decidual.AM, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
#using variable genes
#UE.decidual.susbet <- FindVariableFeatures(UE.decidual.susbet, selection.method = "vst", nfeatures = 2000)
#order.genes <- differentialGeneTest(decidual.cds[expressed_genes, ],
#                                    fullModelFormulaStr = "~time",
#                                    cores = 20)
#order.genes.top2000 <- rownames(order.genes)[order(order.genes$qval)][1:2000]
decidual.cds <- setOrderingFilter(decidual.cds, ordering_genes = unique(AM.DEGs$gene))
decidual.cds <- reduceDimension(decidual.cds, method = "DDRTree")
decidual.cds <- orderCells(decidual.cds)

plot_cell_trajectory(decidual.cds, color_by = "time", cell_size = 0.1)
plot_cell_trajectory(decidual.cds, color_by = "cellname", cell_size = 0.1)

#利用空间数据做pseudotime monocle3
#merge所有的bin50空间数据
Decidual.bin50 <- merge(E65.2.bin50, list(E75.2.bin50, B6.bin50, E95.bin50), add.cell.ids = c("E65.bin50", "E75.bin50", "E85.bin50", "E95.bin50"))
Decidual.bin50@active.ident <- factor(Decidual.bin50$celltype)
Decidual.bin50.subsample <- subset(Decidual.bin50, downsample = 1000)
Decidual.bin50.subsample <- subset(Decidual.bin50.subsample, subset = celltype %in% c("Lum-D", "Sfrp4-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "Top2a-D"))
save(Decidual.bin50.subsample, file = "20210821-decidual-pseudotime/Decidual-ST-downsample-1000-seurat.RData")

#subset MET的annotation的细胞来做
MET.cells <- c(paste("E65.bin50", E65.2.bin50.celltype.df.annotation.MET$cell_name, sep = "_"), 
               paste("E75.bin50", E75.2.bin50.celltype.df.annotation.MET$cell_name, sep = "_"), 
               paste("E85.bin50", E85.B6.bin50.celltype.df.annotation.MET$cell_name, sep = "_"), 
               paste("E95.bin50", E95.bin50.celltype.df.annotation.MET$cell_name, sep = "_"))
Decidual.bin50.subsample <- subset(Decidual.bin50, cells = MET.cells)
Decidual.bin50.subsample <- subset(Decidual.bin50.subsample, subset = celltype %in% c("Lum-D", "Sfrp4-D", "Gatm-D", "S100a8-D", "Prl8a2-D", "Top2a-D"))
save(Decidual.bin50.subsample, file = "20210821-decidual-pseudotime/Decidual-ST-downsample-1000-seurat.RData")

