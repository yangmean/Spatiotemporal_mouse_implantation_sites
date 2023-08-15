#Author: Min Yang
#Using spatial data to construct the pseudotime analyses
#Merge all the bin50 spatial data and downsample 1000 cells for each cell type, saved as DSC.bin50.ST.object
library("monocle3")
library("ggplot2")
library("dplyr")
DSC.data <- as(as.matrix(DSC.bin50.ST.object@assays$Spatial@counts), 'sparseMatrix')
pd <- DSC.bin50.ST.object@meta.data
fd <- data.frame(gene_short_name = row.names(DSC.data), row.names = row.names(DSC.data))
DSC.cds <- new_cell_data_set(DSC.data,
                             cell_metadata = pd,
                             gene_metadata = fd)
DSC.cds <- preprocess_cds(DSC.cds, num_dim = 50)
DSC.cds <- align_cds(DSC.cds, alignment_group = "orig.ident", residual_model_formula_str = "~ nCount_Spatial")
DSC.cds <- reduce_dimension(DSC.cds)
plot_cells(DSC.cds, label_groups_by_cluster = FALSE, color_cells_by = "celltype")

DSC.cds <- cluster_cells(DSC.cds)
plot_cells(DSC.cds, color_cells_by = "cluster", label_cell_groups = T, graph_label_size = 6)
DSC.cds <- learn_graph(DSC.cds)
#We choose the cells in eSF/pre-DSC as the root cells
save(DSC.cds, file = "DSC.cds.monocle3.RData")

DSC.cds.psdeudotime <- data.frame(DSC.cds@reducedDims$UMAP,
                                  cellname = colnames(DSC.cds), 
                                  pseudotime = DSC.cds@principal_graph_aux$UMAP$pseudotime, 
                                  celltype = DSC.cds@colData$celltype, 
                                  time = DSC.cds@colData$orig.ident)
colnames(DSC.cds.psdeudotime)[1:2] <- c("UMAP_1", "UMAP_2")
write.table(DSC.cds.psdeudotime, file = "DSC.cds.psdeudotime.dataframe.csv", sep = ",")

#load DSC.cdsobject and DSC.cds.pseudotime data.frame to monocle2 to run DEGs align the trajectory
DSC.cds.psdeudotime <- DSC.cds.psdeudotime[!is.infinite(DSC.cds.psdeudotime$pseudotime), ]
DSC.cds.psdeudotime.value <- DSC.cds.psdeudotime$pseudotime
names(DSC.cds.psdeudotime.value) <- rownames(DSC.cds.psdeudotime)
DSC.cds.subsample <- subset(DSC.cds, cells = rownames(DSC.cds.psdeudotime))
DSC.cds.subsample@active.ident <- factor(DSC.cds.subsample$celltype)
DSC.cds.subsample <- AddMetaData(DSC.cds.subsample, metadata = DSC.cds.psdeudotime.value, col.name = "pseudotime")
DSC.subsample.data <- as(as.matrix(DSC.cds.subsample@assays$SCT@data), 'sparseMatrix')
pd <- new("AnnotatedDataFrame", DSC.cds.subsample@meta.data)
fd <- data.frame(gene_short_name = row.names(DSC.subsample.data), row.names = row.names(DSC.subsample.data))
fd <- new("AnnotatedDataFrame", data = fd)
DSC.monocle2.cds <- newCellDataSet(DSC.subsample.data,
                                   phenoData = pd,
                                   featureData = fd, 
                                   lowerDetectionLimit = 0.5,
                                   expressionFamily = negbinomial.size())
DSC.monocle2.cds <- estimateSizeFactors(DSC.monocle2.cds)
DSC.monocle2.cds <- estimateDispersions(DSC.monocle2.cds)
DSC.monocle2.cds <- detectGenes(DSC.monocle2.cds, min_expr = 0.1)

#Find all the DEGs among the DSC subclusters, including D1, D2, D3, D5, D6, D7. We excluded D4 due to limit numbers in spatial data
DSC.MET <- subset(DSC.seurat.object, subset = subname %in% c("D1_eSF", "D2_pre-DSC", "D3_proliferating DSC", "D5_Nourishing DSC", "D6_Angiogenic DSC", "D7_Postmature DSC"))
DSC.MET.markers <- FindAllMarkers(DSC.MET, min.pct = 0.25, logfc.threshold = 0.5, only.pos = F)
DSC.MET.DEGs <- intersect(rownames(DSC.monocle2.cds), unique(DSC.MET.markers$gene))
set.seed(10000)
DSC.MET.DEGs.heatmap <- plot_pseudotime_heatmap(DSC.monocle2.cds[DSC.MET.DEGs, ],
                                                num_clusters = 7,
                                                cores = 10,
                                                show_rownames = T,
                                                return_heatmap = T,
                                                use_gene_short_name = T, 
                                                hmcols = c(colorRampPalette(colors = rev(brewer.pal(11, "RdBu")))(60)))


DSC.MET.DEGs.heatmap.cluster <- cutree(DSC.MET.DEGs.heatmap$tree_row, 7)
DSC.MET.DEGs.heatmap.cluster.df <- data.frame(gene = names(DSC.MET.DEGs.heatmap.cluster), cluster = DSC.MET.DEGs.heatmap.cluster)
DSC.MET.DEGs.heatmap.cluster.df$cluster <- factor(DSC.MET.DEGs.heatmap.cluster.df$cluster, levels = c(2,1,6,7,5,3,4))
write.table(DSC.MET.DEGs.heatmap.cluster.df, file = "DSC-ST-pseudotime-heatmap-920genes-7clusters-gene.csv", sep = ",", col.names = T, row.names = F, quote = F)










